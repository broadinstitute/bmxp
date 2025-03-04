"""
blueshift is a LCMS dataset drift correction package, designed to correct for drift or
sudden shifts in signal intensity, as measured in pooled reference samples or internal
standards.
"""

from warnings import warn
import re
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import interpolate
from bmxp import FMDATA, IMDATA, POOL_INJ_TYPES

__version__ = "0.2.1"


class DriftCorrection:  # pylint: disable=too-many-instance-attributes
    """
    DriftCorrection stores a dataset and corresponding sample information sheet and
    corrects sample signals using pooled reference sample or internal standard-based
    correction methods.
    """

    def __init__(
        self,
        raw_data,
        sample_information,
        default_met_col=None,
        schema_labels=None,
        pool_inj_types=None,
    ):
        """
        Initializes the DriftCorrection object with a dataset and corresponding sample
        information sheet.
        :param raw_data: DataFrame or str, dataset or filepath to .csv dataset
        :param sample_information: DataFrame or str, sample information or filepath to
        .csv sample information
        :param pool_names: tuple[str], pooled reference sample names
        """
        fmdata = FMDATA.copy()
        imdata = IMDATA.copy()
        if schema_labels is not None:
            fmdata.update(schema_labels)
            imdata.update(schema_labels)
        if pool_inj_types is None:
            pool_inj_types = POOL_INJ_TYPES

        self.data = pd.DataFrame()
        self.rdata = pd.DataFrame()
        self.metadata = pd.DataFrame()
        self.sample_info = pd.DataFrame()
        self.pools = {}
        self.batches = {}
        self.correction_params = {}

        # column labels
        self.injection_id = imdata["injection_id"]
        self.injection_type = imdata["injection_type"]
        self.column_number = imdata["column_number"]
        self.injection_order = imdata["injection_order"]
        self.batches_label = imdata["batches"]
        self.QCRole = imdata["QCRole"]
        self.Non_Quant = fmdata["Non_Quant"]
        self.Batches_Skipped = fmdata["Batches Skipped"]

        if default_met_col:
            self.met_col = default_met_col
        else:
            self.met_col = fmdata["Metabolite"]

        if not isinstance(raw_data, pd.DataFrame):
            raw_data = pd.read_csv(raw_data)

        if not isinstance(sample_information, pd.DataFrame):
            sample_information = pd.read_csv(sample_information)

        raw_data, sample_information, self.first_sample = self._validate_data(
            raw_data, sample_information
        )

        self.data, self.metadata = self._parse_data(raw_data)
        self.rdata = self.data.copy()
        self.cvs = pd.DataFrame(index=self.rdata.index)
        self.sample_info = sample_information
        self.pools, self.not_used_pools = self._find_pools(pool_inj_types)
        self.batches["default"], self.batches["override"] = self._gen_batches()

    def _validate_data(self, raw_data, sample_information):
        """
        Checks that the dataset and sample information have valid headers and contents
        :param raw_data: DataFrame, raw dataset
        :param sample_information: DataFrame or str, sample information or filepath to
        .csv sample information
        """
        # required columns in sample information
        required_columns = [
            self.injection_id,
            self.injection_type,
            self.column_number,
            self.injection_order,
            self.batches_label,
        ]
        for col in required_columns:
            if col not in sample_information.columns:
                raise ValueError(f"Column '{col}' not found in sample information.")

        # drop fully blank rows in injection info or columns without names in data
        sample_information = sample_information.dropna(axis="index", how="all")
        has_col_name = ~raw_data.columns.str.startswith("Unnamed: ")
        raw_data = raw_data.loc[:, has_col_name]
        sample_information[self.batches_label] = sample_information[
            self.batches_label
        ].fillna("")
        sample_information = sample_information.astype(
            {self.injection_type: str, self.batches_label: str}
        )

        try:
            sample_information = sample_information.astype({self.injection_order: int})
        except ValueError as e:
            raise TypeError(
                "There are missing or non-numeric values in the injection_order."
            ) from e

        if not sample_information[self.injection_order].is_unique:
            raise ValueError("There are duplicate values in the injection_order.")

        if not sample_information[self.injection_id].is_unique:
            raise ValueError("There are duplicate injection_ids.")
        sample_information = sample_information.set_index(self.injection_id, drop=False)

        if not (
            sample_information.sort_values(by=self.injection_order).index
            == sample_information.index
        ).all():
            raise ValueError("Samples must be sorted in the order they were injected.")

        # check all data samples are in info, and get first sample
        in_data = np.isin(sample_information[self.injection_id], raw_data.columns)
        if not (in_data).any():
            raise ValueError(
                "There are no overlapping samples in your sample information and "
                "sample data sheets."
            )
        first_sample = sample_information.loc[in_data, self.injection_id].iloc[0]
        samples = raw_data.loc[:, first_sample:].columns

        # check that batch ends aren't on dropped samples
        non_dc = sample_information[self.injection_type].str.startswith(
            tuple(["not_used", "brpp", "mm", "blank", "other", "tube_blank", "ms2"])
        )
        will_drop = ~in_data & non_dc
        end_drops = sample_information.loc[will_drop, self.batches_label] == "batch end"
        if sum(end_drops) > 0:
            raise ValueError(
                "The 'batch end' label cannot be applied to injections that are "
                "not present in the dataset. Please move the following batch ends "
                "to valid injections: "
                + ", ".join(
                    sample_information.loc[
                        end_drops.index[end_drops], self.injection_id
                    ].values
                )
            )

        # flag unused samples
        sample_information.loc[will_drop, self.QCRole] = "NA"

        if len(samples) > len(sample_information.loc[~will_drop]):
            raise ValueError(
                "The number of samples in sample information does not match the number "
                "of data columns. Check for excess samples at the end of the data."
            )

        # check there aren't 'used' samples missing from data sheet
        missing = ~np.isin(
            sample_information.loc[~will_drop, self.injection_id], samples
        )
        if missing.any():
            raise ValueError(
                "The following samples are marked as usable but missing from your "
                "data sheet: "
                + ", ".join(
                    sample_information.loc[~will_drop, self.injection_id][missing]
                )
            )

        # check they are in the same order
        same = sample_information.loc[~will_drop, self.injection_id] == samples
        if not same.all():
            bad_sample = sample_information.loc[
                ~will_drop & ~same, self.injection_id
            ].iloc[0]
            raise ValueError(
                f"Your usable samples do not match, starting with '{bad_sample}' "
                f"on your information sheet."
            )

        return raw_data, sample_information, first_sample

    def _parse_data(self, raw_data):
        """
        Separates non-data columns from sample data columns, and cleans the data as
        needed before calculations are performed.
        :param raw_data: DataFrame, raw dataset
        :return:
        """

        raw_data = raw_data.reset_index(drop=True)
        first_sample_idx = raw_data.columns.get_loc(self.first_sample)
        metadata = raw_data.iloc[:, :first_sample_idx]
        data = raw_data.iloc[:, first_sample_idx:]

        try:
            data = data.astype("float64")
        except ValueError as e:
            character = re.search("'.+'", str(e))
            if character:
                raise TypeError(
                    f"Data contains a non-numeric character: " f"{character.group()}."
                ) from None
            raise TypeError("Data contains a non-numeric character.") from e
        data = data.replace(0, np.nan)

        if self.met_col in metadata.columns:
            metadata = metadata.astype({self.met_col: str})
        metadata = metadata.replace("nan", np.nan)

        return data, metadata

    def _find_pools(self, pool_names):
        """
        Locates file names of pooled reference samples.
        :param pool_names: tuple[str], pooled reference sample names
        :return: two params
            dict: pool names as keys, list of file names of pools as values
            pd series: list of not-used pools
        """
        pools = {}
        a, b = pool_names
        info = self.sample_info.sort_values(by=self.injection_order)

        # identify pools based on sample type
        # non-case sensitive, for the time being
        pools[a] = info.loc[
            info[self.injection_type].str.lower() == a.lower(), self.injection_id
        ]
        pools[b] = info.loc[
            info[self.injection_type].str.lower() == b.lower(), self.injection_id
        ]

        # not used -- if it starts with "not_used" and contains a pool name
        not_used_pools = info.loc[
            info[self.injection_type].str.contains(f"{a}|{b}", case=False)
            & info[self.injection_type].str.startswith("not_used"),
            self.injection_id,
        ]

        return pools, not_used_pools

    def _gen_batches(self):
        """
        Splits data into default batches based on column breaks and override batches
        based on user-supplied batch start/end.
        :return:
            default_batches - list, 2D nested containing sample file names for each
            batch
            override_batches - list, 2D nested containing sample file names for each
            batch
        """

        default_batches = []
        override_batches = []

        # accept "batch_end" and "batch end" as equivalent
        self.sample_info[self.batches_label] = self.sample_info[
            self.batches_label
        ].str.replace("_", " ")
        for item in set(self.sample_info[self.batches_label]):
            if item.lower() not in ("batch start", "batch end", ""):
                raise ValueError(
                    f"An invalid label is present in the batches column: '{item}'. "
                    f"Valid batch labels are 'batch start' and 'batch end'."
                )

        valid = self.sample_info[self.QCRole] != "NA"

        # default batches split by column number
        for col in pd.unique(self.sample_info[self.column_number]):
            batch = self.sample_info.loc[
                valid & (self.sample_info[self.column_number] == col), self.injection_id
            ]
            if len(batch) > 0:
                default_batches.append(batch)
                # label column ends as batch ends, whether the user has done so or not
                self.sample_info.loc[batch.index[-1], self.batches_label] = "batch end"

        # label very first injection as batch start, whether the user has done so or not
        self.sample_info.loc[self.sample_info.index[0], self.batches_label] = (
            "batch start"
        )

        # override batches split by user-indicated batch end
        sample_info = self.sample_info[valid]
        batch_info = sample_info[self.batches_label].str.lower()
        ends = batch_info == "batch end"
        starts = np.roll(ends, 1)

        for start, end in zip(np.where(starts)[0], np.where(ends)[0]):
            this_batch = sample_info.iloc[start : end + 1, :]
            override_batches.append(this_batch[self.injection_id])

            if (
                sample_info.loc[sample_info.index[start], self.batches_label]
                != "batch start"
            ):
                warn(
                    "The following sample is being used as a batch start but lacks a "
                    "'batch start' label: "
                    + sample_info.loc[sample_info.index[start], self.injection_id]
                )

            contains_starts = np.where(this_batch[self.batches_label] == "batch start")[
                0
            ]
            if (contains_starts > 0).any():
                warn(
                    "The following samples are NOT being used as a batch start but "
                    "have a 'batch start' label: "
                    + ", ".join(
                        [
                            this_batch.loc[this_batch.index[i], self.injection_id]
                            for i in contains_starts
                            if i != 0
                        ]
                    )
                )

        return default_batches, override_batches

    def internal_standard_correct(self, internal_standards):
        """
        Performs drift correction using one or more internal standards as a reference.
        :param internal_standards: str or array-like, single internal standard or
        list of internal standards to use for drift correction
        """
        if isinstance(internal_standards, str):
            internal_standards = [internal_standards]
        if self.met_col not in self.metadata.columns:
            raise ValueError(
                f"'{self.met_col}' was not found in your columns. Cannot find internal "
                "standard(s) for drift correction."
            )

        is_nonquant = (
            np.array([False] * len(self.metadata))
            if self.Non_Quant not in self.metadata.columns
            else self.metadata[self.Non_Quant].fillna(False)
        )
        if not set(internal_standards).issubset(
            set(self.metadata.loc[~is_nonquant, self.met_col])
        ):
            missing = set(internal_standards) - set(
                self.metadata.loc[~is_nonquant, self.met_col]
            )
            raise ValueError(
                "Internal standard(s) " + str(missing)[1:-1] + " not found in "
                "Metabolite list."
            )

        standards_index = (
            np.isin(self.metadata[self.met_col], internal_standards) & ~is_nonquant
        )
        scalers = self.data.loc[standards_index]
        scalers = scalers.transpose() / scalers.median(axis=1)
        scalers = scalers.mean(axis=1)

        self.cvs.loc[standards_index, "Internal Standards"] = "Internal Standard"

        if scalers.isna().any():
            raise ValueError(
                "The following samples cannot be corrected due to missing values in "
                f"the internal standard(s): {self.data.columns[scalers.isna()].values}"
            )

        self.data = self.data / scalers

        # store drift correction parameters for reporting later
        self.correction_params["IS"] = list(internal_standards)

    def pool_correct(
        self, interpolation, pool, override=True, max_missing_percent=30
    ):  # pylint: disable=too-many-branches
        """
        Performs pool-based drift correction.
        :param interpolation: str, interpolation method: NN or linear
        :param pool: str, pool name to use for correction (typically PREFA or PREFB)
        :param override: boolean, uses override batches in place of default batches
        if true
        :param max_missing_percent: int (0-100), skip drift correction for any row
        that has this percent (or more) drift correction pools missing values
        """

        if pool not in self.pools:
            raise ValueError(
                f"'{pool}' is not a valid pool for this dataset. Valid pools: "
                f"{', '.join(self.pools)}."
            )

        methods = {
            "nn": self._nn_interpolation,
            "linear": self._linear_interpolation,
        }
        self.correction_params[interpolation.lower()] = (pool, override)
        interpolation = methods[interpolation.lower()]
        pool_names = self.pools[pool]
        pool_medians = self.data[pool_names].apply(
            lambda row: 1 if row.isnull().all() else row.median(), axis=1
        )
        batches = self.batches["override"] if override else self.batches["default"]
        self._set_pool_roles(pool)

        # begin the scaling
        scalers = pd.DataFrame(index=self.data.index)
        data = self.data.values
        self.cvs[self.Batches_Skipped] = [[] for _ in self.cvs.index]

        self._flag_skipped_rows(pool, max_missing_percent, len(batches))

        for this_batch in batches:
            batch_pools = this_batch.loc[this_batch.isin(pool_names)]
            batch_pools = pd.DataFrame(batch_pools).join(
                self.sample_info[self.injection_order]
            )
            this_batch = pd.DataFrame(this_batch).join(
                self.sample_info[self.injection_order]
            )

            if len(batch_pools) == 0:
                raise ValueError(
                    "The following batch does not contain any drift "
                    "correction pools: "
                    f"{this_batch.iloc[0, 0]} to "
                    f"{this_batch.iloc[-1, 0]}"
                )

            batch_scalers = []
            for row in tqdm(range(len(data))):
                pool_values = data[row, np.isin(self.data.columns, batch_pools.index)]
                pool_order = batch_pools[self.injection_order].values[pool_values > 0]
                pool_values = pool_values[pool_values > 0]

                # skip if previously flagged
                if self.cvs.loc[row, self.Batches_Skipped] == "Not Pool Corrected":
                    batch_scalers.append([pool_medians[row]] * len(this_batch))

                # skip batch if there are no pools
                elif len(pool_order) < 1:
                    batch_scalers.append([pool_medians[row]] * len(this_batch))
                    self.cvs.loc[row, self.Batches_Skipped].append(
                        f"{this_batch.iloc[0,0]} to {this_batch.iloc[-1,0]}"
                    )

                # "nearest neighbor" if there is 1 pool and this isn't the only batch
                elif len(pool_order) < 2:
                    batch_scalers.append([pool_values[0]] * len(this_batch))

                else:
                    batch_scalers.append(
                        interpolation(
                            this_batch, pool_order, pool_values, self.injection_order
                        )
                    )

            scalers = scalers.join(
                pd.DataFrame(batch_scalers, columns=this_batch[self.injection_id])
            )
        scalers = scalers.divide(pool_medians, axis="index")

        self.data = self.data.divide(scalers)

    def _set_pool_roles(self, pool):
        dc_pools = set(self.pools[pool].index)
        qc_pools = set()
        for pool_list in self.pools.values():
            qc_pools.update(pool_list.index)
        qc_pools -= dc_pools
        self.sample_info.loc[self.sample_info[self.QCRole] != "NA", self.QCRole] = (
            "Sample"
        )
        self.sample_info.loc[list(dc_pools), self.QCRole] = "QC-drift_correction"
        self.sample_info.loc[list(qc_pools), self.QCRole] = "QC-pooled_ref"
        self.sample_info.loc[self.not_used_pools.index, self.QCRole] = "QC-not_used"

    def _flag_skipped_rows(self, pool, max_missing_percent, batch_count):
        pools = np.isin(self.data.columns, self.pools[pool])
        for i, row in enumerate(self.data.values):
            pool_count = sum(row[pools] > 0)
            if (
                pool_count == 0
                or pool_count / len(self.pools[pool])
                < (100 - max_missing_percent) / 100
                or (batch_count < 2 and pool_count < 2)
            ):
                self.cvs.loc[i, self.Batches_Skipped] = "Not Pool Corrected"

    @staticmethod
    def _nn_interpolation(
        batch, pool_order, pool_values, order_label="injection_order"
    ):
        """
        :param batch: DataFrame, sample file names and injection order for the current
        batch
        :param pool_order: array-like, injection order for each pool in the batch,
        for a single metabolite
        :param pool_values: array-like, pool value for each pool in the batch,
        for a single metabolite
        :return: numpy ndarray, value of nearest neighbor for each sample in the batch
        """

        nn_x = batch[order_label].values
        nn_interp = interpolate.interp1d(
            pool_order,
            pool_values,
            kind="nearest",
            bounds_error=False,
            fill_value=(pool_values[0], pool_values[-1]),
        )
        nn_y = nn_interp(nn_x)

        return nn_y

    @staticmethod
    def _linear_interpolation(
        batch, pool_order, pool_values, order_label="injection_order"
    ):
        """
        :param batch: DataFrame, sample file names and injection order for the current
        batch
        :param pool_order: array-like, injection order for each pool in the batch,
        for a single metabolite
        :param pool_values: array-like, pool value for each pool in the batch,
        for a single metabolite
        :return: numpy ndarray, linearly interpolated value for each sample in the batch
        """

        linear_x = batch[order_label].values
        linear_y = np.interp(x=linear_x, xp=pool_order, fp=pool_values)

        return linear_y

    def calculate_cvs(self):
        """
        Calculates coefficients of variation for each row of data.
        """
        pools = list(self.pools.keys())
        all_data = {  # explicitly define key order, this will be the column order later
            "Raw": self.rdata,
            f"Raw {pools[0]}": self.rdata[self.pools[pools[0]]],
            f"Raw {pools[1]}": self.rdata[self.pools[pools[1]]],
            "Final": self.data,
            pools[0]: self.data[self.pools[pools[0]]],
            pools[1]: self.data[self.pools[pools[1]]],
        }

        for name, data in all_data.items():
            col_name = f"{name} CVs" if name != "Final" else "CVs"
            self.cvs[col_name] = data.std(axis=1) / data.mean(axis=1)

    def to_csv(self, filepath="results.csv", to_string=False):
        """
        Assembles and downloads a csv of the metadata, data, and calculated CVs
        :param filepath: base/beginning of the results filepath
        """
        filepath = filepath.replace(".csv", "")
        for method, params in self.correction_params.items():
            if True in params:
                filepath += f"_{method}_{params[0]}_override"
            elif False in params:
                filepath += f"_{method}_{params[0]}"
            else:
                filepath += f"_{method}"

        # sort CV columns so raw columns come first
        headers = [i for i in self.cvs.columns if "CVs" in i]
        headers.sort(key=lambda x: x.replace("Raw", "AAA"), reverse=True)
        for col in headers:
            self.cvs.insert(0, col, self.cvs.pop(col))

        # make skipped batches lists more readable
        if self.Batches_Skipped in self.cvs.columns:
            self.cvs[self.Batches_Skipped] = self.cvs[self.Batches_Skipped].map(
                lambda x: ", ".join(x) if isinstance(x, list) else x
            )

        output = self.data.apply(np.floor)
        output = output.replace(0, np.nan)
        output = output.astype("Int64")

        output = self.metadata.loc[:, self.metadata.columns != self.met_col].join(
            [self.cvs, self.metadata.loc[:, self.met_col], output]
        )
        if to_string:
            return (
                output.to_csv(lineterminator="\r\n", index=False),
                filepath,
            )
        return output.to_csv(filepath + ".csv", index=False)

    def report(self, filepath="report.csv", to_string=False):
        """
        Writes a report
        :param filepath: report filepath name
        """
        filepath = filepath.replace(".csv", "") + ".csv"
        if to_string:
            return (
                self.sample_info.to_csv(lineterminator="\r\n", index=False),
                filepath,
            )
        return self.sample_info.to_csv(filepath, index=False)
