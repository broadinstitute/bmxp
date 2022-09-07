"""
MSEclipse is a LCMS dataset alignment package, designed to combine datasets which have
been collected on the same method but perhaps on different instruments or columns.
"""
# pylint: disable=too-many-lines
import copy
import logging
import collections
import io
import networkx as nx
import numpy as np
import pandas as pd
from scipy import interpolate
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

np.random.seed(0)

lowess = sm.nonparametric.lowess
__version__ = "0.0.2"


def dataset_loops(attr=None):
    """
    Wraps a function with a double For-loop for iterating over the datasets.
    It converts the supplied datasets, which may be None, 1 dataset, or multiple,
    to a list.
    Also handles creation of a dict key for an attribute and makes it if it doesn't
     exist. E.g. self.scalers[DS1] -- if that key doesn't exist, it will create it
     so that self.scalers[DS1][DS2] can be assigned.
    """

    def attr_handler(func):
        def loops(self, datasets_1=None, datasets_2=None, *args, **kwargs):
            # convert datasets to a list of datasets for the loop
            if datasets_1 is None:
                datasets_1 = list(self.datasets)
            if isinstance(datasets_1, str):
                datasets_1 = [datasets_1]

            if datasets_2 is None:
                datasets_2 = list(self.datasets)
            if isinstance(datasets_2, str):
                datasets_2 = [datasets_2]

            for ds1 in datasets_1:
                if attr and ds1 not in self.__getattribute__(attr):
                    self.__getattribute__(attr)[ds1] = {}
                for ds2 in datasets_2:
                    if ds1 == ds2:
                        continue
                    func(self, ds1, ds2, *args, **kwargs)

        return loops

    return attr_handler


def ds2_loop(func):
    """
    Wraps a function with an inner ds2 loop. Does not check for attribute
     dictionary keys (handled in functions)
    """

    def wrapper(self, ds1=None, ds2=None, *args, **kwargs):
        ds2 = list(self.datasets) if ds2 is None else ds2
        if not isinstance(ds2, str):
            for item in ds2:
                if item == ds1:
                    continue
                func(self, ds1, item, *args, **kwargs)
        else:
            func(self, ds1, ds2, *args, **kwargs)

    return wrapper


class MSAligner:  # pylint: disable=too-many-instance-attributes
    """
    MSAligner stores 2 or more MS Datasets and performs matching, scaling, and
    alignments between them.
    """

    def __init__(self, *args, names=None):
        """
        Initializes the MSAligner object with a collection of datasets and their names.

        :param args: tuple[str], filepaths of .csv datasets
        :param names: list[str], names of datasets
        """
        self.descriptors = {"RT": "linear", "MZ": "ppm", "Intensity": "log10"}
        self.default_cutoffs = {"RT": 6, "Intensity": 6, "MZ": 6}
        self.default_weights = {"RT": 1, "Intensity": 1, "MZ": 1}
        self.default_coarse_params = {
            "RT": {"upper": 0.5, "lower": -0.5},
            "MZ": {"upper": 15, "lower": -15},
            "Intensity": {"upper": 2, "lower": -2},
        }
        self.feature_name = "Compound_ID"
        self.anchor_priority = "Intensity"
        self.datasets = collections.OrderedDict()
        self.anchors = {}
        self.coarse_matches = {}
        self.scalers = {}
        self.scaled_values = {}
        self.stds = {}
        self.matches = {}
        self.cutoffs = {}
        self.weights = {}
        self.coarse_params = {}
        self.multipliers = {}
        self.graph = nx.DiGraph()
        self.remove_all = True
        # convert names to a list if not supplied
        if names is None:
            names = [None] * len(args)

        # Check there are the same number of datasets as names
        if len(names) != len(args):
            raise ValueError(
                "The number of names you provided does not match the "
                "number of datasets."
            )

        # Add datasets using the add_dataset method
        for arg, name in zip(args, names):
            self.add_dataset(arg, name)

    def add_dataset(self, dataset, name=None):
        """
        Adds a dataset to the existing alignment object.

        :param dataset: str, filepath of .csv dataset
        :param name: str, name of dataset
        """

        if name is None:
            # check if there are any dataframes present; if so, raise an exception
            if isinstance(dataset, pd.DataFrame):
                raise ValueError(
                    "Names must be provided if importing DataFrames directly."
                )
            name = dataset

        if not isinstance(dataset, pd.DataFrame):
            dataset = pd.read_csv(dataset)

        # check that the name isn't already present
        if name in self.datasets.keys():
            raise ValueError("You have duplicate names for your datasets.")

        self.datasets[name] = dataset
        self.datasets[name].index = self.datasets[name].index.map(str)
        self._validate_data(name)

    def _validate_data(self, names=None):
        """
        Performs a check to see if a dataset has any compatibility errors.

        :param names: list, name of all datasets to check
        """

        # If names are not provided, check all
        if names is None:
            names = self.datasets.keys()

        # If one name is provided as a string, convert to a 1 member list
        elif isinstance(names, str):
            names = [names]

        # Check all datasets included in the names
        for name in names:
            ds = self.datasets[name]
            # calculate intensity if needed
            if "Intensity" not in ds.columns and "Intensity" in self.descriptors:
                gen_intensity(
                    ds,
                    ignore=[descr for descr in self.descriptors if descr != "Intensity"]
                    + [self.feature_name],
                )

            # Check required columns
            required_columns = [*self.descriptors, self.feature_name]
            for req_col in required_columns:
                if req_col not in ds.columns:
                    raise ValueError(f"Column {req_col} not found in dataset {name}.")

    def set_defaults(self, params):
        """
        Updates default values for weights and cutoffs.
        """

        for key, item in params.items():
            if key.lower() == "cutoffs":
                self.default_cutoffs.update(item)
            elif key.lower() == "weights":
                self.default_weights.update(item)
            elif key.lower() == "coarse_params":
                self.default_coarse_params.update(item)
            else:
                raise KeyError(
                    f"Unknown value: {key}. Acceptable values are 'cutoffs', 'weights',"
                    " and 'coarse_params'."
                )

    def set_params(self, params, rec=True):
        """
        Sets weights, cutoffs, and scalers for any dataset pairs
        For example:
        {"scalers": {
            "DS1": {"DS2": {"RT": 1}}
            },
        "cutoffs": { ...}
        }
        """
        for key, item in params.items():
            if key.lower() == "scalers":
                self.set_scalers(item, rec)
            elif key.lower() == "cutoffs":
                self.set_cutoffs(item, rec)
            elif key.lower() == "weights":
                self.set_weights(item, rec)
            elif key.lower() == "coarse_params":
                self.set_coarse_params(item, rec)
            else:
                raise KeyError(
                    f"Unknown value: {key}. Acceptable values are 'cutoffs', 'weights',"
                    " and 'coarse_params'."
                )

    def set_scalers(self, scalers, rec=True):
        """
        Function to set custom scalers. 
         If rec is true, will calculate the reciprocal scalers
         scalers is a nested dict with the dataframe names 
         and scalers to be used. For example:
        {
            'DS1': {
                'DS2':{
                    'RT': {'x': [1, 2, 3], 'y': [.1, .2, .3]},
                    'MZ': {'x': [100, 200, 300], 'y' : [2, 4, 3]},
                    'Intensity': 1.001
                }
            }
        }
        While it seems clunky, it allows us to set multiple scalers at once. 
         If rec is True, and recipricol scalers are specified, they will be overwritten.
        :param scalers: dict
        :param rec: bool
        """
        for ds1 in scalers:
            # make scaler dicts if they don't exist
            if ds1 not in self.scalers:
                self.scalers[ds1] = {}
            for ds2 in scalers[ds1]:
                if ds2 not in self.scalers[ds1]:
                    self.scalers[ds1][ds2] = {}
                if rec:
                    if ds2 not in self.scalers:
                        self.scalers[ds2] = {}
                    if ds1 not in self.scalers[ds2]:
                        self.scalers[ds2][ds1] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in scaler names
                if not set(scalers[ds1][ds2].keys()).issubset(
                    set(self.descriptors.keys())
                ):
                    raise IndexError(
                        "The only acceptable columns for scalers are "
                        f"{self.descriptors.keys()}. "
                        f"You provided {list(scalers.keys())}."
                    )

                self._set_scalers(scalers[ds1][ds2], ds1, ds2, rec)

    def _set_scalers(self, scalers, ds1, ds2, rec):
        """
        Takes a list of scalers in a dictionary form 
         and sets them as scalers in the object.
        For example:
        {
            'RT': {'x': [1, 2, 3], 'y': [.1, .2, .3]},
            'MZ': {'x': [100, 200, 300], 'y' : [2, 4, 3]},
            'Intensity': {'x': [3, 4, 5], 'y' : [.1, .3, .4]}
        }
        Scalers should already be in their proper space

        :param scalers: Dict(str: Dict(str:list), 
         containing the scaler names as Keys and values
        :param ds1: str, name of Source dataset
        :param ds2: str, name of Target dataset
        :param rec: bool, Flag to create reciprocal or not (DS1->DS2, DS2->DS1)
            Convenience for manually specified scalers, not for use autogenerated
        """
        for k, v in scalers.items():
            values = list(zip(*v))
            df = pd.DataFrame({"x": values[0], "y": values[1]})
            self.scalers[ds1][ds2][k] = df
            if rec:
                if self.descriptors[k] in ["linear", "log10"]:
                    self.scalers[ds2][ds1][k] = pd.DataFrame(
                        {"x": df["x"] + df["y"], "y": -df["y"]}
                    )
                elif self.descriptors[k] == "ppm":
                    new_xs = df["x"] + df["x"] * df["y"] / 1000000
                    new_ys = -df["y"] * df["x"] / new_xs
                    self.scalers[ds2][ds1][k] = pd.DataFrame({"x": new_xs, "y": new_ys})
                else:
                    raise NotImplementedError(
                        "Not sure how to reverse scale the mode ", self.descriptors[k]
                    )

    def set_cutoffs(self, cutoffs, rec=True):
        """
        Function to set custom cutoffs. If rec is true, apply reciprocal cutoffs.
        {
            'DS1': {
                'DS2':{
                    'RT': 5.5,
                    'MZ': 5.5,
                    'Intensity': 5.0
                }
            }
        }
        :param cutoffs: dict
        :param rec: bool
        """
        cutoffs = cutoffs.copy()
        for ds1 in cutoffs:
            # make cutoff dicts if they don't exist
            if ds1 not in self.cutoffs:
                self.cutoffs[ds1] = {}
            for ds2 in cutoffs[ds1]:
                if ds2 not in self.cutoffs[ds1]:
                    self.cutoffs[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in cutoff names
                if not set(cutoffs[ds1][ds2].keys()).issubset(set(self.descriptors)):
                    raise IndexError(
                        f"The only acceptable columns for cutoffs are {self.descriptors.keys()}, and"
                        f" . You provided {cutoffs.keys()}."
                    )

                # replace 0's with inf, negating their use in cutoffs
                for key, value in cutoffs[ds1][ds2].items():
                    if value == 0:
                        cutoffs[ds1][ds2][key] = float("inf")

                self.cutoffs[ds1][ds2].update(cutoffs[ds1][ds2])
                if rec:
                    if ds2 not in self.cutoffs:
                        self.cutoffs[ds2] = {}
                    if ds1 not in self.cutoffs[ds2]:
                        self.cutoffs[ds2][ds1] = {}
                    self.cutoffs[ds2][ds1].update(cutoffs[ds1][ds2])

    def set_weights(self, weights, rec=True):
        """
        Function to set custom weights. If rec is true, apply reciprocal weights
        {
            'DS1': {
                'DS2':{
                    'RT': .5,
                    'MZ': 1,
                    'Intensity': 0
                }
            }
        }
        :param weights: dict
        :param rec: bool
        """
        weights = weights.copy()
        for ds1 in weights:
            # make weight dicts if they don't exist
            if ds1 not in self.weights:
                self.weights[ds1] = {}
            for ds2 in weights[ds1]:
                if ds2 not in self.weights[ds1]:
                    self.weights[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in weight names
                if not set(weights[ds1][ds2].keys()).issubset(set(self.descriptors)):
                    raise IndexError(
                        f"The only acceptable columns for weights are "
                        f"{self.descriptors.keys()}.  You provided {weights.keys()}."
                    )

                self.weights[ds1][ds2].update(weights[ds1][ds2])
                if rec:
                    if ds2 not in self.weights:
                        self.weights[ds2] = {}
                    if ds1 not in self.weights[ds2]:
                        self.weights[ds2][ds1] = {}
                    self.weights[ds2][ds1].update(weights[ds1][ds2])

    def set_coarse_params(
        self, coarse_params, rec=True
    ):  # pylint: disable=too-many-branches
        """
                Function to set custom cutoffs. If rec is true, apply reciprocal cutoffs
                {
                    'DS1': {
                        'DS2':{
                            'rt_minus': -.7,
                            'ppm_plus': 4
                        }
                    }
                }
                :param coarse_params: dict
                :param rec: bool
                """
        coarse_params = coarse_params.copy()
        for ds1 in coarse_params:
            # make coarse_params dicts if they don't exist
            if ds1 not in self.coarse_params:
                self.coarse_params[ds1] = {}
            for ds2 in coarse_params[ds1]:
                if ds2 not in self.coarse_params[ds1]:
                    self.coarse_params[ds1][ds2] = {}

                # Check if dataframes are in the aligner object
                if ds1 not in self.datasets or ds2 not in self.datasets:
                    raise IndexError(
                        f"Either {ds1} or {ds2} was not found in the datasets."
                    )

                # Check for typos in coarse_param names
                if not set(coarse_params[ds1][ds2].keys()).issubset(
                    set(key for key in self.descriptors)
                ):
                    raise IndexError(
                        f"The only acceptable columns for coarse_params are"
                        f" {set(key for key in self.descriptors)}. You"
                        f" provided {coarse_params.keys()}."
                    )

                self.coarse_params[ds1][ds2].update(coarse_params[ds1][ds2])
                if rec:
                    rec_dict = {}
                    for descriptor in coarse_params[ds1][ds2]:
                        new_feature = {}
                        if "upper" in coarse_params[ds1][ds2][descriptor]:
                            new_feature["lower"] = -coarse_params[ds1][ds2][descriptor][
                                "upper"
                            ]
                        if "lower" in coarse_params[ds1][ds2][descriptor]:
                            new_feature["upper"] = -coarse_params[ds1][ds2][descriptor][
                                "lower"
                            ]
                        rec_dict.update({descriptor: new_feature})

                    if ds2 not in self.coarse_params:
                        self.coarse_params[ds2] = {}
                    if ds1 not in self.coarse_params[ds2]:
                        self.coarse_params[ds2][ds1] = {}
                    self.coarse_params[ds2][ds1].update(rec_dict)

    def to_csv(self, *args, filepath=None, union_only=True, to_bytes=False):
        """
        Condenses network into a flat CSV file
        :param filepath, str, location to be saved
        :param *args, str, names of datasets that you want in the export file
        :param, union_only, bool, only cliques with all datasets should be returned
        """

        if len(args) < 1:
            targ_datasets = list(self.datasets)
        else:
            targ_datasets = list(args)

        df = self.results(*targ_datasets, union_only=union_only)
        if to_bytes:
            bytes_csv = df.to_csv().encode()
            return io.BytesIO(bytes_csv)
        if filepath is None:
            filepath = "results.csv"
        df.to_csv(filepath)
        return None

    def results(self, *args, union_only=True):
        """
        Returns a results dataframe
        """
        if len(args) < 1:
            targ_datasets = list(self.datasets)
        else:
            targ_datasets = list(args)

        df = self.produce_clique_table(targ_datasets, union_only)
        for ds in df.columns:
            df = df.merge(
                self.datasets[ds],
                how="left",
                left_on=ds,
                right_index=True,
                suffixes=(None, f"_{ds}"),
            )
        return df

    def compress_graph(self):
        """
        Compresses network into single highest rank non-directed graph
        """
        graph = self.graph.copy()
        selected_edges = [(u, v) for u, v, e in graph.edges(data=True) if e["rank"] > 0]
        for edge in selected_edges:
            graph.remove_edge(*edge)
        return graph.to_undirected(reciprocal=True)

    def produce_clique_table(self, targ_datasets, union_only):
        """
        Makes a dataframe where each row is a feature, and each column is a dataset
        Entry is null if there is no matching clique
        """
        df_indices = {ds: [] for ds in self.datasets}

        for clique in nx.find_cliques(self.compress_graph()):
            clique_ds = [feature_name.split("__")[0] for feature_name in clique]
            if union_only:
                if len(clique) >= len(targ_datasets):
                    # Here targ_datasets must be subset of clique_ds
                    if not all(ds in clique_ds for ds in targ_datasets):
                        continue
                else:
                    continue
            else:
                # Here we want any value in the clique to be in targ_datasets
                if not any(ds in targ_datasets for ds in clique_ds):
                    continue

            for df in df_indices.keys():
                for feature in clique:
                    ds_name, feature_name = feature.split("__")
                    if ds_name == df:
                        df_indices[df].append(feature_name)
                        break
                else:
                    df_indices[df].append(None)
        return pd.DataFrame.from_dict(df_indices)

    def align(self, ds1=None, ds2=None):
        """
        Performs all steps for generating the network
        """
        self.gen_anchors(ds1, ds2)
        self.gen_coarse(ds1, ds2)
        self.gen_scalers(ds1, ds2)
        self.gen_scaled_values(ds1, ds2)
        self.gen_stds(ds1, ds2)
        self.gen_matches(ds1, ds2)
        self.gen_graph(ds1, ds2)

    def gen_anchors(self, ds1=None, ds2=None):
        """
        Creates and sets anchors for all datasets specified.
        If either ds1 or ds2 is None, calculates anchors fall all datasets.
        """
        if ds1 is None or ds2 is None:
            datasets = self.datasets
        else:
            datasets = set()
            for ds in [ds1, ds2]:
                if isinstance(ds, str):
                    datasets.add(ds)
                else:
                    datasets.update(ds)
        coarse_params = copy.deepcopy(
            {k: self.default_coarse_params[k] for k in self.descriptors}
        )

        # copy the mode into the params
        for descr_param in coarse_params:
            coarse_params[descr_param].update({"mode": self.descriptors[descr_param]})

        for ds in datasets:
            LOGGER.info(f"Generating anchors for {ds}...")
            self.anchors[ds] = anchors(
                self.datasets[ds], coarse_params, self.remove_all
            )

    @dataset_loops("coarse_matches")
    def gen_coarse(self, ds1, ds2):
        """
        Generates all coarse anchors between datasets
        """
        LOGGER.info(f"Generating coarse matches for {ds1} -> {ds2}...")
        coarse_params = self.default_coarse_params.copy()
        try:
            coarse_params.update(self.coarse_params[ds1][ds2])
        except KeyError:
            pass
        # remove coarse params which aren't in the descriptors
        coarse_params = copy.deepcopy({k: coarse_params[k] for k in self.descriptors})
        for descr_param in coarse_params:
            coarse_params[descr_param].update({"mode": self.descriptors[descr_param]})
        # extract anchors, match, and convert results to a dataframe
        res = simple_match(
            self.datasets[ds1].loc[self.anchors[ds1], :],
            self.datasets[ds2].loc[self.anchors[ds2], :],
            coarse_params,
        )
        self.coarse_matches[ds1][ds2] = pd.DataFrame(res)
        self.coarse_matches[ds1][ds2].columns = [ds1, ds2]

    @dataset_loops("scalers")
    def gen_scalers(self, ds1, ds2):
        """
        Creates and saves scalers for rt, mz, and intensity.
        These autogenerated scalers will not overwrite scalers already in place.
        """
        LOGGER.info(f"Generating scalers for {ds1} -> {ds2}...")
        coarse_match = self.coarse_matches[ds1][ds2]
        ds1_matches = self.datasets[ds1].loc[coarse_match[ds1], :]
        ds2_matches = self.datasets[ds2].loc[coarse_match[ds2], :]
        if ds2 not in self.scalers[ds1]:
            self.scalers[ds1][ds2] = {}

        scalers = {}
        for (descr, mode) in self.descriptors.items():
            if descr not in self.scalers[ds1][ds2]:
                scaled = calc_scalers(
                    ds1_matches[descr].values,
                    ds2_matches[descr].values,
                    smoothing="lowess",
                    mode=mode,
                )
                scalers[descr] = zip(*scaled)
        self._set_scalers(scalers, ds1, ds2, rec=False)

    @dataset_loops("scaled_values")
    def gen_scaled_values(self, ds1, ds2):
        """
        Scales the rt, mz, and intensity values for all datasets
        """
        LOGGER.info(f"Scaling {ds1} -> {ds2}...")
        # extract anchors, match, and convert results to a dataframe
        scaled = pd.DataFrame(index=self.datasets[ds1].index)
        scalers = self.scalers[ds1][ds2]
        for (descr, mode) in self.descriptors.items():
            scaled[descr] = scale_to(
                self.datasets[ds1][descr],
                scalers[descr]["x"],
                scalers[descr]["y"],
                mode=mode,
            )
        self.scaled_values[ds1][ds2] = scaled

    @dataset_loops("stds")
    def gen_stds(self, ds1, ds2):
        """
        Calculates standard deviations to be used for matching.
        """
        LOGGER.info(f"Calculating variance from {ds1} -> {ds2}...")
        coarse_matches = self.coarse_matches[ds1][ds2]
        scalers = self.scalers[ds1][ds2]
        self.stds[ds1][ds2] = {}
        for (descr, mode) in self.descriptors.items():
            ds1_vals = self.datasets[ds1].loc[coarse_matches[ds1], descr]
            ds2_vals = self.datasets[ds2].loc[coarse_matches[ds2], descr]
            s_ds1_vals = scale_to(
                ds1_vals, scalers[descr]["x"], scalers[descr]["y"], mode=mode
            )
            if mode == "linear":
                res = (s_ds1_vals - ds2_vals).values
            if mode == "ppm":
                res = ((s_ds1_vals - ds2_vals) / s_ds1_vals * 1_000_000).values
            if mode == "log10":
                res = (np.log10(s_ds1_vals) - np.log10(ds2_vals)).values

            self.stds[ds1][ds2][descr] = np.std(res[~is_outlier(res)])

    @dataset_loops("matches")
    def gen_matches(self, ds1, ds2):
        """
        Generates best matches for each dataset
        """
        LOGGER.info(f"Matching {ds1} -> {ds2}...")
        # generate the weights and cutoffs
        cutoffs = self.default_cutoffs.copy()
        try:
            cutoffs.update(self.cutoffs[ds1][ds2])
        except KeyError:
            pass

        weights = self.default_weights.copy()
        try:
            weights.update(self.weights[ds1][ds2])
        except KeyError:
            pass

        self.matches[ds1][ds2] = score_match(
            self.scaled_values[ds1][ds2],
            self.datasets[ds2],
            self.stds[ds1][ds2],
            self.descriptors,
            cutoffs,
            weights,
        )

    @dataset_loops()
    def gen_graph(self, ds1, ds2):
        """
        Generates the network graph to be used for feature selection
        """
        LOGGER.info(f"Loading feature graph with {ds1} -> {ds2}...")
        # adding all features to the graph
        # If the node already exists (because it was added from a previous edge)\
        # it does not overwrite based on my tests -JGP
        nodes = []
        nodes = [f"{ds1}__{str(i)}" for i in self.datasets[ds1].index]
        self.graph.add_nodes_from(nodes)

        matches = self.matches[ds1][ds2]
        edges = []
        for col in matches:
            this_match = matches[col][~pd.isnull(matches[col])]
            sources = this_match.index
            sources = [ds1 + "__" + str(name) for name in sources]
            targets = this_match.values
            targets = [ds2 + "__" + str(name) for name in targets]
            edges.extend(
                [tuple([s, t, {"rank": col}]) for s, t in zip(sources, targets)]
            )
        self.graph.add_edges_from(edges)

    def report(
        self, datasets_1=None, datasets_2=None, filepath="report.pdf", to_bytes=False,
    ):
        """
        Creates a report of the alignment parameters.
        """
        if datasets_1 is None:
            datasets_1 = list(self.datasets)
        if isinstance(datasets_1, str):
            datasets_1 = [datasets_1]

        if datasets_2 is None:
            datasets_2 = list(self.datasets)
        if isinstance(datasets_2, str):
            datasets_2 = [datasets_2]
        if to_bytes:
            file_handle = io.BytesIO()
        else:
            file_handle = filepath

        # create scaling report
        with PdfPages(file_handle) as pdf:
            for ds1 in datasets_1:
                for ds2 in datasets_2:
                    if ds1 == ds2:
                        continue
                    eval_chart(self, ds1, ds2, show=False)
                    pdf.savefig()
                    plt.close("all")
                    eval_chart(self, ds1, ds2, show=False, results=True)
                    pdf.savefig()
                    plt.close("all")
        if to_bytes:
            file_handle.seek(0)
            return file_handle
        return None


def score_match(df1, df2, stds, descriptors, cutoffs, weights, top_n=3):
    """
    Returns scored best matches for a dataframe to another
    :param df1: DataFrame, first dataset
    :param df2: DataFrame, second dataset
    :param stds: Dict, standard deviations for "RT", "MT", and "Intensity"
    :param weights: Dict, Multipliers for weights
    :param cutoffs: Dict, Multipliers for finding cutoffs
    :param top_n: Int, number of potential matches to record
    :return: Dataframe, df1 index as index, and matching df2 index in top_n columns
    """
    results = []

    np1 = df1[list(descriptors)].values
    np2 = df2[list(descriptors)].values

    for row in np1:
        potential_hits = np.full(len(np2), True, dtype=bool)
        ds1_values = {}

        # generate values and bounds
        for i, (descr, mode) in enumerate(descriptors.items()):
            ds1_values[descr] = row[i]
            if mode == "linear":
                bounds = (
                    ds1_values[descr] - cutoffs[descr] * stds[descr],
                    ds1_values[descr] + cutoffs[descr] * stds[descr],
                )
            elif mode == "ppm":
                bounds = (
                    ds1_values[descr]
                    - (ds1_values[descr] / 1000000 * cutoffs[descr]) * stds[descr],
                    ds1_values[descr]
                    + (ds1_values[descr] / 1000000 * cutoffs[descr]) * stds[descr],
                )
            elif mode == "log10":
                bounds = (
                    10 ** (np.log10(ds1_values[descr]) - cutoffs[descr] * stds[descr]),
                    10 ** (np.log10(ds1_values[descr]) + cutoffs[descr] * stds[descr]),
                )
            else:
                raise NotImplementedError("Not sure how to handle the mode ", mode)

            # whittle hits down to list
            potential_hits = potential_hits & (
                (np2[:, i] > bounds[0]) & (np2[:, i] < bounds[1])
            )

        hits_index = df2.index[potential_hits]

        # score and sort results
        scores = []
        for hit in np2[potential_hits]:
            ds2_values = {descr: hit[i] for i, descr in enumerate(list(descriptors))}
            scores.append(
                calc_score(ds1_values, ds2_values, descriptors, stds, weights)
            )

        # append the top 3
        results.append([score for _, score in sorted(zip(scores, hits_index))][:top_n])
    results = pd.DataFrame(results, index=df1.index)
    return results


def calc_scalers(x1, x2, smoothing=None, mode="linear", **kwargs):
    """
    Creates and returns scalers from X1 and X2 coordinates
    i.e. Dataset 1 RTs - [7.1, 8.43, 9.5]
         Dataset 2 RTs - [7.24, 8.74, 9.7]
    Returns scalers in the proper space, e.g. log10 space for intensity, ppm for mz
         Dataset 1 RTs -     [7.1, 8.43, 9.5]
         Dataset 2 Scalers - [0.14, 0.31, 0.2]

    :param x1: Array or List
    :param x2: Array or List
    :param smoothing: str, Smoothing mode (None or "lowess" currently)
    :param mode: str, "linear" or "ppm"
    :return: Numpy array
    """
    x1 = np.array(x1)
    x2 = np.array(x2)

    if mode == "linear":
        y = x2 - x1
    elif mode == "ppm":
        y = (x2 - x1) / x1 * 1000000
    elif mode == "log10":
        y = np.log10(x2) - np.log10(x1)
        x1 = np.log10(x1)
    else:
        raise NotImplementedError(f"We do not have a scaling method for {mode}.")
    if smoothing is None:
        return (x1, y)

    # +/- <0.5 ppt random noise to break ties
    np.random.seed(0)
    noise = (np.random.random(size=len(x1)) - 0.5) * x1 / 10 ** 12

    # allow us to override frac via kwargs
    if "frac" not in kwargs:
        kwargs["frac"] = 0.1

    # If we get NaN for all our answers, add a bit to the frac and try again.
    # Or give up. It will be obvious in the report that it failed
    delta = 0.01 * (max(x1) - min(x1))
    scalers = lowess(y, x1 + noise, return_sorted=False, delta=delta, **kwargs)
    while np.isnan(scalers).all():
        kwargs["frac"] += 0.01
        if kwargs["frac"] >= 1.0:
            break
        scalers = lowess(y, x1 + noise, return_sorted=False, delta=delta, **kwargs,)
    return (x1, scalers)


def anchors(ds, match_params, remove_all=True, priority="Intensity"):
    """
    Creates a list of anchors for a given dataset.
    :param ds: DataFrame
    :param match_params: dict, parameters for matching
    :param remove_all: bool, flag to indicate removing all confusable feature
    :param priority: str, the column indicating the which features are highest quality
        i.e. higher intensity is a strong feature
    :return: list, the index names of the anchors
    """
    ds = ds[[*match_params.keys()]]
    if not remove_all:
        ds.sort_values(priority, ascending=False)

    # convert to numpy for performance
    index = ds.index
    ds = ds.values

    inclusion = set()
    exclusion = set()
    for i, row in enumerate(ds):

        # if a feature has been marked for removal, skip it.
        # important for non-remove all algorithm
        if not remove_all and index[i] in exclusion:
            continue
        bool_array = np.full(len(index), True, dtype=bool)
        for j, param in enumerate(match_params.values()):
            descr_val = row[j]
            if param["mode"] == "linear":
                lower = descr_val + param["lower"]
                upper = descr_val + param["upper"]
            elif param["mode"] == "ppm":
                lower = descr_val + (descr_val / 1_000_000 * param["lower"])
                upper = descr_val + (descr_val / 1_000_000 * param["upper"])
            elif param["mode"] == "log10":
                lower = descr_val * 10 ** param["lower"]
                upper = descr_val * 10 ** param["upper"]
            bool_array = bool_array & (ds[:, j] > lower) & (ds[:, j] < upper)

        results = ds[bool_array]
        result_indices = index[bool_array]
        if remove_all:
            if len(results) == 1:
                inclusion.add(index[i])
            elif len(results) > 1:
                # if there is more than one result, exclude all of them
                exclusion.update(set(result_indices))
        else:
            if len(results) == 1:
                # if the only result is itself add it and remove the others
                inclusion.add(index[i])
            elif len(results) > 1:
                if result_indices[0] == index[i]:
                    inclusion.add(index[i])
                    not_max = np.delete(result_indices, 0, 0)
                    exclusion.update(set(not_max))

    ds_anchors = list(inclusion - exclusion)
    ds_anchors.sort()
    return ds_anchors


def simple_match(df1, df2, match_params):
    """
    Simply iterates through rows and finds features that fall within a
     threshold in a reference dataset.
    In the future this function will be much more generalized
    :param df1: DataFrame, first dataset
    :param df2: DataFrame, second dataset
    :param match_params: dict, match parameters
    :return: list, 2D Nested containing indices for df1 and df2 matches
    """
    results = []

    np1 = df1[[*match_params.keys()]].values
    np2 = df2[[*match_params.keys()]].values
    for i, row in enumerate(np1):
        bool_array = np.full(len(df2), True, dtype=bool)
        for j, param in enumerate(match_params.values()):
            descr_val = row[j]
            if param["mode"] == "linear":
                lower = descr_val + param["lower"]
                upper = descr_val + param["upper"]
            elif param["mode"] == "ppm":
                lower = descr_val + (descr_val / 1_000_000 * param["lower"])
                upper = descr_val + (descr_val / 1_000_000 * param["upper"])
            elif param["mode"] == "log10":
                lower = descr_val * 10 ** param["lower"]
                upper = descr_val * 10 ** param["upper"]
            bool_array = bool_array & (np2[:, j] > lower) & (np2[:, j] < upper)

        for result in df2.index[bool_array]:
            results.append([df1.index[i], result])
    return results


def gen_intensity(dataset, ignore=None):
    """
    calculate an intensity column, using all columns that are not named '
     RT','MZ','Compound_ID'
    overwrites the existing dataframe with another one where Intensity is on the far end

    :param dataset: pandas Dataframe where there are intensities to calculate
    """
    if ignore:
        calc_df = dataset.drop(columns=ignore)
    calc_df.fillna(0, inplace=True)
    dataset["Intensity"] = calc_df.mean(axis=1)


def scale_to(x, x_scalers, y_scalers, mode="linear"):
    """
    Adjusts the x values to the scalers provided.

    :param x: Numpy Array, unscaled values
    :param x_scalers: Numpy Array, x values for the scalers
    :param y_scalers: Numpy Array, y values for the scalers
    :return: Numpy Array, scaled values for x
    """
    if mode == "linear":
        pass
    elif mode == "ppm":
        y_scalers = y_scalers * x_scalers / 1_000_000
    elif mode == "log10":
        x = np.log10(x)
    else:
        raise NotImplementedError(f"We do not have a scaling method for {mode}.")

    x_min_y = y_scalers.iloc[np.argmin(x_scalers)]
    x_max_y = y_scalers.iloc[np.argmax(x_scalers)]
    scaled = x + interpolate.interp1d(
        x_scalers, y_scalers, bounds_error=False, fill_value=(x_min_y, x_max_y)
    )(x)
    if isinstance(scaled, pd.Series):
        scaled = scaled.values
    if mode == "log10":
        scaled = 10 ** scaled
    return scaled


def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.
    Borrowed from https://stackoverflow.com/questions/22354094/pythonic-way-of-detecting
    -outliers-in-one-dimensional-observation-data

    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

    :param points: An numobservations by numdimensions array of observations
    :param thresh: The modified z-score to use as a threshold. Observations with
    :return: Numpy Array, A numobservations-length boolean array.
    """
    try:
        points = points.values
    except AttributeError:
        pass
    median = np.median(points, axis=0)
    diff = (points - median) ** 2
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh


def calc_score(v1, v2, descriptors, stds, weights):
    """
    Returns the score of two features, provided their parameters, standard deviations,
     and scoring weights
    :return: float, score
    """
    penalty = 0
    for (descr, mode) in descriptors.items():

        if mode == "linear":
            penalty += ((v1[descr] - v2[descr]) / stds[descr]) ** 2 * weights[descr]
        elif mode == "ppm":
            penalty += (
                (((v1[descr] - v2[descr]) / v1[descr]) * 1_000_000) / stds[descr]
            ) ** 2 * weights[descr]
        elif mode == "log10":
            penalty += (np.log10(v1[descr] / v2[descr]) / stds[descr]) ** 2 * weights[
                descr
            ]
        else:
            raise NotImplementedError("Not sure how to score mode ", mode)
    return penalty


def eval_chart(self, key1, key2, show=True, results=False):
    """
    Description
    """

    ds1 = self.datasets[key1]
    ds2 = self.datasets[key2]

    if not results:
        coarse_df = self.coarse_matches[key1][key2]
        ds1 = ds1.loc[coarse_df[key1]]
        ds2 = ds2.loc[coarse_df[key2]]
        title = "Coarse"
    else:
        # filter our datasets to just the matched features
        has_match = pd.notnull(self.matches[key1][key2][0])
        matches = self.matches[key1][key2].loc[has_match]
        ds1_matches = matches.index
        ds2_matches = matches[0]
        ds1 = ds1.loc[ds1_matches]
        ds2 = ds2.loc[ds2_matches]
        title = "Matched"

    plt.figure(figsize=(20, 12))
    for i, (descr, mode) in enumerate(self.descriptors.items()):
        descr_1 = ds1[descr].values
        descr_2 = ds2[descr].values
        scalers = self.scalers[key1][key2][descr].sort_values(by="x")

        # calculate coarse scalers, or pull matched scaled values
        if not results:
            s_descr_1 = scale_to(descr_1, scalers["x"], scalers["y"], mode=mode)
        else:
            s_descr_1 = self.scaled_values[key1][key2].loc[ds1_matches][descr]
        if mode == "linear":
            y_vals = descr_2 - descr_1
            s_y_vals = descr_2 - s_descr_1
        if mode == "ppm":
            y_vals = 1e6 * (descr_2 - descr_1) / descr_1
            s_y_vals = 1e6 * (descr_2 - s_descr_1) / s_descr_1
        if mode == "log10":
            y_vals = np.log10(descr_2) - np.log10(descr_1)
            s_y_vals = np.log10(descr_2) - np.log10(s_descr_1)
            descr_1 = np.log10(descr_1)

        # plot the deltas
        plt.subplot(len(self.descriptors), 3, 3 * i + 1)
        plt.scatter(
            descr_1, y_vals, alpha=0.4, s=3, rasterized=True,
        )
        plt.title(f"{descr} {key1} vs {key2}")
        plt.xlabel(f"{title} {descr} {key1}")
        plt.ylabel(f"{title} {descr} ({key2} - {key1}) : {mode}")
        plt.plot(scalers["x"], scalers["y"], color="red", rasterized=True)

        # plot the scaled
        plt.subplot(len(self.descriptors), 3, 3 * i + 2)
        plt.scatter(
            descr_1, s_y_vals, alpha=0.4, s=3, rasterized=True,
        )
        plt.title(f"{descr} {key1} vs {key2}")
        plt.xlabel(f"{title} {descr} {key1}")
        plt.ylabel(f"{title} {descr} ({key2} - {key1}) : {mode}")
        plt.axhline(y=0, color="red", rasterized=True)

        # plot histogram
        plt.subplot(len(self.descriptors), 3, 3 * i + 3)
        plt.hist(
            x=(s_y_vals), bins="auto", alpha=0.7, rwidth=0.85,
        )
        plt.title(f"{descr} Histograms")
        plt.xlabel(f"{title} {descr} ({key2} - {key1})")
        plt.tight_layout()

    if show:
        plt.show()
