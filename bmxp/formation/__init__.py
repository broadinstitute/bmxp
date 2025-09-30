from .formation_plots import *
import logging
import re
import io
import copy
from scipy.stats import zscore
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
import xlsxwriter
from bmxp.gravity import spearman, pearson
from bmxp import FMDATA, IMDATA

matplotlib.use("agg")

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

__version__ = "0.2.15"


def parse_formatted(dataset):
    """
    Given a formatted dataset, parses it into the abundances, injection metadata,
    feature metadta, and the sample names
    """
    if isinstance(dataset, pd.DataFrame):
        df = dataset
    else:
        try:
            df = pd.read_csv(dataset, header=None)
        except Exception as e:  # pylint: disable=broad-except
            try:
                df = pd.read_excel(dataset, engine="openpyxl", header=None)
            except:
                raise ValueError(
                    "Your file does not appear to be a csv or xlsx."
                ) from e

    # Find the index of the first non-missing value
    first_row = df.iloc[0, :].copy()
    first_row.replace("", np.nan, inplace=True)
    pivot_col = df.columns.get_loc(first_row.first_valid_index())

    first_col = df.iloc[:, 0].copy()
    first_col.replace("", np.nan, inplace=True)
    pivot_row = df.index.get_loc(first_col.first_valid_index())

    if pivot_row == 0 and pivot_col == 0:
        # there's no imdata, use the last recognized fmdata header as the pivot
        for header in FMDATA.values():
            if header in first_row.values:
                pivot_col = max(pivot_col, first_row[first_row == header].index[-1])

    pivot_value = df.iloc[pivot_row, pivot_col].lower()

    if pivot_col < 0 or pivot_row < 0:
        raise ValueError("Your dataset is blank on either the first row or column.")

    feature_metadata = df.iloc[pivot_row + 1 :, : pivot_col + 1].copy()
    feature_metadata.columns = df.iloc[pivot_row, : pivot_col + 1]

    inj_metadata = df.iloc[: pivot_row + 1, pivot_col:].copy()
    inj_metadata = inj_metadata.transpose().reset_index(drop=True)
    inj_metadata.columns = inj_metadata.iloc[0]
    inj_metadata = inj_metadata.drop(0)
    inj_metadata.columns = inj_metadata.columns.str.lower()
    # makes column name unique - https://stackoverflow.com/questions/24685012
    cols = pd.Series(inj_metadata.columns)
    for dup in inj_metadata.columns[inj_metadata.columns.duplicated(keep=False)]:
        cols[inj_metadata.columns.get_loc(dup)] = [
            dup + "." + str(d_idx) if d_idx != 0 else dup
            for d_idx in range(inj_metadata.columns.get_loc(dup).sum())
        ]
    inj_metadata.columns = cols
    # check and create necessary columns

    sample_data = df.iloc[pivot_row + 1 :, pivot_col + 1 :]
    sample_names = df.iloc[pivot_row, pivot_col + 1 :].astype(str)
    sample_data.columns = inj_metadata.index.astype(str) + " (" + sample_names + ")"
    sample_data.index = feature_metadata.index.copy()
    sample_data = sample_data.replace("", "nan")
    sample_data = sample_data.astype(float)
    inj_metadata = inj_metadata.drop(pivot_value, axis=1)
    return sample_data, inj_metadata, feature_metadata, sample_names


def report(
    sample_data,
    smdata,
    fmdata,
    dataset_name,
    sample_names=None,
    out_filepath=None,
    write_pdf=True,
):
    """
    Handles /formation endpoint
    """
    fmdata = fmdata.copy()
    smdata = smdata.copy()
    sample_data = sample_data.copy()

    # raise ValueError("'Compound_ID' was not found in your dataset.")
    if "Compound_ID" in fmdata:
        fmdata.set_index("Compound_ID", drop=True, inplace=True)
        sample_data.index = fmdata.index

    if sample_names is None:
        if "reporting_name" in smdata:
            sample_names = smdata["reporting_name"].copy()
        elif "program_id" in smdata:
            sample_names = smdata["program_id"].copy()
        else:
            sample_names = pd.Series(smdata.index, index=smdata.index.copy())
    # check and create necessary columns

    fill_values = {
        "Annotation_ID": "",
        "HMDB_ID": "",
        "Metabolite": "",
    }
    for fill_key in fill_values:
        if fill_key not in fmdata.columns:
            fmdata[fill_key] = fill_values[fill_key]

    fmdata.columns.name = None

    if "raw_file_name" not in smdata.columns:
        smdata["raw_file_name"] = smdata.index
        # raise ValueError("'raw_file_name' was not found in the sample metadata.")

    fill_values = {
        "sample_type": "unknown_type",
        "injection_type": "unknown_type",
        "injection_order": list(range(len(smdata))),
        "column_number": 1,
    }
    for fill_key in fill_values:
        if fill_key not in smdata.columns:
            smdata[fill_key] = fill_values[fill_key]

    # fill necessary missing values
    for sample_type in ("injection_type", "sample_type"):
        smdata[sample_type] = smdata[sample_type].fillna("")

    smdata.set_index("raw_file_name", drop=True, inplace=True)
    smdata = smdata.replace({"": float("nan"), "NA": float("nan")})
    smdata = smdata.convert_dtypes()
    # date_extracted usually does not convert to datetime correctly
    for date_type in ("date_extracted", "date_injected"):
        if date_type not in smdata.columns:
            continue
        try:
            smdata[date_type] = pd.to_datetime(smdata[date_type])
        except Exception as e:  # pylint: disable=broad-except
            raise ValueError(f"{date_type} contains one or more invalid dates.") from e

    try:
        smdata["injection_order"] = pd.to_numeric(smdata["injection_order"])
    except ValueError as e:
        character = re.search('".+"', str(e))
        if character:
            raise ValueError(
                f"Injection_order contains a non-numeric character: "
                f"{character.group()}."
            ) from e
        raise ValueError("Injection_order contains a non-numeric character.") from e

    transposed_s_data = sample_data.T
    transposed_s_data.index.name = "rawfile"
    transposed_s_data = transposed_s_data.apply(lambda x: x.fillna(x.min() / 2))
    # drop empty columns
    transposed_s_data = transposed_s_data.dropna(how="all", axis=1)
    s_data_zscores = zscore(transposed_s_data.to_numpy())
    s_data_zscores = np.nan_to_num(s_data_zscores)  # zscore is NaN if stdev is 0
    sample_pca = PCA(n_components=2)
    sample_pca_data = sample_pca.fit_transform(s_data_zscores)

    pca_df = pd.DataFrame(data=sample_pca_data, columns=["pc1", "pc2"])

    is_indices = fmdata.loc[fmdata["HMDB_ID"].str.lower() == "internal standard"].index
    internal_standard_s_data = sample_data.loc[is_indices]
    file_handle = io.BytesIO()
    palette = [
        # basel, from the R library 'yarrr'
        # https://cran.r-project.org/web/packages/yarrr/vignettes/piratepal.html
        (12, 91, 176, 0.7),  # blue
        (238, 0, 17, 0.7),  # red
        (21, 152, 61, 0.7),  # green
        (236, 87, 154, 0.7),  # pink
        (250, 107, 9, 0.7),  # orange
        (20, 155, 237, 0.7),  # light blue
        (161, 199, 32, 0.7),  # light green
        (254, 183, 11, 0.7),  # yellow
        (22, 160, 140, 0.7),  # teal
        (154, 112, 62, 0.7),  # brown
    ]
    nan_color = (0.55, 0.55, 0.55, 0.3)
    with PdfPages(file_handle) as pdf:
        plt.figure(figsize=(20, 12))
        plt.axis("off")
        plt.text(
            0.5,
            1,
            f"QC Report for {dataset_name}",
            ha="center",
            va="top",
            fontsize=24,
        )
        plt.text(
            0.5,
            0.9,
            f"Number of samples: {len(sample_data.columns)}",
            ha="center",
            va="top",
            fontsize=18,
        )
        plt.text(
            0.5,
            0.8,
            f"Number of features: {len(fmdata.index)}",
            ha="center",
            va="top",
            fontsize=18,
        )
        if "date_injected" in smdata.columns:
            try:
                date = str(smdata["date_injected"].min().strftime("%B %d, %Y"))
            except:
                date = "Unknown"
            plt.text(
                0.5,
                0.70,
                ("Date of first injection: " + date),
                ha="center",
                va="top",
                fontsize=18,
            )
        pdf.savefig()
        plt.close("all")
        # create zoomable PCA plot
        plot_formation_zoomable_plot(
            pca_df,
            "PCA of metabolites labeled by sample name",
            f"PC1 ({round(sample_pca.explained_variance_ratio_[0] * 100)}"
            "% variance explained)",
            f"PC2 ({round(sample_pca.explained_variance_ratio_[1] * 100)}"
            "% variance explained)",
            sample_names,
            pdf,
        )

        components_df = pd.DataFrame(
            data=sample_pca.components_.transpose(),
            columns=["pc1", "pc2"],
            index=transposed_s_data.columns,
        )
        # create zoomable loadings plots
        plot_formation_zoomable_plot(
            components_df,
            "Loadings plot labeled by Annotation_ID",
            "Principal Component 1",
            "Principal Component 2",
            fmdata.loc[transposed_s_data.columns, "Metabolite"],
            pdf,
            rasterized=True,
        )
        plot_formation_zoomable_plot(
            components_df,
            "Loadings plot labeled by Compound_ID",
            "Principal Component 1",
            "Principal Component 2",
            transposed_s_data.columns,
            pdf,
            rasterized=True,
        )

        # create PCAs colored by inj metadata
        num_columns = smdata.select_dtypes(include=["float64", "int64"]).columns
        date_columns = smdata.select_dtypes(include=["datetime"]).columns
        str_columns = smdata.select_dtypes(include=["string"]).columns
        cmap = copy.copy(plt.cm.get_cmap("magma"))
        cmap.set_bad(color=nan_color)
        graph_index = 0
        for i, column_name in enumerate(smdata.columns):
            if i % 6 == 0 and i != 0:
                plt.tight_layout()
                pdf.savefig()
                plt.close("all")
            if i % 6 == 0:
                plt.figure(figsize=(20, 12))

            plt.subplot(2, 3, i % 6 + 1)
            plt.xlabel(
                f"PC1 ({round(sample_pca.explained_variance_ratio_[0] * 100)}"
                "% variance explained)"
            )
            plt.ylabel(
                f"PC2 ({round(sample_pca.explained_variance_ratio_[1] * 100)}"
                "% variance explained)"
            )
            plt.title("PCA of metabolites - colored by " + str(column_name))

            if column_name in num_columns:
                col = smdata[column_name].astype("float").replace({pd.NA: np.nan})
                plt.scatter(
                    pca_df["pc1"],
                    pca_df["pc2"],
                    c=col,
                    s=50,
                    cmap=cmap,
                    plotnonfinite=True,
                    alpha=0.7,
                )
                if col.isnull().values.any():
                    plt.legend(["NA"])
                    axes = plt.gca()
                    leg = axes.get_legend()
                    leg.legend_handles[0].set_color(nan_color)
                plt.colorbar(label=column_name)

            elif column_name in str_columns:
                col = smdata[column_name].replace({pd.NA: "NA"})
                targets = col.unique()
                colors = {
                    targets[k]: tuple(x / 255 for x in palette[k % len(palette)][:-1])
                    + (0.7,)
                    for k in range(len(targets))
                }
                plt.scatter(
                    pca_df["pc1"],
                    pca_df["pc2"],
                    color=col.map(colors),
                    s=50,
                )
                legend_colors = {k: colors[k] for k in list(colors)[:20]}
                handles = [
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor=v,
                        label=k,
                        markersize=8,
                    )
                    for k, v in legend_colors.items()
                ]
                leg = plt.legend(title="color", handles=handles)
                plt.gca().add_artist(leg)

            elif column_name in date_columns:
                plt.scatter(
                    pca_df["pc1"],
                    pca_df["pc2"],
                    c=mdates.date2num(smdata[column_name]),
                    s=50,
                    cmap=cmap,
                    alpha=0.7,
                )
                colorbar = plt.colorbar(label=column_name)
                loc = mdates.AutoDateLocator()
                colorbar.ax.yaxis.set_major_locator(loc)
                colorbar.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
            graph_index = i
        # reset graph_index to create page break between PCAs and other plots
        graph_index = 0
        if not internal_standard_s_data.empty:
            is_ids = fmdata.loc[internal_standard_s_data.index, "Annotation_ID"]
            internal_standard_s_data.reset_index().apply(
                lambda row: plot_formation_line_plot(
                    row,
                    graph_index + 3 * row.name,
                    smdata,
                    palette,
                    pdf,
                    "injection_type",
                    ann_id=is_ids.iloc[row.name],
                ),
                axis=1,
            )
            graph_index += 3 * len(internal_standard_s_data.index)
            norm = internal_standard_s_data.apply(
                lambda row: row / row.median(), axis=1
            )
            line_colors = []
            num_inj_types = len(smdata["injection_type"].unique())
            for i in range(num_inj_types, num_inj_types + len(norm.index)):
                line_colors.append(
                    tuple(x / 255 for x in palette[i % len(palette)][:-1]) + (0.7,)
                )

            norm.reset_index().apply(
                lambda row: plot_formation_line_plot(
                    row,
                    graph_index,
                    smdata,
                    palette,
                    pdf,
                    "injection_type",
                    line_colors[row.name],
                    True,
                ),
                axis=1,
            )
            plt.ylim(bottom=0)
            leg_labels = [f"{k} - {v}" for k, v in zip(norm.index, is_ids)]
            handles_two = [
                Line2D([0], [0], color=v, label=k)
                for k, v in zip(leg_labels, line_colors)
            ]
            leg2 = plt.legend(title="color", handles=handles_two, loc=4)
            plt.gca().add_artist(leg2)
            graph_index += 3
        sample_medians = sample_data.fillna(0).median()
        plot_formation_line_plot(
            sample_medians,
            graph_index,
            smdata,
            palette,
            pdf,
            "sample_type",
            sample_names=sample_names,
        )
        graph_index += 3
        median_data = sample_data.apply(lambda col: col.fillna(col.min() / 2))
        median_data = median_data.apply(np.log)
        median_data = median_data.T
        unique_samp_types = list(smdata["sample_type"].unique())
        sample_colors = {
            key: tuple(x / 255 for x in value[:-1]) + (0.7,)
            for key, value in zip(unique_samp_types, palette)
        }
        colors = smdata["sample_type"].copy().astype(object).map(sample_colors)
        if len(median_data.index) <= 1500:
            for i in range(0, len(median_data.index), 150):
                plot_formation_quartile(
                    median_data.iloc[i : i + 150].T,
                    graph_index,
                    sample_colors,
                    colors[i : i + 150],
                    pdf,
                    sample_names[i : i + 150],
                )
                graph_index += 3
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close("all")
        inj_types = smdata["injection_type"].unique()
        pools = [
            inj_type.upper() for inj_type in inj_types if inj_type.startswith("pref")
        ]
        create_pools_cv_table(fmdata, pools, pdf)
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close("all")
    if out_filepath is None:
        out_filepath = f"{dataset_name}_formation_report.pdf"
    if write_pdf:
        with open(out_filepath, "wb") as f:
            f.write(file_handle.getbuffer())
    return file_handle


def report_from_formatted(
    dataset, dataset_name=None, out_filepath=None, write_pdf=True
):
    if dataset_name is None:
        if isinstance(dataset, pd.DataFrame):
            dataset_name = "Dataset"
        else:
            dataset_name = dataset.split("\\")[-1].split("/")[-1]

    sample_data, imdata, fmdata, sample_names = parse_formatted(dataset)

    return report(
        sample_data,
        imdata,
        fmdata,
        dataset_name,
        sample_names,
        out_filepath,
        write_pdf,
    )


def _sort_dataset(final_dataset):
    """Sort annotated metabolites by order, then, class MZ, RT;
    sort nontargeted by RT, MZ"""
    # make sure all annotated features are above nontargeted features
    final_dataset = final_dataset.sort_values(by=["Metabolite"])
    annotated = ~pd.isnull(final_dataset["Metabolite"])

    # force internal standards sort to top, missing orderNum next
    ordernum_present = "orderNum" in final_dataset
    final_dataset.loc[
        final_dataset["HMDB_ID"].str.lower() == "internal standard", "orderNum"
    ] = -99
    final_dataset.loc[annotated, "orderNum"] = final_dataset.loc[
        annotated, "orderNum"
    ].fillna(-98)

    # sort annotated if present
    try:
        final_dataset.loc[annotated, :] = (
            final_dataset.loc[annotated, :]
            .sort_values(
                by=[
                    "orderNum",
                    "superClass",
                    "mainClass",
                    "subClass",
                    "MZ_Calculated",
                    "RT",
                ],
                key=lambda col: (
                    col
                    if not pd.api.types.is_string_dtype(col.fillna(""))
                    else col.str.lower()
                ),
                na_position="first",
            )
            .values
        )
    except KeyError:
        pass

    # sort nontargeted separately
    final_dataset.loc[~annotated, :] = (
        final_dataset.loc[~annotated, :].sort_values(by=["RT", "MZ"]).values
    )
    if not ordernum_present:
        del final_dataset["orderNum"]

    return final_dataset


def harmonize_metadata(data, injectionset, sampleset):
    """Generates combined Sampleset and Injectionset metadata with columns matching
    data. Also removes not_used. If there are samples in data not present in either
    sampleset or injectionset, they are preserved and missing metadata is null.
    """
    warnings = []
    # rename column for older datasets
    if "reporting_name" not in injectionset.columns:
        injectionset = injectionset.rename(columns={"program_id": "reporting_name"})
    combined = injectionset.loc[data.columns].copy()
    combined["broad_id"] = combined["broad_id"].fillna("NA")

    s_broad_ids = combined.loc[combined["injection_type"] == "sample", "broad_id"]
    duplicates = set(s_broad_ids.loc[s_broad_ids.duplicated()]) - set([""])
    if len(duplicates) > 0:
        warnings.append(f"There are duplicate Broad_IDs: {', '.join(duplicates)}")

    sampleset_copy = sampleset.copy()
    sampleset_copy.loc["", :] = np.nan  # blank broad_id means no metadata
    combined = combined.join(sampleset_copy, on="broad_id", lsuffix="_inj")

    additional_meta = [
        col for col in sampleset.columns if col not in ["program_id", "broad_id"]
    ]

    # use injection set reporting_name for pools or missing program_ids if present
    combined["program_id"] = combined["program_id"].astype(object)
    pools = combined.index[combined["broad_id"].str.lower() == "pref"]
    missing = pd.isnull(combined["program_id"])
    if "reporting_name" in combined.columns:
        mismatches = combined.loc[
            (combined["program_id"] != combined["reporting_name"])
            & ~pd.isnull(combined["program_id"])
            & ~pd.isnull(combined["reporting_name"])
        ]
        if len(mismatches) > 0:
            warnings.append(
                "There are mismatches between program_ids and reporting_names. Please "
                + "review and confirm program_ids for the following samples before "
                + f"sharing results: {', '.join(mismatches.index)}"
            )
        combined.loc[missing, "program_id"] = combined.loc[missing, "reporting_name"]
    # name pools if possible
    elif pools.str.contains("PREF").all():
        combined.loc[pools, "program_id"] = pools.str.extract(r"(PREF\w*)").values
        warnings.append(
            "PREF program_id names were determined from injection_id names. Please "
            + "review and confirm program_id names before sharing results."
        )

    if pd.isnull(combined.loc[:, "program_id"]).any():
        warnings.append(
            "Not all program_ids could be filled. Please add any missing program_ids "
            "before sharing results."
        )
    combined[additional_meta] = combined[additional_meta].astype(object).fillna("NA")

    qcrole = IMDATA["QCRole"]
    qcrole_map = {
        "tube_blank": "QC-Tube_Blank",
        "blank": "QC-Blank",
        "not_used": "QC-Not_Used",
        "mm": "QC-Master_Mix",
        "ms2": "QC-MS2",
        "brpp": "QC-BRPP",
        "bridge_pref": "QC-Scaling_Pool",
    }

    if qcrole in combined.columns:
        for prefix, replacement in qcrole_map.items():
            mask = combined["injection_type"].str.startswith(prefix, na=False)
            combined.loc[mask, qcrole] = replacement
    else:
        combined[qcrole] = combined["injection_type"]

    # re-order
    injection_meta = combined.reindex(
        columns=[
            "date_extracted",
            "date_injected",
            "column_number",
            "injection_order",
            "injection_type",
            qcrole,
            "broad_id",
        ]
        + additional_meta
        + [
            "raw_file_name",
            "program_id",
        ]
    )
    injection_meta["raw_file_name"] = injection_meta.index
    # rename sample metadata and QCRole columns
    names = {col: col.lower().replace(" ", "_") for col in additional_meta}
    names[qcrole] = "sample_type"
    injection_meta = injection_meta.rename(columns=names)

    injection_meta.loc[:, "date_extracted"] = pd.to_datetime(
        injection_meta["date_extracted"],
        errors="coerce",
    )
    injection_meta.loc[:, "date_injected"] = pd.to_datetime(
        injection_meta["date_injected"],
        errors="coerce",
    )

    return injection_meta, warnings


def fill_f_mdata(metadata, annotations):
    """combine feature metadata from dataset and database queries for final results"""
    # request annotations and their HMBD IDs
    main_cols = ["__annotation_id", "Compound_ID", "MZ", "RT"]
    additional_meta = [col for col in metadata.columns if col not in main_cols]
    feature_meta = metadata.loc[:, main_cols + additional_meta]
    feature_meta["superClass"] = None
    feature_meta["mainClass"] = None
    feature_meta["subClass"] = None
    feature_meta["Adduct_Priority"] = None
    feature_meta["HMDB_ID"] = None
    feature_meta["HMDB_specificity (1=match; 2=representative)"] = None
    feature_meta["Annotation_ID"] = None
    feature_meta["Metabolite"] = None

    # move method to front
    method_col = ""
    if "Method" in feature_meta:
        method_col = feature_meta.pop("Method")
    feature_meta.insert(0, "Method", method_col)

    for key, value in annotations.items():
        if key not in feature_meta["__annotation_id"].unique():
            continue
        to_fill = feature_meta["__annotation_id"] == key
        for anno_key in value:
            feature_meta.loc[to_fill, anno_key] = value[anno_key]

    # format labels for each column; None = default formatting
    # feature_meta = feature_meta.reset_index(drop=True)
    return feature_meta


def label_most_abundant(data, metadata, corr_method="spearman"):
    """Creates Primary and Corr_To_Primary columns given a clustered dataset. "Primary"
    feature is annotated or has the highest mean abundance.

    :param data: DataFrame, data columns that were used for clustering
    :param metadata: DataFrame, contains columns "Annotation_ID", "Cluster_Num", and
    "Cluster_Size" at minimum
    :param corr_method: str, correlation method used in clustering, "spearman" or
    "pearson", defaults to "spearman"
    :return: DataFrame, Primary and Corr_To_Primary columns
    """
    metadata = metadata.copy()
    # explicitly cast as str in case the Annotation_ID column is empty
    metadata["Annotation_ID"] = metadata["Annotation_ID"].fillna("")
    # identify the primary adduct
    for cluster_num in metadata["Cluster_Num"].unique():
        if not cluster_num and cluster_num != 0:
            continue
        cluster_index = metadata.index[metadata["Cluster_Num"] == cluster_num]
        key = (metadata["Cluster_Num"] == cluster_num).values & (
            metadata["Annotation_ID"] > ""
        ).values

        if key.any():
            metadata.loc[key, "Primary"] = True
            max_idx = metadata.index[key][0]
        else:
            max_idx = data.loc[cluster_index, :].fillna(0).mean(axis="columns").idxmax()
            metadata.loc[max_idx, "Primary"] = True

        for feature in cluster_index:
            arr1 = data.loc[feature, :].values.astype(np.float64)
            arr2 = data.loc[max_idx, :].values.astype(np.float64)
            if corr_method == "spearman":
                metadata.loc[feature, "Corr_To_Primary"] = spearman(
                    arr1, arr2, "drop", legacy_mode=True
                )
            else:
                metadata.loc[feature, "Corr_To_Primary"] = pearson(arr1, arr2)

    # label the singletons
    metadata.loc[metadata["Cluster_Size"] == 1, "Primary"] = True

    return metadata[["Primary", "Corr_To_Primary"]]


def apply_pool_instructions(data, feature_instructions, inplace=False):
    """Apply all pool instructions"""
    if feature_instructions is None:
        return data
    if not inplace:
        data = data.copy()
    for feature, instructions in feature_instructions.items():
        if feature not in data.index:
            continue
        for inst in instructions:
            if inst["op1"] not in data.columns or (
                inst["instruction"] != "delete" and inst["op2"] not in data.columns
            ):
                # we're probably doing final formatting and already removed some unused
                # pools, so just move on
                continue
            if inst["instruction"] == "swap":
                op1_value = data.loc[feature, inst["op1"]]
                data.loc[feature, inst["op1"]] = data.loc[feature, inst["op2"]]
                data.loc[feature, inst["op2"]] = op1_value
            elif inst["instruction"] == "copy":
                # copy op1 to op2
                data.loc[feature, inst["op2"]] = data.loc[feature, inst["op1"]]
            elif inst["instruction"] == "delete":
                data.loc[feature, inst["op1"]] = np.nan

    return data


def filter_samples_mask(data, smdata, formation_params):
    """
    Returns Boolean Mask of samples to keep
    the data columns must match the sample metadata index
    handles the old way, with 'only_prefs_samples' and 'keep_not_used'
    or the new way, with 'include_types'
    """
    warnings = []
    if not smdata.index.equals(data.columns):
        raise IndexError("Your Sample Metadata does not match your data columns.")

    if "include_types" in formation_params:
        if not formation_params["include_types"]:
            return np.array([True] * len(smdata)), warnings
        to_keep = np.array([False] * len(smdata))
        # Update to_keep for matching injection types
        for inj_type in formation_params["include_types"]:
            to_keep |= smdata["injection_type"].str.startswith(inj_type).values
    else:
        to_keep = np.array([True] * len(smdata))
        # filter out anything not labeled a pref or sample
        if formation_params["only_prefs_samples"]:
            to_keep = (
                smdata["injection_type"].str.startswith("pref").values
                | smdata["injection_type"].str.startswith("sample").values
            )

        # filter out not_used
        if not formation_params["keep_not_used"]:
            not_used = smdata["injection_type"].str.startswith("not_used").values
            to_keep = to_keep & ~not_used
    return to_keep, warnings


def filter_features_mask(data, smdata, fmdata, form_params):
    """
    Generates a filtering boolean key for data
    """
    warnings = []
    to_keep = np.array([True] * len(fmdata))

    if not fmdata.index.equals(data.index):
        raise IndexError("Your Feature Metadata does not match your data index.")

    pref_as = smdata.index[smdata["injection_type"] == "prefa"]
    pref_bs = smdata.index[smdata["injection_type"] == "prefb"]

    if form_params["filter_by_pref_missing"]:
        miss_column = f"PREF Missing (of {len(pref_as) + len(pref_bs)})"
        if form_params["missing_as_percent"]:
            missing_cutoff = np.ceil(
                (len(pref_as) + len(pref_bs)) * form_params["missing_cutoff"] / 100
            )
        else:
            missing_cutoff = form_params["missing_cutoff"]
        to_keep = to_keep & (fmdata[miss_column] <= missing_cutoff)

    if form_params["filter_by_pref_cv"]:
        if len(pref_as) > 0:
            to_keep = to_keep & (fmdata["PREFA CVs"] <= form_params["cv_cutoff"] / 100)
        if len(pref_bs) > 0:
            to_keep = to_keep & (fmdata["PREFB CVs"] <= form_params["cv_cutoff"] / 100)

    if form_params["filter_by_clusters"] and "Primary" in fmdata.columns:
        to_keep = to_keep & fmdata["Primary"].fillna(False)

    # annotated compounds don't get filtered by CV or missing PREFs
    to_keep = to_keep | pd.notnull(fmdata["Annotation_ID"])

    # but everything by non_quant
    if form_params["filter_by_nonquant"]:
        non_quant = fmdata["Non_Quant"].fillna(False).astype(bool)
        to_keep = to_keep & ~non_quant

    # find all duplicated annotations that weren't filtered as non_quant
    annotated = to_keep & pd.notnull(fmdata["Annotation_ID"])
    is_duplicated = fmdata.loc[annotated, "Annotation_ID"].duplicated(keep=False)
    duplicates = fmdata.loc[is_duplicated.index[is_duplicated]]
    if "Adduct" in fmdata.columns and not form_params["keep_adducts"] == "all":
        if form_params["keep_adducts"] == "top":
            # only keep annotation if top priority adduct
            fmdata_adduct_priority_null_vals = (
                fmdata["Adduct_Priority"].fillna(False).astype(bool)
            )
            fmdata_adduct_priority_defined = fmdata[fmdata_adduct_priority_null_vals]
            fmdata_top_adduct = fmdata_adduct_priority_defined["Adduct_Priority"].map(
                lambda x: x.split(",")[0]
            )
            to_drop = fmdata_adduct_priority_defined["Adduct"] != fmdata_top_adduct
            to_drop_idx = to_drop.index[to_drop]
            to_keep[to_drop_idx] = False
        else:
            for anno in set(duplicates["Annotation_ID"]):
                group = duplicates.loc[duplicates["Annotation_ID"] == anno]
                # pick highest priority adduct, if possible
                priorities = group.loc[group.index[0], "Adduct_Priority"]
                if not priorities:
                    continue
                priorities = priorities.split(",")
                for adduct in priorities:
                    if adduct in group["Adduct"].values:
                        to_drop = group["Adduct"] != adduct
                        to_drop_idx = to_drop.index[to_drop]
                        to_keep[to_drop_idx] = False
                        break

    # find annotations that are still duplicated
    if form_params["filter_by_extraction_method"]:
        annotated = to_keep & pd.notnull(fmdata["Annotation_ID"])
        is_duplicated = fmdata.loc[annotated, "Annotation_ID"].duplicated(keep=False)
        duplicates = fmdata.loc[is_duplicated.index[is_duplicated]]

        for anno in set(duplicates["Annotation_ID"]):
            group = duplicates.loc[duplicates["Annotation_ID"] == anno]
            # drop QI annotation(s) if there is a TF annotation
            if (
                form_params["feature_priority"][0]
                in group["__extraction_method"].values
            ):
                to_drop = (
                    group["__extraction_method"] == form_params["feature_priority"][1]
                )
                to_drop_idx = to_drop.index[to_drop]
                to_keep[to_drop_idx] = False

    # finally, if there are still duplicates, warn the user
    annotated = to_keep & pd.notnull(fmdata["Annotation_ID"])
    cols_to_check = ["Annotation_ID"]
    if "Adduct" in fmdata.columns and form_params["keep_adducts"] == "all":
        cols_to_check.append("Adduct")
    is_duplicated = fmdata.loc[annotated].duplicated(subset=cols_to_check, keep=False)
    if sum(is_duplicated) > 0:
        warnings.append(
            "There are duplicate annotations that could not be filtered by adduct "
            + "or extraction method. Please remove duplicate annotations before "
            + "sharing results."
        )

    return to_keep.values, warnings


def feature_qc(data, smdata, fmdata):
    """
    Generates QC data and adds to feature metadata
    """
    pref_as = smdata.index[smdata["injection_type"] == "prefa"]
    pref_bs = smdata.index[smdata["injection_type"] == "prefb"]
    miss_column = f"PREF Missing (of {len(pref_as) + len(pref_bs)})"

    if len(pref_as) > 0:
        fmdata["PREFA CVs"] = data.loc[:, pref_as].std(axis=1) / data.loc[
            :, pref_as
        ].mean(axis=1)
    else:
        fmdata["PREFA CVs"] = None
    if len(pref_bs) > 0:
        fmdata["PREFB CVs"] = data.loc[:, pref_bs].std(axis=1) / data.loc[
            :, pref_bs
        ].mean(axis=1)
    else:
        fmdata["PREFB CVs"] = None
    fmdata[miss_column] = data.loc[:, pref_as.union(pref_bs)].isnull().sum(axis=1)
    return fmdata, miss_column


def combine(data, smdata, fmdata, fmdata_formats):
    """Returns a formatted dataset and any warnings"""
    warnings = []
    abundance_threshold = 1
    fmdata = fmdata.copy()
    smdata = smdata.copy()

    # check mdata and data have same index
    dataset_index = data.index
    if set(fmdata.index) != set(data.index):
        warnings.append(
            "Your metadata and data have mismatched features, likely because this"
            " is an old dataset. You might have missing features..."
        )
        dataset_index = fmdata.index.intersection(data.index).copy()
        data = data.loc[dataset_index, :]
        fmdata = fmdata.loc[dataset_index, :]

    # sort and drop index (which was Compound_ID); Compound_ID should be in feature metadata
    fmdata = fmdata.loc[dataset_index, :].reset_index(drop=True)
    data.reset_index(drop=True, inplace=True)
    data[data < abundance_threshold] = np.nan

    # fill non_quants and cast as object
    if "Non_Quant" not in fmdata.columns:
        fmdata["Non_Quant"] = False
    fmdata.fillna({"Non_Quant": False}, inplace=True)
    fmdata["Non_Quant"] = fmdata["Non_Quant"].astype("object")

    # fill values less than 1 with 1 and round
    data = data.clip(lower=1).fillna(0)
    data = data.round(0).astype(np.int64)
    data = data.replace(0, np.nan)

    # keep additional columns but drop columns that are for internal use only
    additional_cols = [
        col for col in fmdata if col not in fmdata_formats and not col.startswith("__")
    ]
    fmdata = fmdata.reindex(columns=additional_cols + list(fmdata_formats))

    final_dataset = pd.concat([fmdata, data], axis="columns")
    final_dataset = _sort_dataset(final_dataset)

    # set headers as top row
    injection_meta = smdata.T
    final_dataset = final_dataset.reset_index(drop=True).reset_index()
    final_dataset = final_dataset.astype("object")
    final_dataset.loc[-1, final_dataset.columns] = final_dataset.columns
    final_dataset.loc[-1, injection_meta.columns] = injection_meta.loc["program_id"]
    final_dataset = final_dataset.sort_index()
    final_dataset = final_dataset.rename(index={-1: "Metabolite"})
    injection_meta = injection_meta.drop(index="program_id")
    injection_meta.loc[:, "Metabolite"] = injection_meta.index.str.capitalize()
    original_columns = final_dataset.columns
    final_dataset = pd.concat([injection_meta, final_dataset], axis="index")
    final_dataset = final_dataset.loc[:, original_columns]
    final_dataset = final_dataset.fillna("")
    return final_dataset, warnings


def to_excel(data, fmdata_formats, sdmata_formats, miss_column, method_name="Default"):
    """
    Returns a BytesIO object of a formatted excel sheet
    """

    output = io.BytesIO()
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet(method_name)
    formats = {
        "date": workbook.add_format(
            {
                "font_name": "Arial",
                "font_size": 9,
                "align": "center",
                "num_format": "m/d/yyyy",
            }
        ),
        "int": workbook.add_format(
            {
                "font_name": "Arial",
                "font_size": 9,
                "num_format": "0",
            }
        ),
        "center": workbook.add_format(
            {"font_name": "Arial", "font_size": 9, "align": "center"}
        ),
        "float2": workbook.add_format(
            {"font_name": "Arial", "font_size": 9, "num_format": "0.00"}
        ),
        "float4": workbook.add_format(
            {"font_name": "Arial", "font_size": 9, "num_format": "0.0000"}
        ),
        "label": workbook.add_format(
            {
                "font_name": "Arial",
                "font_size": 9,
                "italic": True,
                "align": "right",
            }
        ),
        "header": workbook.add_format(
            {"font_name": "Arial", "font_size": 9, "bold": True}
        ),
        "header_center": workbook.add_format(
            {"font_name": "Arial", "font_size": 9, "bold": True, "align": "center"}
        ),
        "default": workbook.add_format({"font_name": "Arial", "font_size": 9}),
        "percent": workbook.add_format(
            {"font_name": "Arial", "font_size": 9, "num_format": "0%"}
        ),
    }

    column_formats = {
        header: formats[label] for header, label in fmdata_formats.items()
    }
    row_formats = {index: formats[label] for index, label in sdmata_formats.items()}

    # column widths
    worksheet.set_column_pixels(
        data.columns.get_loc("Compound_ID"),
        data.columns.get_loc("Compound_ID"),
        86,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("superClass"),
        data.columns.get_loc("superClass"),
        80,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("mainClass"),
        data.columns.get_loc("mainClass"),
        80,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("subClass"),
        data.columns.get_loc("subClass"),
        80,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc(miss_column),
        data.columns.get_loc(miss_column),
        80,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("Annotation_ID"),
        data.columns.get_loc("Annotation_ID"),
        90,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("Non_Quant"),
        data.columns.get_loc("Non_Quant"),
        80,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("HMDB_ID"),
        data.columns.get_loc("HMDB_ID"),
        102,
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("HMDB_specificity (1=match; 2=representative)"),
        data.columns.get_loc("HMDB_specificity (1=match; 2=representative)"),
        270,
    )
    worksheet.set_column(
        data.columns.get_loc("Metabolite"),
        data.columns.get_loc("Metabolite"),
        data["Metabolite"].str.len().max(),
    )
    worksheet.set_column_pixels(
        data.columns.get_loc("Metabolite") + 1,
        len(data.columns) - 1,
        150,
    )

    for i, row in enumerate(data.values):
        index_label = data.index[i]
        if index_label in row_formats.keys():  # injection metadata
            worksheet.set_row_pixels(i, 16, row_formats[index_label])
        else:  # default to int
            worksheet.set_row_pixels(i, 16, formats["int"])

        for j, val in enumerate(row):
            column_label = data.columns[j]
            # metadata headers
            if index_label == "Metabolite" and column_label in column_formats.keys():
                worksheet.write(i, j, val, formats["header"])
            # data headers
            elif index_label == "Metabolite":
                worksheet.write(i, j, val, formats["header_center"])
            # row labels
            elif index_label in row_formats.keys() and column_label == "Metabolite":
                worksheet.write(i, j, val, formats["label"])
            # feature metadata
            elif column_label in column_formats.keys():
                worksheet.write(i, j, val, column_formats[column_label])
            else:
                worksheet.write(i, j, val)
    workbook.close()
    return output
