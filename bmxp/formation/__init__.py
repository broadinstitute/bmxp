import logging
import re
import io
import copy
import pandas as pd
import numpy as np
from scipy.stats import zscore
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import xlsxwriter
from bmxp.gravity import spearman, pearson

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

__version__ = "0.0.4"


def report(
    sample_data,
    inj_metadata,
    feature_metadata,
    dataset_name,
    sample_names=None,
    out_filepath=None,
    write_pdf=True,
):
    """
    Handles /formation endpoint
    """
    feature_metadata = feature_metadata.copy()
    inj_metadata = inj_metadata.copy()
    sample_data = sample_data.copy()

    # raise ValueError("'Compound_ID' was not found in your dataset.")
    if "Compound_ID" in feature_metadata:
        feature_metadata.set_index("Compound_ID", drop=True, inplace=True)
        sample_data.index = feature_metadata.index

    if sample_names is None:
        if "reporting_name" in inj_metadata:
            sample_names = inj_metadata["reporting_name"].copy()
        elif "program_id" in inj_metadata:
            sample_names = inj_metadata["program_id"].copy()
        else:
            sample_names = pd.Series(
                inj_metadata.index, index=inj_metadata.index.copy()
            )
    # check and create necessary columns

    fill_values = {
        "Annotation_ID": "",
        "HMDB_ID": "",
        "Metabolite": "",
    }
    for fill_key in fill_values:
        if fill_key not in feature_metadata.columns:
            feature_metadata[fill_key] = fill_values[fill_key]

    feature_metadata.columns.name = None

    if "raw_file_name" not in inj_metadata.columns:
        inj_metadata["raw_file_name"] = inj_metadata.index
        # raise ValueError("'raw_file_name' was not found in the sample metadata.")

    fill_values = {
        "sample_type": "unknown_type",
        "injection_type": "unknown_type",
        "injection_order": list(range(len(inj_metadata))),
        "column_number": 1,
    }
    for fill_key in fill_values:
        if fill_key not in inj_metadata.columns:
            inj_metadata[fill_key] = fill_values[fill_key]

    # fill necessary missing values
    for sample_type in ("injection_type", "sample_type"):
        inj_metadata[sample_type] = inj_metadata[sample_type].fillna("")

    inj_metadata.set_index("raw_file_name", drop=True, inplace=True)
    inj_metadata = inj_metadata.replace({"": float("nan"), "NA": float("nan")})
    inj_metadata = inj_metadata.convert_dtypes()
    # date_extracted usually does not convert to datetime correctly
    for date_type in ("date_extracted", "date_injected"):
        if date_type not in inj_metadata.columns:
            continue
        try:
            inj_metadata[date_type] = pd.to_datetime(inj_metadata[date_type])
        except Exception as e:  # pylint: disable=broad-except
            raise ValueError(f"{date_type} contains one or more invalid dates.") from e

    try:
        inj_metadata["column_number"] = pd.to_numeric(inj_metadata["column_number"])
    except ValueError as e:
        character = re.search('".+"', str(e))
        if character:
            raise ValueError(
                f"Column_number contains a non-numeric character: "
                f"{character.group()}."
            ) from e
        raise ValueError("Column_number contains a non-numeric character.") from e

    try:
        inj_metadata["injection_order"] = pd.to_numeric(inj_metadata["injection_order"])
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
    transposed_s_data = zscore(transposed_s_data)
    sample_pca = PCA(n_components=2)
    sample_pca_data = sample_pca.fit_transform(transposed_s_data)

    pca_df = pd.DataFrame(data=sample_pca_data, columns=["pc1", "pc2"])

    is_indices = feature_metadata.loc[
        feature_metadata["HMDB_ID"].str.lower() == "internal standard"
    ].index
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
            f"Number of features: {len(feature_metadata.index)}",
            ha="center",
            va="top",
            fontsize=18,
        )
        if "date_injected" in inj_metadata.columns:
            plt.text(
                0.5,
                0.70,
                (
                    "Date of first injection: "
                    f"{inj_metadata['date_injected'].min().strftime('%B %d, %Y')}"
                ),
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
            feature_metadata.loc[transposed_s_data.columns, "Metabolite"],
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
        num_columns = inj_metadata.select_dtypes(include=["float64", "int64"]).columns
        date_columns = inj_metadata.select_dtypes(include=["datetime"]).columns
        str_columns = inj_metadata.select_dtypes(include=["string"]).columns
        cmap = copy.copy(plt.cm.get_cmap("magma"))
        cmap.set_bad(color=nan_color)
        graph_index = 0
        for i, column_name in enumerate(inj_metadata.columns):
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
                col = inj_metadata[column_name].astype("float").replace({pd.NA: np.nan})
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
                    leg.legendHandles[0].set_color(nan_color)
                plt.colorbar(label=column_name)

            elif column_name in str_columns:
                col = inj_metadata[column_name].replace({pd.NA: "NA"})
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
                    c=mdates.date2num(inj_metadata[column_name]),
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
        is_ids = feature_metadata.loc[internal_standard_s_data.index, "Annotation_ID"]
        internal_standard_s_data.reset_index().apply(
            lambda row: plot_formation_line_plot(
                row,
                graph_index + 3 * row.name,
                inj_metadata,
                palette,
                pdf,
                "injection_type",
                ann_id=is_ids[row.name],
            ),
            axis=1,
        )
        graph_index += 3 * len(internal_standard_s_data.index)
        norm = internal_standard_s_data.apply(lambda row: row / row.median(), axis=1)
        line_colors = []
        num_inj_types = len(inj_metadata["injection_type"].unique())
        for i in range(num_inj_types, num_inj_types + len(norm.index)):
            line_colors.append(
                tuple(x / 255 for x in palette[i % len(palette)][:-1]) + (0.7,)
            )

        norm.reset_index().apply(
            lambda row: plot_formation_line_plot(
                row,
                graph_index,
                inj_metadata,
                palette,
                pdf,
                "injection_type",
                line_colors[row.name],
                True,
            ),
            axis=1,
        )
        if not internal_standard_s_data.empty:
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
            inj_metadata,
            palette,
            pdf,
            "sample_type",
            sample_names=sample_names,
        )
        graph_index += 3
        median_data = sample_data.apply(lambda col: col.fillna(col.min() / 2))
        median_data = median_data.apply(np.log)
        median_data = median_data.T
        unique_samp_types = list(inj_metadata["sample_type"].unique())
        sample_colors = {
            key: tuple(x / 255 for x in value[:-1]) + (0.7,)
            for key, value in zip(unique_samp_types, palette)
        }
        colors = inj_metadata["sample_type"].map(sample_colors)
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
        inj_types = inj_metadata["injection_type"].unique()
        pools = [
            inj_type.upper() for inj_type in inj_types if inj_type.startswith("pref")
        ]
        create_pools_cv_table(feature_metadata, pools, pdf)
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
    if isinstance(dataset, pd.DataFrame):
        if dataset_name is None:
            dataset_name = "Dataset"
        df = dataset
    else:
        if dataset_name is None:
            dataset_name = dataset.split("\\")[-1].split("/")[-1]
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
    row = df.iloc[0, :].copy()
    row.replace("", np.nan, inplace=True)
    pivot_col = df.columns.get_loc(row.first_valid_index())

    row = df.iloc[:, 0].copy()
    row.replace("", np.nan, inplace=True)
    pivot_row = df.index.get_loc(row.first_valid_index())

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

    return report(
        sample_data,
        inj_metadata,
        feature_metadata,
        dataset_name,
        sample_names,
        out_filepath,
        write_pdf,
    )


def plot_formation_zoomable_plot(
    pc_df, title, xaxis, yaxis, labels, pdf, rasterized=False
):
    """
    Creates large zoomable scatter plots for formation
    """
    fig = plt.figure(figsize=(20, 20), dpi=400)
    fig.set_rasterized(rasterized)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.title(title)

    plt.scatter(pc_df["pc1"], pc_df["pc2"], c="#000", s=1, alpha=0.5)
    # add jitter to labels to prevent overlapping
    np.random.seed(0)
    pc1_stdev = 0.0005 * (max(pc_df["pc1"]) - min(pc_df["pc1"]))
    pc1_jitter = pc_df["pc1"] + np.random.randn(len(pc_df["pc1"])) * pc1_stdev
    pc2_stdev = 0.0005 * (max(pc_df["pc2"]) - min(pc_df["pc2"]))
    pc2_jitter = pc_df["pc2"] + np.random.randn(len(pc_df["pc2"])) * pc2_stdev

    for i, label in enumerate(labels):
        if pd.isnull(label):
            continue
        plt.annotate(
            label,
            xy=(pc_df["pc1"][i], pc_df["pc2"][i]),
            xycoords="data",
            xytext=(pc1_jitter[i], pc2_jitter[i]),
            textcoords="data",
            fontsize=2,
            ha="left",
            rotation=30,
        )
    plt.tight_layout()
    pdf.savefig()
    plt.close("all")


def plot_formation_line_plot(
    row,
    count,
    inj_metadata,
    palette,
    pdf,
    types,
    line_color="grey",
    is_aggregate=False,
    ann_id="",
    sample_names=None,
):
    """
    Creates line plot for formation
    """
    unique_types = list(inj_metadata[types].unique())
    colors = {}
    if not pd.isnull(unique_types[0]):
        if "sample" in unique_types:
            unique_types.remove("sample")
        if "Sample" in unique_types:
            unique_types.remove("Sample")
        colors = {
            key: tuple(x / 255 for x in value[:-1]) + (0.7,)
            for key, value in zip(unique_types, palette)
        }
        if "sample" in list(inj_metadata[types].unique()):
            colors["sample"] = (1, 1, 1, 0)
        if "Sample" in list(inj_metadata[types].unique()):
            colors["Sample"] = (1, 1, 1, 0)
    row_index = 2
    if (count % 6 == 0 or count % 6 > 3) and (not is_aggregate or row.name == 0):
        row_index = 1
        plt.tight_layout()
        pdf.savefig()
        plt.close("all")
        plt.figure(figsize=(20, 12))

    if not is_aggregate:
        plt.subplot(2, 1, row_index)
        if types == "injection_type":
            plt.title(
                ("Internal Standards Line Graph for: ")
                + (f"{row['Compound_ID']} - {ann_id}")
            )
            plt.ylabel("Raw Abundance")
        else:
            plt.title("Sample Abundance Median Plot")
            plt.ylabel("Median Abundance")

    elif row.name == 0:
        plt.subplot(2, 1, row_index)
        plt.title("Internal Standards Aggregate Line Graph")
        plt.ylabel("Normalized Abundance")

    plt.xlabel("Injection Order")

    if types == "injection_type":
        # remove internal standard name from values
        row_vals = row.values[1:]
    else:
        row_vals = row.values
    plt.scatter(
        inj_metadata["injection_order"],
        row_vals,
        c=inj_metadata[types].map(colors),
        zorder=2,
    )

    if "sample" in colors:
        del colors["sample"]
    if "Sample" in colors:
        del colors["Sample"]
    if not is_aggregate or (is_aggregate and row.name == 0):
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                color="grey",
                markerfacecolor=v,
                label=k,
                markersize=8,
            )
            for k, v in colors.items()
        ]
        if types != "sample_type":
            leg_loc = 3 if is_aggregate else 0
            leg = plt.legend(title="color", handles=handles, loc=leg_loc)
            plt.gca().add_artist(leg)

    plt.plot(inj_metadata["injection_order"], row_vals, c=line_color, zorder=1)
    if not is_aggregate:
        plt.ylim(bottom=0)
    if types == "sample_type":
        handle_median_plot_labels(row, inj_metadata, sample_names, handles)


def handle_median_plot_labels(values, inj_metadata, sample_names, leg_handles):
    """Handles outliers and column break labels on Formation sample median plot"""
    combined_median_df = pd.concat(
        [
            inj_metadata.reset_index(drop=True)["injection_order"],
            values.reset_index(drop=True),
            sample_names,
        ],
        axis=1,
    )
    combined_median_df.columns = [
        "injection_order",
        "median_intensities",
        "sample_names",
    ]
    max_features = combined_median_df.nlargest(10, "median_intensities", keep="all")
    min_features = combined_median_df.nsmallest(10, "median_intensities", keep="all")
    min_max_features = pd.concat(
        [min_features, max_features], ignore_index=True, axis=0
    )
    np.random.seed(0)
    plt.scatter(
        min_max_features["injection_order"],
        min_max_features["median_intensities"],
        marker="*",
        zorder=3,
        color=(0, 0, 0, 0.5),
    )
    min_max_features.apply(
        lambda row: plt.annotate(
            row["sample_names"],
            xy=(row["injection_order"], row["median_intensities"]),
            xycoords="data",
            xytext=(np.random.randn() * 3, np.random.randn() * 5),
            textcoords="offset points",
            fontsize=5,
            ha="left",
            rotation=30,
        ),
        axis=1,
    )
    handles = leg_handles + [
        Line2D(
            [0],
            [0],
            marker="*",
            color="grey",
            markerfacecolor=(0, 0, 0, 0.5),
            label="10 Min and Max Samples",
            markersize=8,
        )
    ]
    leg = plt.legend(title="color", handles=handles)
    plt.gca().add_artist(leg)
    plt.ylim(bottom=0, top=1.05 * values.max())
    col_change_vals = inj_metadata.reset_index()["column_number"].diff()
    col_change = inj_metadata.reset_index()[
        (col_change_vals != 0) & (~col_change_vals.isna())
    ]
    col_change_x = col_change["injection_order"]
    plt.vlines(
        x=col_change_x,
        ymin=0,
        ymax=1.05 * values.max(),
        color="blue",
        linestyle="dashed",
    )


def plot_formation_quartile(data, count, sample_colors, colors, pdf, sample_names):
    """Creates formation sample boxplot"""
    row_index = 2
    if count % 6 == 0 or count % 6 > 3:
        row_index = 1
        plt.tight_layout()
        pdf.savefig()
        plt.close("all")
        plt.figure(figsize=(20, 12))
    axes = plt.subplot(2, 1, row_index)
    axes.set_rasterized(True)
    plt.title("Sample Feature Abundances (1/2 min imputed)")
    plt.xlabel("Sample")
    plt.ylabel("Log(Abundances)")

    medianprops = {"linestyle": "-.", "linewidth": 2.5, "color": "black"}
    bplot = plt.boxplot(data, widths=0.7, medianprops=medianprops, patch_artist=True)
    plt.xticks(range(1, len(data.columns) + 1), sample_names, rotation=90)
    for patch, color in zip(bplot["boxes"], colors):
        patch.set_facecolor(color)
    handles = [mpatches.Patch(color=v, label=k) for k, v in sample_colors.items()]
    axes.legend(
        title="sample_type", handles=handles, loc="center left", bbox_to_anchor=(1, 0.5)
    )


def create_pools_cv_table(feature_metadata, pools, pdf):
    """Creates formation Pools CV Table"""
    ann_features = feature_metadata.loc[~feature_metadata["Annotation_ID"].isnull()]
    bounds = [
        (0, 0.05),
        (0.05, 0.1),
        (0.1, 0.15),
        (0.15, 0.20),
        (0.20, 0.25),
        (0.25, 0.30),
        (0.30, 0.35),
        (0.35, 0.4),
        (0.4, ""),
    ]
    for i, pool in enumerate(pools):
        col_name = f"{pool} CVs"
        if col_name not in feature_metadata.columns:
            continue

        ann_features_col = ann_features[col_name].replace("", float("nan")).dropna()
        feature_met_col = feature_metadata[col_name].replace("", float("nan")).dropna()
        try:
            ann_features_col = ann_features_col.astype("float64")
            feature_met_col = feature_met_col.astype("float64")
        except ValueError:
            ann_indexer = ann_features_col.str.endswith("%")
            all_indexer = feature_met_col.str.endswith("%")
            ann_features_col.loc[ann_indexer] = ann_features_col.loc[
                ann_indexer
            ].str.strip("%")
            feature_met_col.loc[all_indexer] = feature_met_col.loc[
                all_indexer
            ].str.strip("%")
            ann_features_col = ann_features_col.astype("float64")
            feature_met_col = feature_met_col.astype("float64")
            ann_features_col.loc[ann_indexer] /= 100
            feature_met_col.loc[all_indexer] /= 100

        cv_vals = []
        rows_names = []
        for bound in bounds:
            if bound[1] == "":
                known_cvs = ann_features_col >= bound[0]
                all_cvs = feature_met_col >= bound[0]
                row_label = f">{round(bound[0] * 100)}%"
            else:
                known_cvs = (ann_features_col >= bound[0]) & (
                    ann_features_col < bound[1]
                )
                all_cvs = (feature_met_col >= bound[0]) & (feature_met_col < bound[1])
                row_label = f"{round(bound[0] * 100)}% - {round(bound[1] * 100)}%"
            known_count = known_cvs.sum()
            all_count = all_cvs.sum()

            cv_vals.append([known_count, all_count])
            rows_names.append(row_label)

        cv_vals.append(
            [
                ann_features[col_name].shape[0],
                feature_metadata[col_name].shape[0],
            ]
        )
        rows_names.append("Total")
        if i == 0:
            plt.figure(figsize=(20, 12))
        elif i % 4 == 0:
            plt.tight_layout()
            pdf.savefig()
            plt.close("all")
            plt.figure(figsize=(20, 12))
        plt.subplot(2, 2, i % 4 + 1)
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.box(on=None)
        plt.title(col_name)
        plt.table(
            cellText=cv_vals,
            rowLabels=rows_names,
            colLabels=["Knowns", "Knowns + Unknowns"],
            loc="center",
            cellLoc="center",
        )


def _sort_dataset(final_dataset):
    """Sort annotated metabolites by order, then, class MZ, RT;
    sort nontargeted by RT, MZ"""
    # force internal standards sort to top
    final_dataset.loc[
        final_dataset["HMDB_ID"].str.lower() == "internal standard", "superClass"
    ] = "AAAAAAAAAAAAAA"

    # make sure all annotated features are above nontargeted features
    final_dataset = final_dataset.sort_values(by=["Metabolite"])

    annotated = ~pd.isnull(final_dataset["Metabolite"])
    # sort annotated
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
            ]
        )
        .values
    )

    # sort nontargeted separately
    final_dataset.loc[~annotated, :] = (
        final_dataset.loc[~annotated, :].sort_values(by=["RT", "MZ"]).values
    )

    return final_dataset


def _combine_injection_metadata(data, sampleset, injectionset):
    """Generates combined Sampleset and Injectionset metadata with columns matching
    data. Also removes not_used. If there are samples in data not present in either
    sampleset or injectionset, they are preserved and missing metadata is null.
    "not_used" samples are removed as well."""
    warnings = []
    # rename column for older datasets
    if "reporting_name" not in injectionset.columns:
        injectionset = injectionset.rename(columns={"program_id": "reporting_name"})
    combined = injectionset.loc[data.columns].copy()
    combined["broad_id"] = combined["broad_id"].fillna("")
    sampleset_copy = sampleset.copy()
    sampleset_copy.loc["", :] = np.nan  # blank broad_id means no metadata
    combined = combined.join(sampleset_copy, on="broad_id", lsuffix="_inj")

    additional_meta = [
        col for col in sampleset.columns if col not in ["program_id", "broad_id"]
    ]

    # use injection set reporting_name for pools or missing program_ids if present
    pools = combined.index[combined["broad_id"].str.lower() == "pref"]
    missing = pd.isnull(combined["program_id"])
    if "reporting_name" in combined.columns:
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
    combined.loc[:, additional_meta] = combined.loc[:, additional_meta].fillna("NA")

    not_used = combined["injection_type"].str.startswith("not_used")
    if not_used.sum() > 0:
        not_used = combined.index[not_used]
        data.drop(not_used, axis="columns", inplace=True)
        combined = combined.drop(not_used, axis="index")
        warnings.append(
            "The following injections were marked as not_used and were removed from the"
            f" dataset: {', '.join(not_used)}"
        )

    # until filtering can check or create QCRole column
    if "QCRole" in combined.columns:
        combined.rename(columns={"QCRole": "qcrole"}, inplace=True)

    if "qcrole" not in combined.columns:
        warnings.append(
            "No QCRole column is present in the injection info sheet. Please manually "
            "fill the sample_type metadata row if you need this."
        )

    # re-order
    injection_meta = combined.reindex(
        columns=[
            "date_extracted",
            "date_injected",
            "column_number",
            "injection_order",
            "injection_type",
            "qcrole",
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
    names["qcrole"] = "sample_type"
    injection_meta = injection_meta.rename(columns=names)

    injection_meta.loc[:, "date_extracted"] = pd.to_datetime(
        injection_meta["date_extracted"], errors="ignore", infer_datetime_format=True
    )
    injection_meta.loc[:, "date_injected"] = pd.to_datetime(
        injection_meta["date_injected"], errors="ignore", infer_datetime_format=True
    )

    # format labels for each column (for each row, after returning transposed dataframe)
    formats = {
        "date_extracted": "date",
        "date_injected": "date",
        "column_number": "center",
        "injection_order": "center",
        "injection_type": "center",
        "sample_type": "center",
        "raw_file_name": "default",
    }
    formats.update({names[header]: "center" for header in additional_meta})

    return injection_meta, data, formats, warnings


def _combine_feature_metadata(metadata, annotations):
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
    feature_meta = feature_meta.reset_index(drop=True)
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


def combine(
    data,
    f_mdata,
    inj_mdata,
    sample_mdata,
    form_params,
    annotations=None,
    pool_instructions=None,
    include_report=True,
    results_filename=None,
):
    """Given the above information, returns a formatted dataset,
    an optional report, and a filename"""
    method_name = "DefaultMethod"
    if results_filename is None:
        results_filename = "Results"
    abundance_threshold = 1
    warnings = []

    # check mdata and data have same index
    dataset_index = data.index
    if set(f_mdata.index) != set(data.index):
        warnings.append(
            "Your metadata and data have mismatched features, likely because this"
            " is an old dataset. You might have missing features..."
        )
        dataset_index = f_mdata.index.intersection(data.index).copy()
        data = data.loc[dataset_index, :]
        f_mdata = f_mdata.loc[dataset_index, :]

    # sort and drop index; Compound_ID should be in feature metadata
    f_mdata = f_mdata.loc[dataset_index, :].reset_index(drop=True)
    data.reset_index(drop=True, inplace=True)
    data[data < abundance_threshold] = np.nan

    # fill non_quants and cast as object
    if "Non_Quant" not in f_mdata.columns:
        f_mdata["Non_Quant"] = False
    f_mdata["Non_Quant"].fillna(False, inplace=True)
    f_mdata["Non_Quant"] = f_mdata["Non_Quant"].astype("object")

    injection_meta, data, row_formats, warnings = _combine_injection_metadata(
        data, sample_mdata, inj_mdata
    )
    feature_meta = _combine_feature_metadata(f_mdata, annotations)
    if "Method" in feature_meta:
        method_name = feature_meta["Method"][0]

    # identify pooled references, just assume two for now
    pref_as = injection_meta.index[injection_meta["injection_type"] == "prefa"]
    pref_bs = injection_meta.index[injection_meta["injection_type"] == "prefb"]
    miss_column = f"PREF Missing (of {len(pref_as) + len(pref_bs)})"

    column_formats = {
        "superClass": "default",
        "mainClass": "default",
        "subClass": "default",
        "PREFA CVs": "percent",
        "PREFB CVs": "percent",
        miss_column: "default",
        "Annotation_ID": "default",
        "Adduct_Priority": "default",
        "Adduct": "default",
        "MZ_Calculated": "float4",
        "Non_Quant": "default",
        "Cluster_Num": "default",
        "Cluster_Size": "default",
        "Corr_To_Primary": "float4",
        "Primary": "default",
        "Method": "default",
        "Compound_ID": "default",
        "MZ": "float4",
        "RT": "float2",
        "HMDB_ID": "default",
        "HMDB_specificity (1=match; 2=representative)": "center",
        "Metabolite": "default",
    }

    if "Cluster_Num" not in feature_meta:
        for key in ["Cluster_Num", "Cluster_Size", "Primary", "Corr_To_Primary"]:
            column_formats.pop(key)

    # apply pool modifications to raw data to match drift corrected results
    if not form_params["drift_corrected"] and pool_instructions is not None:
        data_index = data.index.copy()
        data.index = feature_meta["Compound_ID"]
        data = apply_pool_instructions(data, pool_instructions, inplace=True)
        data.index = data_index

    if len(pref_as) > 0:
        feature_meta["PREFA CVs"] = data.loc[:, pref_as].std(axis=1) / data.loc[
            :, pref_as
        ].mean(axis=1)
    else:
        feature_meta["PREFA CVs"] = None
    if len(pref_bs) > 0:
        feature_meta["PREFB CVs"] = data.loc[:, pref_bs].std(axis=1) / data.loc[
            :, pref_bs
        ].mean(axis=1)
    else:
        feature_meta["PREFB CVs"] = None
    feature_meta[miss_column] = data.loc[:, pref_as.union(pref_bs)].isnull().sum(axis=1)

    # filter, if requested
    if form_params["filter_features"]:
        if form_params["missing_as_percent"]:
            missing_cutoff = np.ceil(
                (len(pref_as) + len(pref_bs)) * form_params["missing_cutoff"] / 100
            )
        else:
            missing_cutoff = form_params["missing_cutoff"]
        to_filter = feature_meta[miss_column] > missing_cutoff

        to_filter = to_filter | (
            feature_meta["PREFA CVs"] > form_params["cv_cutoff"] / 100
        )
        to_filter = to_filter | (
            feature_meta["PREFB CVs"] > form_params["cv_cutoff"] / 100
        )

        if "Primary" in feature_meta.columns:
            to_filter = to_filter | ~feature_meta["Primary"].fillna(False)

        # annotated compounds don't get filtered by CV or missing PREFs
        to_filter = to_filter & pd.isnull(feature_meta["Annotation_ID"])

        to_filter = to_filter | feature_meta["Non_Quant"].fillna(False)
        # find all duplicated annotations that weren't filtered as non_quant
        annotated = ~to_filter & ~pd.isnull(feature_meta["Annotation_ID"])
        is_duplicated = feature_meta.loc[annotated, "Annotation_ID"].duplicated(
            keep=False
        )
        duplicates = feature_meta.loc[is_duplicated.index[is_duplicated]]

        if "Adduct" in feature_meta.columns:
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
                        to_filter[to_drop_idx] = True
                        break

            # find annotations that are still duplicated
            annotated = ~to_filter & ~pd.isnull(feature_meta["Annotation_ID"])
            is_duplicated = feature_meta.loc[annotated, "Annotation_ID"].duplicated(
                keep=False
            )
            duplicates = feature_meta.loc[is_duplicated.index[is_duplicated]]

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
                to_filter[to_drop_idx] = True

        data = data.loc[~to_filter]
        feature_meta = feature_meta.loc[~to_filter]

        # finally, if there are still duplicates, warn the user
        annotated = ~to_filter & ~pd.isnull(feature_meta["Annotation_ID"])
        is_duplicated = feature_meta.loc[annotated, "Annotation_ID"].duplicated(
            keep=False
        )
        if sum(is_duplicated) > 0:
            warnings.append(
                "There are duplicate annotations that could not be filtered by adduct "
                + "or extraction method. Please remove duplicate annotations before "
                + "sharing results."
            )

    # strip lone _As, _Bs, and _Cs
    grouped_annotations = feature_meta.loc[
        ~pd.isnull(feature_meta["Metabolite"])
        & feature_meta["Metabolite"].str.endswith(("_A", "_B", "_C")),
        "Metabolite",
    ]
    for i in grouped_annotations.index:
        if (grouped_annotations.str[:-2] == grouped_annotations[i][:-2]).sum() == 1:
            feature_meta.loc[i, "Metabolite"] = feature_meta.loc[i, "Metabolite"][:-2]

    # normalize
    if form_params["to_normalize"]:
        to_norm_cols = injection_meta.index[
            injection_meta["injection_type"].str.startswith("pref").values
            | injection_meta["injection_type"].str.startswith("sample").values
        ].intersection(data.columns)
        to_norm = data.loc[:, to_norm_cols]
        medians = to_norm.median()
        to_norm = to_norm.mul(((1 / medians) * medians.mean()), axis="columns")
        data[to_norm_cols] = to_norm[to_norm_cols]

    # fill values less than 1 with 1 and round
    data = data.clip(lower=1).fillna(0)
    data = data.round(0).astype(np.int64)
    data = data.replace(0, np.nan)

    # filter out any prefs that mess with the PCAs
    if "prefs_to_remove" in form_params:
        data = data.drop(columns=form_params["prefs_to_remove"])
        injection_meta = injection_meta.drop(form_params["prefs_to_remove"])

    # keep additional columns but drop columns that are for internal use only
    additional_cols = [
        col
        for col in feature_meta
        if col not in column_formats and not col.startswith("__")
    ]
    feature_meta = feature_meta.reindex(columns=additional_cols + list(column_formats))

    # run report
    report_pdf = None
    if include_report:
        report_pdf = report(
            data,
            injection_meta,
            feature_meta,
            dataset_name=results_filename,
            write_pdf=False,
        )

    final_dataset = pd.concat([feature_meta, data], axis="columns")
    final_dataset = _sort_dataset(final_dataset)

    # set headers as top row
    injection_meta = injection_meta.T
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

    if form_params["csv_fallback"]:
        return final_dataset, report_pdf, warnings

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
        header: formats[label] for header, label in column_formats.items()
    }
    row_formats = {index: formats[label] for index, label in row_formats.items()}

    # column widths
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("Compound_ID"),
        final_dataset.columns.get_loc("Compound_ID"),
        86,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("superClass"),
        final_dataset.columns.get_loc("superClass"),
        80,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("mainClass"),
        final_dataset.columns.get_loc("mainClass"),
        80,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("subClass"),
        final_dataset.columns.get_loc("subClass"),
        80,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc(miss_column),
        final_dataset.columns.get_loc(miss_column),
        80,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("Annotation_ID"),
        final_dataset.columns.get_loc("Annotation_ID"),
        90,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("Non_Quant"),
        final_dataset.columns.get_loc("Non_Quant"),
        80,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("HMDB_ID"),
        final_dataset.columns.get_loc("HMDB_ID"),
        102,
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("HMDB_specificity (1=match; 2=representative)"),
        final_dataset.columns.get_loc("HMDB_specificity (1=match; 2=representative)"),
        270,
    )
    worksheet.set_column(
        final_dataset.columns.get_loc("Metabolite"),
        final_dataset.columns.get_loc("Metabolite"),
        final_dataset["Metabolite"].str.len().max(),
    )
    worksheet.set_column_pixels(
        final_dataset.columns.get_loc("Metabolite") + 1,
        len(final_dataset.columns) - 1,
        150,
    )

    for i, row in enumerate(final_dataset.values):
        index_label = final_dataset.index[i]
        if index_label in row_formats.keys():  # injection metadata
            worksheet.set_row_pixels(i, 16, row_formats[index_label])
        else:  # default to int
            worksheet.set_row_pixels(i, 16, formats["int"])

        for j, val in enumerate(row):
            column_label = final_dataset.columns[j]
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
    return output, report_pdf, warnings
