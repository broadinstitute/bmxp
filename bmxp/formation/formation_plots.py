import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches


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
            xy=(pc_df["pc1"].iloc[i], pc_df["pc2"].iloc[i]),
            xycoords="data",
            xytext=(pc1_jitter.iloc[i], pc2_jitter.iloc[i]),
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
            sample_names.reset_index(drop=True),
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
    col_vals = inj_metadata.reset_index()["column_number"]
    col_change_vals = col_vals == col_vals.shift().fillna(col_vals)
    col_change = inj_metadata.reset_index()[~col_change_vals]
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
