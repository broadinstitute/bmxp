"""
Gravity clusters LCMS datasets by RT and Correlation


NaN Policy - Zeroes, Drop, Fill(spearman)

Method, NaN

Options:
Spearman - Fill NaNs with Zeroes (Fill in Python, call by_corr with ties=whatever)
Spearman - Drop NaNs (Don't fill, call by_rt with Ties=whatever)
Spearman - Backfill NaNs (Don't fill, call "backfill")

Pearson - Fill NaNs with Zeroes (Fill in Python, call by_corr)
Pearson - Drop NaNs (Don't fill, call by_rt)


"""

import logging
import math
from ctypes import CDLL, c_double, c_int32, c_void_p
import os
import platform
import pandas as pd
import numpy as np
import networkx as nx
from bmxp.gravity import fallback
from bmxp import FMDATA

dir_path = os.path.dirname(os.path.realpath(__file__))
correlation = None
if platform.system() == "Windows" and os.path.exists(
    os.path.join(dir_path, "correlation.dll")
):
    correlation = CDLL(os.path.join(dir_path, "correlation.dll"))
elif platform.system() == "Linux" and os.path.exists(
    os.path.join(dir_path, "correlation.so")
):
    correlation = CDLL(os.path.join(dir_path, "correlation.so"))

c_array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
int32_array = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")

if correlation:
    correlation.free_p.argtypes = [c_void_p]

    correlation.pearson_array.argtypes = [c_array, c_int32, c_int32, c_int32]
    correlation.pearson_array.restype = c_double

    correlation.spearman_array.argtypes = [
        c_array,  # arr
        c_int32,  # r1_start
        c_int32,  # r2_start
        c_int32,  # size
        c_int32,  # dropNan
        c_int32,  # legacyMode
    ]
    correlation.spearman_array.restype = c_double
    correlation.pearson.argtypes = [c_array, c_array, c_int32]
    correlation.pearson.restype = c_double

    correlation.spearman.argtypes = [c_array, c_array, c_int32, c_int32, c_int32]
    correlation.spearman.restype = c_double
else:
    correlation = fallback

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

__version__ = "0.2.0"


def free_p(p):
    correlation.free_p(p)


def pearson_array(arr1, r1_start, r2_start):
    r_p = correlation.pearson_array(
        arr1,
        r1_start,
        r2_start,
        len(arr1[0]),
    )
    return r_p


def spearman_array(arr1, r1_start, r2_start, nan_policy="fill", legacy_mode=False):
    drop_nan = nan_policy != "fill"
    r_p = correlation.spearman_array(
        arr1, r1_start, r2_start, len(arr1[0]), drop_nan, legacy_mode
    )
    return r_p


def pearson(x, y):
    x = np.array(x).astype(np.float64)
    y = np.array(y).astype(np.float64)
    return correlation.pearson(x, y, len(x))


def spearman(x, y, nan_policy="fill", legacy_mode=False):
    x = np.array(x).astype(np.float64)
    y = np.array(y).astype(np.float64)
    drop_nan = nan_policy != "fill"
    return correlation.spearman(x, y, len(x), drop_nan, legacy_mode)


def deconvolute(df_index, graph, fmdata):
    num_label = fmdata["Cluster_Num"]
    size_label = fmdata["Cluster_Size"]
    cluster_df = pd.DataFrame({num_label: None, size_label: 1}, index=df_index)
    i = 0
    while True:
        # delete singletons
        graph.remove_nodes_from(list(nx.isolates(graph)))

        # build list of cliques of the largest size
        cliques = []
        max_size = 0
        for clique in nx.find_cliques(graph):
            if len(clique) > max_size:
                max_size = len(clique)
                cliques = [[set(clique), 0]]
            elif len(clique) == max_size:
                cliques.append([set(clique), 0])
        if not cliques:
            break
        # calculate correlation of largest cliques
        for clique in cliques:
            clique_graph = graph.subgraph(clique[0])
            temp = 0
            for edge in clique_graph.edges(data=True):
                temp += edge[2]["coor"]
            clique[1] = temp

        # sort by best corr
        cliques.sort(key=lambda x: x[1], reverse=True)
        to_remove = set()
        while len(cliques) > 0:
            # add the best one
            best_clique = cliques[0]
            features = list(best_clique[0])
            cluster_df.loc[features, num_label] = i
            cluster_df.loc[features, size_label] = max_size
            if i % 50 == 0:
                LOGGER.info(f"On cluster {i}, which contains {max_size} features.")
            i += 1

            # search the clique for ones to remove
            cliques = [
                clique
                for clique in cliques
                if not clique[0].intersection(best_clique[0])
            ]
            # add to our "to remove" queue
            to_remove.update(best_clique[0])
        graph.remove_nodes_from(to_remove)
        if max_size == 2:
            break
    return cluster_df


def corr_array(
    i, j, rt_series, df, corr_value, rt_thresh, method, nan_policy, legacy_mode
):
    # long winded way to find RT differences and the indices
    # similar to pd.stack
    np_batch_a = rt_series.iloc[i[0] : i[1]].values
    np_batch_b = rt_series.iloc[j[0] : j[1]].values
    rt_flattened = np.absolute(np.subtract.outer(np_batch_b, np_batch_a)).flatten()
    a_index = np.tile(np.fromiter(range(i[0], i[1]), int), (j[1] - j[0]))
    b_index = np.tile(np.fromiter(range(j[0], j[1]), int), ((i[1] - i[0]), 1)).flatten(
        order="F"
    )
    to_keep = (rt_flattened < rt_thresh) & (a_index != b_index)
    a_index = a_index[to_keep]
    b_index = b_index[to_keep]
    corr_results = np.empty(len(a_index))
    for k, a_val in enumerate(a_index):
        b_val = b_index[k]
        if method == "spearman":
            corr_results[k] = spearman_array(df, a_val, b_val, nan_policy, legacy_mode)
        else:
            corr_results[k] = pearson_array(df, a_val, b_val)
    to_return = corr_results > corr_value
    return a_index[to_return], b_index[to_return], corr_results[to_return]


def gen_batches(num_features, batch_size):
    """
    Generates batches for clustering, based on number of features and batch size
    """
    num_batches = math.ceil(num_features / batch_size)
    for i in range(num_batches):
        for j in range(i, num_batches):
            i_indices = [
                i * batch_size,
                min((i + 1) * batch_size, num_features),
            ]
            j_indices = [
                j * batch_size,
                min((j + 1) * batch_size, num_features),
            ]
            yield {
                "i": i,
                "j": j,
                "i_indices": i_indices,
                "j_indices": j_indices,
                "num_batches": num_batches,
            }


def cluster(
    df,
    rt_thresh=0.02,
    corr_value=0.8,
    batch_size=1000,
    method="spearman",
    nan_policy="fill",
    legacy_mode=False,
    schema_labels=None,
):
    """
    Cluster aggregates LCMS features into groups based on sample-correlation and
        retention time. It builds a network based on retention time difference and
        correlation. Then it labels the largest clique as a cluster, deletes the
        features from the network, and identifies the next largest clique. Ties are
        broken but a summation of correlations in the cluster.

    df: Dataframe, Eclipse compatible
    rt_thresh: float, maximum Rt different for clustering
    corr_value: float, correlation cutoff
    batch_size: int, batch size for calculating correlations
    method: string, "Spearman" or "Pearson" correlation
    legacy_mode: boolean, when True, NaNs are filled with 0s when correlating features
    with very few overlapping non-NaN values (Spearman only)
    Returns
    Dataframe, containing cluster number and number of members for each feature. -1
       indicates a single, unclustered feature.
    """
    fmdata = FMDATA.copy()
    if schema_labels is not None:
        fmdata.update(schema_labels)
    rt_series = df[fmdata["RT"]]
    rt_series.index.name = None  # otherwise it crashes during stack
    sample_df = df.copy().drop(list(fmdata.keys()), axis=1, errors="ignore")

    num_features = len(sample_df.index)
    graph = nx.Graph()
    sample_df = sample_df.values.copy()
    if nan_policy == "zeroes":
        sample_df = np.nan_to_num(sample_df)
    for batch_info in gen_batches(num_features, batch_size):
        num_batches = batch_info["num_batches"]
        LOGGER.info(f"Correlating Batch { batch_info['i'] + 1}/{num_batches}...")
        a_index, b_index, corrs = corr_array(
            batch_info["i_indices"],
            batch_info["j_indices"],
            rt_series,
            sample_df,
            corr_value,
            rt_thresh,
            method,
            nan_policy,
            legacy_mode,
        )
        a_index = df.index[a_index]
        b_index = df.index[b_index]
        edges = [(s, t, {"coor": coor}) for s, t, coor in zip(a_index, b_index, corrs)]
        graph.add_nodes_from(a_index)
        graph.add_nodes_from(b_index)
        graph.add_edges_from(edges)
    return deconvolute(df.index, graph, fmdata)


if __name__ == "__main__":
    data = pd.read_csv(
        r"C:\Users\danie\work\beta-muricholate.csv", index_col="Compound_ID"
    )
    data = data.replace(0, np.nan)
    # fmt: off
    arr1 = [
        [np.nan, np.nan],  # 0
        [1, np.nan],  # 1
        [1, np.nan],
        [1, np.nan],
        [1, np.nan],
        [2, np.nan],  # 5
        [2, np.nan],
        [2, 1],  # 7
        [np.nan, 2],
        [3, np.nan],
        [3, np.nan],  # 10
        [3, 5],
        [np.nan, np.nan],
        [4, 7],
        [4, np.nan],
        [4, 9],  # 15
        [4, 10],  # 16
        [np.nan, np.nan],  # 17
    ]

    arr = np.array(arr1).T.astype(np.float64).copy()
    # fmt: on

    # print(pearson_array_m(arr1, 0, 1))
    # spearman zeroes
    # spearman fill
    # spearman drop
    # pearson zeroes
    # pearson drop
    print(spearman_array(np.nan_to_num(arr), 0, 1, nan_policy="drop"))
    # print(pearson_array(arr1, 0, 1))
    # print("done")
