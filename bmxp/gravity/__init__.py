"""
Gravity clusters LCMS datasets by RT and Correlation
"""
import logging
import math
import pandas as pd
import networkx as nx

logging.basicConfig()
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)


__version__ = "0.0.1"


def cluster(df, rt_thresh=0.02, corr_value=0.8, batch_size=1000, method="Spearman"):
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
    Returns
    Dataframe, containing cluster number and number of members for each feature. -1
       indicates a single, unclustered feature.
    """
    eclipse_columns = [
        "RT",
        "MZ",
        "Intensity",
        "Non_Quant",
        "Compound_ID",
        "Adduct",
        "Annotation_ID",
        "Metabolite",
    ]
    sample_df = df.copy().drop(eclipse_columns, axis=1, errors="ignore").transpose()
    sample_df.columns.name = None  # otherwise it crashes during stack
    rt_series = df["RT"]
    num_batches = math.ceil(len(sample_df.columns) / batch_size)
    graph = nx.Graph()

    for i in range(num_batches):
        LOGGER.info(f"Correlating Batch {i+1}/{num_batches}...")
        for j in range(i, num_batches):
            # calculate correlation of the two batches
            batch_a = sample_df.iloc[:, i * batch_size : (i + 1) * batch_size]
            batch_b = sample_df.iloc[:, j * batch_size : (j + 1) * batch_size]
            if method == "Spearman":
                batch_a = batch_a.rank(na_option="top")
                batch_b = batch_b.rank(na_option="top")
            a_zs = batch_a - batch_a.mean()
            b_zs = batch_b - batch_b.mean()
            corr = (
                a_zs.T.dot(b_zs)
                .div(len(batch_a))
                .div(b_zs.std(ddof=0))
                .div(a_zs.std(ddof=0), axis=0)
            )
            links = corr.stack().reset_index()
            links.columns = ["var1", "var2", "value"]
            # filter spearman
            links = links.loc[
                (links["value"] > corr_value) & (links["var1"] != links["var2"])
            ]
            # get RT and filter by that, and append to list
            links["RT_diff"] = (
                rt_series[links["var1"]].values - rt_series[links["var2"]].values
            )
            links["RT_diff"] = links["RT_diff"].abs()
            links = links.loc[(links["RT_diff"] < rt_thresh)]
            edges = [
                (s, t, {"coor": coor})
                for s, t, coor in zip(
                    links["var1"].values, links["var2"].values, links["value"]
                )
            ]
            graph.add_nodes_from(links["var1"].values)
            graph.add_nodes_from(links["var2"].values)
            graph.add_edges_from(edges)

    cluster_df = pd.DataFrame({"Cluster_Num": None, "Cluster_Size": 1}, index=df.index)
    i = 0
    while True:
        # delete singletons
        graph.remove_nodes_from(list(nx.isolates(graph)))

        # find the largest cliques
        max_cliques = []
        max_size = 1
        for clique in nx.find_cliques(graph):
            if len(clique) > max_size:
                max_size = len(clique)
                max_cliques = [clique]
            elif len(clique) == max_size:
                max_cliques.append(clique)
        if not max_cliques:
            break
        # calculate correlation
        best_clique = None
        best_score = 0
        for clique in max_cliques:
            temp = 0
            clique_graph = graph.subgraph(clique).copy()
            for edge in clique_graph.edges(data=True):
                temp += edge[2]["coor"]
            if temp > best_score:
                best_score = temp
                best_clique = clique

        # label, delete and iterate
        features = list(best_clique)
        cluster_df.loc[features, "Cluster_Num"] = i
        cluster_df.loc[features, "Cluster_Size"] = max_size
        if i % 50 == 0:
            LOGGER.info(f"On cluster {i}, which contains {max_size} features.")
        i += 1
        graph.remove_nodes_from(features)
    return cluster_df
