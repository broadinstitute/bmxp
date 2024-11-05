"""
A quick example on how to use Gravity
"""

import pandas as pd
from bmxp.gravity import cluster

# load a dataset
df = pd.read_csv("../../tests/test_gravity.csv", header=0, index_col="Compound_ID")

# cluster, with an RT cutoff of 0.025 and required correlation of .9
sample_df = cluster(df, 0.025, 0.9)

# returns a dataframe with Cluster_Num and Cluster_Size
print(df)

# add to our own dataframe
df["Cluster_Num"] = sample_df["Cluster_Num"]
df["Cluster_Size"] = sample_df["Cluster_Size"]

# record results
df.to_csv("results.csv")
