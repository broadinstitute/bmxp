"""Module from __init__.py"""
from bmxp.eclipse import MSAligner


datasets = [
    "../../tests/test1.csv",
    "../../tests/test2.csv",
    "../../tests/example1.csv",
]
names = ["DS1", "DS2", "DS3"]
############
# Example 1 - Unsupervised alignment between all datasets
############

a = MSAligner(*datasets, names=names)
a.align()
# Automatically aligns all datasets to all other datasets with defaults

print(a.anchors["DS1"])
# Output list of features
# ['0', '1', '10', '100', '101', ...]

print(a.coarse_matches["DS1"]["DS2"])
# Output the coarse alignment results from DS1 -> DS2.
# Will contain duplicates.
#      DS1  DS2
# 0      0    0
# 1      1    1
# 2     10   10
# 3    100  100
# 4    101  101


print(a.scalers["DS1"]["DS2"])
# Output dict of scalers. For example:
# [271 rows x 2 columns]
# {'RT':         x         y
# 0    1.68  0.004317
# 1    3.53 -0.008928
# ..    ...       ...
# 269  5.51 -0.021388
# 270  5.45 -0.020587
#
# [271 rows x 2 columns], 'MZ':             x         y
# 0    121.1013  0.828166
# 1    197.0667  1.028160
# ..        ...       ...
# 269  792.5513  1.144283
# 270  776.5563  1.152530
#
# [271 rows x 2 columns], 'Intensity':             x         y
# 0    4.703667 -0.097755
# 1    4.341406 -0.167542
# ..        ...       ...
# 269  5.963024 -0.089592
# 270  6.352471 -0.094823


print(a.stds["DS1"]["DS2"])
# Output standard deviations to DS2
# {
#   'RT': 0.01441794266731159,
#   'MZ': 0.12278117619894874,
#   'Intensity': 0.1113562600674442
# }


print(a.scaled_values["DS1"]["DS2"])
# Output DataFrame of scaled values
#            RT          MZ     Intensity
# 0    1.684317  121.101400  4.035640e+04
# 1    3.521072  197.066903  1.492328e+04
# ..        ...         ...           ...
# 297  7.319298  793.542507  1.170663e+06
# 298  0.345388  158.026900  1.156930e+06


print(a.matches["DS1"]["DS2"])
# Output a DataFrame of DS1 feature names to top n (1? 2?) DS2 matches
# [299 rows x 3 columns]
#         0
# 0       0
# 1    None
# 2       2
# ..    ...
# 296  None
# 297   268
# 298  None


# Saves a report of the alignment params
a.report(filepath="results.pdf")

# Reports only matches in all 3 datasets
a.to_csv(filepath="full_matches.csv")

# Reports all cliques with a match to dataset 1
# (Will likely have multiple matches)
a.to_csv(
    "DS1", filepath="DS1.csv", g_size_or_loss=1, c_size_or_loss=1, remove_rerank=False
)


############
# Example 2 - Supervised all by all alignment with custom initial params
############
a = MSAligner(*datasets[:2], names=["DS1", "DS2"])

# Change the initial search window
# Automatically applies the opposite to DS2->DS1
a.set_coarse_params({"DS1": {"DS2": {"RT": {"upper": 0.5, "lower": -0.01}}}}, rec=True)

# Automatically generates all anchors for every dataset
a.gen_anchors()

# Generate coarse matches
a.gen_coarse()

# Generate the scalers. Will not overwrite custom scalers
a.gen_scalers()

# Find standard deviations
a.gen_stds()

# Scale the feature descriptors
a.gen_scaled_values()

# Generate subalignment matches
a.gen_matches()

# generates a.graph, a network-x graph
a.gen_graph()


############
# Example 3 - One way supervised alignment to another dataset
############
a = MSAligner(*datasets, names=names)

# Automatically generates all anchors for every dataset
a.gen_anchors()

# Creates a coarse result of anchors from DS1 to DS2
a.gen_coarse("DS1", "DS2")

# Create the scalers for DS1 to DS2
a.gen_scalers("DS1", "DS2")

# Calculate the standard deviations
a.gen_stds("DS1", "DS2")

# Scale the feature descriptors
a.gen_scaled_values("DS1", "DS2")

# Scales DS1 to DS2 based on generated scalers
a.gen_scaled_values("DS1", "DS2")

a.gen_matches("DS1", "DS2")

a.report(
    datasets_1="DS1",
    datasets_2="DS2",
    filepath="one-way-results.pdf",
)


############
# Example 4 - Two way unsupervised alignment to one or more other datasets
############
a = MSAligner(*datasets, names=["DS1", "DS2", "DS3"])

# Align the
a.align("DS1", ["DS2", "DS3"])
a.align(["DS2", "DS3"], "DS1")

a.report(datasets_1="DS1", datasets_2=["DS2", "DS3"], filepath="DS1_to.pdf")
a.report(
    datasets_1=["DS2", "DS3"],
    datasets_2="DS1",
    filepath="DS1_from.pdf",
)

a.to_csv("DS1", g_size_or_loss=1, c_size_or_loss=1, remove_rerank=False)
