# MSEclipse

MSEclipse is a library that finds corresponding features in same-method LC-MS datasets. It accomplishes this by breaking the whole alignment own into individual directed subalignments, loading the results into a graph, then deconvoluting this graph to produce a results table.

## Quick use
1. Install with `pip install bmxp` (requires Python 3)
2. Get two (or more!) datasets with the columns `Compound_ID`, `RT`, `MZ`, `Intensity` (provided in `/tests`).
3. Have fun!

```python
from bmxp.eclipse import MSAligner
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3']) # Initialization Phase
a.align() # Subalignment Phase
a.to_csv() # Feature Aggregation Phase
a.report() # Generate a report
```

This will run the data on default settings, and will only produce singular matches (no redundancies). Somtimes, default settings don't work though, so feel free to continue reading.


## Overview
`MSAligner.datasets` : dict[str, DataFrame], a dictionary containing each dataset and its name.

e.g. `print(a.datasets["DS1"]`

Datasets will be loaded into the object and given a name, such as `HP-1`. The datasets should have columns for retention time, mass-to-charge ratio, average feature intensity, and sample values. Also, simplified datasets are created, which are used for anchor generation.

### Initialization Phase

#### Specify Params

Eclipse has two way to specify params. One is through whole-alignment class attributes, using `set_defaults`. The other is through specific subalignments. Only change whole-alignment params if you want this change to be propagated to all suablignments. If you need to change just one or two subalignment pairs, use the subalignment-specific params. For more information on changing the defaults, see the "Changing Defaults" section.

Whole-alignment params, best modified using the `set_defaults` method:
1. `descriptors`: `{"RT": "linear", "MZ": "ppm", "Intensity": "log10"}` - Specified by `descriptors`. The column names and their feature transformation space can be specified. 
2. `default_cutoffs`: `{"RT": 6, "Intensity": 6, "MZ": 6}` - Specified by `cutoffs`. Allows 6 standard deviations of variance when matching
3. `default_weights`: `{"RT": 1, "Intensity": 1, "MZ": 1}` - Specified by `weights`. Equally weighs all descriptors when matching
4. `default_coarse_params`: `{
            "RT": {"upper": 0.5, "lower": -0.5},
            "MZ": {"upper": 15, "lower": -15},
            "Intensity": {"upper": 2, "lower": -2},
        }` - Specified by `coarse_params`. The coarse params for anchor determination and coarse matching.
5. `feature_name`: `Compound_ID` - The name of the rows.

Subalignment-specific params, best modified with `set_params`. You may specify the source and target dataset. By default, those params will be applied to that subalignment as well as the reverse:
1. `coarse_params`:
2. `scalers`
3. `cutoffs`
4. `weights`

#### Datasets are loaded
Any number of datasets can be loaded. Loading over 10 datasets might require changing subalignment from all-to-all to one-to-all.
#### Anchors are determined

`MSAligner.anchors`: dict[str, tuple], a dictionary containing dataset names and a tuple containing the index names of the anchors.


```Python
a.gen_anchors()
print(a.anchors["DS1"])
```
Simplified datasets are used for determining scalers. They are actually calculated during the subalignment phase. There are two modes, where `remove_all` can be True or False:
`True`: Removes all features which fall within the coarse params of eachother.
`False`: Keeps only the most prioritized feature which falls within the coarse params. Uses the `priority` variable, which is currently hardcoded to `Intensity`.

### Subalignment Phase (Scaling, Matching))
Running `.align()` begins the actual process scaling and matching. This conveniently runs all the functions shown below, as well as the `MSAligner.gen_anchors` function above.

The default subalignment queue is and all-by-all alignment. If one-to-all must be run, it requires two runs.

```Python
# all-by-all, default
a.align()

# one-to-all
# a.align('DS1', ['DS2', 'DS3'])
# a.align(['DS2', 'DS3'], 'DS1'])
```

There are several steps, and they can be run individually if desired. Below shows how to run this piecemeal.

#### Coarse Matches (For Scaling) are determined
`MSAligner.coarse_matches`: dict[str, DataFrame], a dictionary containing dataset names and dataframes with matches between anchor datasets.

```Python
a.gen_coarse()
print(a.coarse_matches["DS1"]["DS2"])
```

Coarse matches will be an initial match between two unscaled anchor datasets. At first, the matches are recorded along with their descriptors. There may be duplicates, given the ambiguous nature of this initial match. In a later step, we will scale the feature descriptors for these matches. Each DataFrame will have the columns `index, RT1, MZ1, int1, rt2, mz2, int1, rt1_scaled, mz1_scaled, int1_scaled`. These matches will be used to determine the scalers.

#### Scalers are determined
`MSAligner.scalers`: dict[str, dict[str, str/float]], a dictionary containing scaling constants to every other dataset.

```Python
a.gen_scalers()
print(a.scalers["DS1"]["DS2"])
```
The scalers will be the scaling instructions to perform a non-parametric adjustment. The values will be in the space of the feature -- `linear`, `ppm` or `intensity`. They are determined by comparing the coarse matches between two dataset. To be used in actual scaling tasks, they are linearly interpolated.

As an example, say we have coarse matches between two datasets and it looks like so:
```
rt1: 1.0, 2.0, 3.0, 4.0, 5.0
rt2: 1.1, 2.2, 3.1, 3.8, 5.0
```
In this case, if we used a LOWESS curve, and it fit perfectly, our scalers would be:
```
rtx: 1.0, 2.0, 3.0,  4.0, 5.0
rty: 0.1  0.2, 0.1, -0.2, 0
```
#### Variance (Noise, Standard Deviation) Determination
`MSAligner.stds`: dict[str, dict[str, dict[str, float]]], a dictionary containing the residual standard deviations for the descriptors in a subalignment pair.

```Python
a.gen_stds()
print(a.stds["DS1"]["DS2"])
```

After scaling, there will still be noise on the fit as well as outliers from a wide coarse matching window. Outliers will be removed by a MAD based algorithm, then standard deviations will be calculated. These standard deviations will be used for the actual match step, for both selection as well as scoring.

#### Scaling Full Datasets
`MSAligner.scaled_values`: dict[str, dict[str, DataFrame]], a dictionary containing the scaled descriptors for a subalignment pair. The source dataset (i.e. DS1) is scaled.

```Python
a.gen_scaled_values()
print(a.scaled_values["DS1"]["DS2"])
```
The above will print the scaled values for DS1 when being compared to DS2.

#### Subalignment Feature Matching
`MSAligner.gen_matches`: dict[str, dict[str, Dataframe]], a dictionary containing best matches.

Example:
```Python
a.gen_matchs()
print(a.matches["DS1"]["DS2"])
```

Once the dataset is scaled, we can find best matches to a dataset in which it is being aligned. A window will be calculated (based on standard deviations), and features that fall in that window will be scored (again, based on standard deviations). The top n best matches will be recorded.

### Feature Aggregation Phase

The whole alignment can be aggregated and deconvoluted with a single command:
```Python
# subalignment match report
a.report()

# strict, default, must contain all members
a.to_csv()

# Dataset centric - reports all matches to a specified dataset
# a.to_csv('DS1', union_only=False)
```

#### Network construction
```Python
a.gen_graph()
print(a.graph) # networkx graph
```
 This command builds a directed graph based on the subalignment matches. It can be exported, or deconvoluted

#### Network deconvolution
See the main heading, Feature Aggregation Phase
