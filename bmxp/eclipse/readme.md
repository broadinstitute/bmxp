# MSEclipse

MSEclipse is a library that finds corresponding features in same-method LC-MS datasets. It accomplishes this by breaking the whole alignment own into individual directed subalignments, loading the results into a graph, then deconvoluting this graph to produce a results table.

## Quick use
1. Install with `pip install bmxp` (requires Python 3)
2. Get two (or more!) datasets with the columns `Compound_ID`, `RT`, `MZ`, `Intensity` (provided in `/tests`). *Note: These can be changed either in the BMXP Schema or directly in Eclipse
3. Have fun!

```python
from bmxp.eclipse import MSAligner
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3']) 
a.align() # Subalignment Phase
a.to_csv() # Feature Aggregation Phase
a.report() # Generate a report
```

This will run the data on default settings, and will only produce singular matches (no redundancies). Sometimes, default settings don't work though, so feel free to continue reading.

Don't like our column header names? If you don't want to adjust them globally through the [BMXP Modules](https://github.com/broadinstitute/bmxp/), just do this:
```python
from bmxp.eclipse import MSAligner
schema_labels = {
    "Compound_ID" : "feature_id",
    "RT": "rt (min)",
    "MZ": "m/z",
    "Intensity": "abundance"
}
# be sure all datasets use "feature_id", "rt (min)", etc...
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'], schema_labels=schema_labels) 
a.align() # Subalignment Phase
```

### Overview
`MSAligner.datasets` : dict[str, DataFrame], a dictionary containing each dataset and its name.

e.g. `print(a.datasets["DS1"])`

Datasets will be loaded into the object and given a name, such as `HP-1`. The datasets should have columns for retention time, mass-to-charge ratio, average feature intensity, and sample values. Also, simplified datasets are created, which are used for anchor generation.

### Initialization Phase
#### Specify Params
Sometimes the default values aren't good enough. Eclipse has three ways to modify the alignment params:
1. Global, class-level defaults based on (immutable) schema keys
2. Instance-level defaults based on column names (`schema_labels`), using `set_instance_defaults`
3. Subalignment level params based on column names, using `set_params`, or the individual `set_scalers`, `set_cutoffs`, `set_weights`, `set_coarse_params`, or `set_scaler_params`

When an MSAligner object is initialized, the global class-level defaults are copied into the instance. Then if `set_instance_defaults` is used it will only change the instance. Finally if `set_params` is used to specify subalignments, these are evaluated during `align`.

##### Global, class-level defaults
Global settings must be applied before objects are initialized. It will not affect already initialized MSAligner objects. This can be done by modifying the class attributes, using the schema keys.
The variables are:
1. `MSAligner.descriptors`: `{"RT": "linear", "MZ": "ppm", "Intensity": "log10"}` - Specified by `descriptors`. The column names and their feature transformation space can be specified. 
2. `MSAligner.default_cutoff`: `6` - Allows 6 standard deviations of variance when matching
3. `MSAligner.default_weight`: `1` - Equally weighs all descriptors when matching
4. `MSAligner.default_cutoffs`: `{}` - Overrides the default cutoff for a descriptor
5. `MSAligner.default_weights`: `{}` - Overrides the default weight for a descriptor
6. `MSAligner.default_coarse_params`: `{
            "RT": {"upper": 0.5, "lower": -0.5},
            "MZ": {"upper": 15, "lower": -15},
            "Intensity": {"upper": 2, "lower": -2},
        }` - Specified by `coarse_params`. The coarse params for anchor determination and coarse matching.
7. `MSAligner.default_scaler_params`: `{
        "smoothing_method": "lowess",
        "smoothing_params": {"frac": 0.1},
    }` - Scaling method and params
8. `intensity_col`: `"Intensity"` - Not really used

If we were to change the column labels, remove intensity, and change the coarse params and matching window:
```python
from bmxp.eclipse import MSAligner
MSAligner.descriptors = {"RT": "linear", "MZ": "ppm"} # remove intensity
MSAligner.default_coarse_params['RT'] = {"upper": 1.0, "lower": -1.0} # wider coarse matching
MSAligner.default_cutoff = 10 # very generous +/-10 stdev matching for the descriptors
MSAligner.default_cutoffs = {'MZ': 6} # except MZ, force a +/- 6 stdev cutoff 
schema_labels = {
    "RT": "rt (min)",
}
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'], schema_labels=schema_labels) 
a.align()
```
##### Instance level defaults
If you would like to change a single instance vs globally, it's slightly different. And be sure to use the actual column name, not the schema key.

```python
from bmxp.eclipse import MSAligner

schema_labels = {
    "RT": "rt (min)", # the "RT" column is "rt (min)" in the dataset
}

# warning, if intensity is to be disabled it should be done globally since it will attempt
# to calculate it on initialization
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'], schema_labels=schema_labels)
a.descriptors = {"rt (min)": "linear", "MZ": "ppm"} # remove intensity

a.default_coarse_params["rt (min)"] = {"upper": 1.0, "lower": -1.0} # wider coarse matching
a.default_cutoff = 10 # very generous +/-10 stdev matching
a.default_cutoffs = {'MZ': 6} # except for mz
# or all at once -- this will not clear values not present, just update via dict.update()
# a.set_instance_defaults({
#   "coarse_params" : {"rt (min)": {"upper": 1.0, "lower": -1.0}},
#   "cutoff": 10,
#   "cutoffs": {'MZ': 6}
# })
a.align()
```

##### Granular subalignment-specific params
Subalignment-specific params, best modified with `set_params`. You may specify the source and target dataset. By default, those params will be applied to that subalignment as well as the reverse. The following param are valid: `coarse_params`, `scalers`, `cutoffs`, `weights`, `scaler_params`. 

Let's say DS1 ran a whole minute later than DS2 and DS3, but we don't want to increase the window for DS2 and DS3. We can do it for just that subalignment, like:
```python
from bmxp.eclipse import MSAligner
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'])
a.set_params({
    'coarse_params': {
        'DS1': {
            'DS2':{
                'RT': {'upper': 1.5, 'lower': 0.5}
            },
            'DS3':{
                'RT': {'upper': 1.5, 'lower': 0.5}
            }
        }
    }         
})
a.align()
```
The equivalent reverse (DS2->DS1, DS3->DS1) will be applied by default. It's a bit messy, but the only way for granular control over subalignments.

#### Datasets are loaded
Any number of datasets can be loaded. Loading over 10 datasets might require changing subalignment from all-to-all to one-to-all.
#### Anchors are determined

`MSAligner.anchors`: dict[str, tuple], a dictionary containing dataset names and a tuple containing the index names of the anchors.


```Python
from bmxp.eclipse import MSAligner
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'])
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
from bmxp.eclipse import MSAligner
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'])
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
from bmxp.eclipse import MSAligner
a = MSAligner('DS1.csv', 'DS2.csv', 'DS3.csv', names=['DS1', 'DS2','DS3'])
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

As an example, say we have coarse matches between two datasets (DS1->DS2) and it looks like so:
```
DS1 RT: 1.0, 2.0, 3.0, 4.0, 5.0
DS2 RT: 1.1, 2.2, 3.1, 4.2, 5.1
```
After fitting, the resulting scalers are:
```
DS1 Scalers:
RTx: 1.0, 2.0, 3.0, 4.0, 5.0
Rty: 0.1, 0.2, 0.1, 0.2, 0.1
```
Positive Y-values indicate that the source dataset has lower/smaller values than the target. 

In the reverse subalignment (DS2->DS1), the DS2 scalers determined are similar **but not identical**. 
```
DS2 Scalers:
RTx:  1.1,  2.2,  3.1,  4.2,  5.1
RTy: -0.1, -0.2, -0.1, -0.2, -0.1
``` 

#### Variance (Noise, Standard Deviation) Determination
`MSAligner.stds`: dict[str, dict[str, dict[str, float]]], a dictionary containing the residual standard deviations for the descriptors in a subalignment pair.

```Python
a.gen_stds()
print(a.stds["DS1"]["DS2"])
```

After scaling, there will still be remaining noise and outliers from a wide coarse matching window. Outliers will be removed by a MAD based algorithm, then standard deviations will be calculated. These standard deviations will be used for the actual match step, for both selection as well as scoring.

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
