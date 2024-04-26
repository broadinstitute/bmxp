# Blueshift

Blueshift is a library that corrects for signal intensity drift over the course of a study or between multiple studies. Using internal standard-based drift correction and/or one of three methods of pool-based drift correction, blueshift corrects for changes in signal intensity due to changes in machine sensitivity, differences between LC columns, etc.

## Load Data and Sample Information

`DriftCorrection.rdata` : DataFrame; raw data and corresponding metadata for all samples

`DriftCorretion.sample_information` : DataFrame; sample information for all samples

When a `DriftCorrection` object is created, the constructor accepts a single dataset and corresponding sample information sheet. The dataset should contain columns for retention time, mass-to-charge ratio, sample values, and a unique identifier for each metabolite measured (`Compound_ID`). The sample information sheet should contain file names, sample IDs, injection order, and, optionally, any needed manual overrides to the drift correction algorithm (`breakpoints`; see "Pool-Based Drift Correction").

`DriftCorrection.pool_indices` : list; indices for pooled reference samples

Upon being loaded, the sample information sheet is used to identify the indices of pooled samples, which can be used to locate the corresponding rows in `sample_information` or columns in `rdata`.

## Internal Standard-Based Drift Correction

`DriftCorrection.cdata` : DataFrame; scaled data and corresponding metadata for all samples

The first drift correction method to be performed adjusts signal intensities for each sample based on the measured values for a given internal standard (IS). Scaling values are calculated by dividing the IS intensity in each sample by the average IS intensity across all samples. In the case where more than one internal standard is provided, scaling values are calculated for each and then the mean is taken. These values are then used to scale all non-IS metabolites.

Internal standard-based drift correction is optional. If selected, it can be performed in addition to or instead of reference pool-based drift correction.

## Reference Pool-Based Drift Correction

Reference pool-based drift correction methods can be used to account for signal drift in metabolites that don't follow the same trends as the internal standards, or in cases where no consistent internal standard is available. Unlike internal standard-based drift correction, where a single set of scalers are calculated and applied to the full dataset, reference pool-based drift correction calculates a separate set of scalers for each metabolite.

Two interpolation methods are available: nearest neighbor and linear. In addition, users will have the option to provide manual overrides to the algorithm in the form of breakpoints in the run order. Samples before a breakpoint will only ever be scaled to reference pools before the breakpoint, and samples after a breakpoint will only ever be scaled to reference pools after the breakpoint, regardless of what scalers the algorithm would typically select.

Only one reference pool-based drift correction method can be used on a single dataset. Reference pool-based drift correction is optional. If selected, it can be performed in addition to or instead of internal standard-based drift correction.

### Nearest Neighbor Interpolation

Scalers for each sample are calculated by locating the reference pool which is nearest (with regard to injection order) to the given sample and dividing by the average intensity across all reference pools. Scalers are calculated and applied on a per metabolite basis.

When indicated, breakpoints can be used to override the nearest neighbor selection. A variety of circumstances may require a manual override, including accounting for sudden changes in signal between samples (e.g. due to addition of new mobile phase or between two LC columns). In these cases, samples before the breakpoint will scale to the reference pools before the breakpoint, and samples after the breakpoint will scale to reference pools after the breakpoint, even if there is a reference pool with a closer injection order.

Similarly indicated in the sample information sheet, a user can choose to exclude specific reference pools that otherwise would have been used to correct nearby samples. This is useful to provide accurate drift correction in data where certain reference pools are known to be outliers or otherwise not reflective of the general trends in signal intensity.

### Linear Interpolation

Linear drift correction is the simplest "smooth" curve correction and provides a foundation for creating more complex (and smoother) smooth curve corrections if needed. For each reference pool, a line is identified that runs through that pool and the reference pool that comes next in the injection order. Intermediate values on each line are used to calculate scalars for samples that lie between the reference pools. Line drift correction is most effective when used on datasets where signal increases or decreases steadily and substantially between reference pools.

When indicated, breakpoints can be used to account for sudden drops or increases in signal within a dataset which is otherwise well-suited for linear interpolation. For a given breakpoint, the closest reference pools on either side of the breakpoint are located. Scalers for samples between these two reference pools are calculated using a method akin to nearest neighbor interpolation, while accounting for the breakpoint. Samples before the breakpoint are scaled directly to the nearest reference pool before the breakpoint, and samples after the breakpoint are scaled directly to the nearest reference pool after the breakpoint. Outside of this range, all other samples still use the typical linear interpolation.

## Coefficient of Variation Calculation

As a final step, coefficients of variation (CV) are calculated for each metabolite. CVs are calculated for drift correction reference pools, QC reference pools, and samples, for both raw and post-correction data. Comparison of raw and corrected CVs can be part of the process of assessing the accuracy and usefulness of corrected data.