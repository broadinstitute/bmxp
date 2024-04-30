# BMXP - The Metabolomics Platform at the Broad Institute

`pip install bmxp`

This is a collection of tools we use daily for processing our data. Each tool performs a step in our processing workflow. Currently, we have published

* [Eclipse](https://github.com/broadinstitute/bmxp/tree/main/bmxp/eclipse) - Align two or more same-method nontargeted LCMS datsets.
* [Gravity](https://github.com/broadinstitute/bmxp/tree/main/bmxp/gravity) - Cluster redundant LCMS features based on RT and Correlation (And someday, XIC shape)
* [Blueshift](https://github.com/broadinstitute/bmxp/tree/main/bmxp/blueshift) - Drift Correction via pooled technical replicates and internal standards

In the future we will add:

* Rawfilereader - C-optimized LCMS raw file reader for exploring raw files and formatting MS2

