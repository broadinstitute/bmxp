# BMXP - The Metabolomics Platform at the Broad Institute
`pip install bmxp`

Please cite:
https://www.biorxiv.org/content/10.1101/2023.06.09.544417v1.full

This is a collection of tools for processing our data, which powers our cloud processing workflow. Each tool is meant to be a standalone module that performs a step in our processing pipeline. They are written in Python and C, and designed to be perfomant and cloud-compatible.

* [Eclipse](https://github.com/broadinstitute/bmxp/blob/main/bmxp/eclipse/readme.md) - Align two or more same-method nontargeted LCMS datsets.
* [Gravity](https://github.com/broadinstitute/bmxp/blob/main/bmxp/gravity/readme.md) - Cluster redundant LCMS features based on RT and Correlation (And someday, XIC shape)
* [Blueshift](https://github.com/broadinstitute/bmxp/blob/main/bmxp/blueshift/readme.md) - Drift Correction via pooled technical replicates and internal standards
* [Formation](https://github.com/broadinstitute/bmxp/blob/main/bmxp/formation/readme.md) - Formatting and Final QC
* [Chroma](https://github.com/broadinstitute/bmxp/blob/main/bmxp/chroma/readme.md) - Read .raw and .mzml files

We expect users to be familiar with Python and already have an understanding of LCMS Metabolomics data processing and the specific steps they wish to accomplish.

While the tools are and always will be standalone, we are working on linking them closer together with a shared schema, and eventually may have a pipeline ability to run all steps, given a set of parameters.

We are open to feedback and suggestions, with a focus on performance and application in pipelines.

# Shared Schema
All BMXP modules use a shared schema and file formats with our prefered columns headers. These files are (along with their labels):
* Feature Metadata `bmxp.FMDATA` - Describes the feature. Index default is `Compound_ID`
* Injection Metadata `bmxp.IMDATA` - Describes the Injection. Index default is `injection_id`
* Sample Metadata `bmxp.SMDATA` - Describes the biospecimen from which the Injection is derived. Index default is `broad_id` 
* Feature Abundances - Pivot table of Feature x Injection (`Compound_ID` x `injection_id`) containing the abundances.

Some modules (Blueshift, Eclipse) require merging Feature Metadata + Feature Abundances.
 
These can be changed globally so that all packages will use the same terminology.
To update the schema, modify the dictionary objects in the module directly prior to running code. For example:
```python
import bmxp
from bxmp.eclipse import MSAligner
from bxmp.blueshift import DriftCorrection
from bmxp.gravity import cluster
bmxp.FMDATA['Compound_ID'] = 'Feature_ID'
bmxp.IMDATA['injection_id'] = 'Filename'

# continue with work...
```
With those changes above, Eclipse, Blushift and Gravity will use "Feature_ID" and "Filename" as column headers instead of "Compound_ID" and "injection_id".

## Feature Metadata - bmxp.FMDATA
Feature Metadata describes the LCMS feature. This is a mixture of fundamental nontargeted feature information, annotation info, and anything else.

### Feature Specific
* `Compound_ID` - Index, Project-unique feature ID (a bit of a misnomer)
* `RT` - Unitless retention time, may or may not be scaled
* `MZ` - Unsigned mass-to-charge ratio
* `Intensity` - Average feature intensity
* `Method` - Human Readable name of LCMS method used
* `__extraction_method` - Name of extraction method/software used. Used to denote mixed Targeted/Nontargeted

### Annotation
* `Annotation_ID` - Method-unique annotation label
* `Adduct` - Adduct form of the annotation
* `__annotation_id` - Globally unique annotation identifier
* `Metabolite` - Preferred display/reporting name of metabolite
* `Non_Quant` - Boolean denoting that a feature is not quanitifiable

### Generated by Gravity
* `Cluster_Num` - Cluster number assigned during Gravity clustering
* `Cluster_Size` -  Number of members in the cluster

### Generated by Blueshift
* `Batches Skipped` - Batches that were skipped due to lack of PREFs

## Injection Metadata - bmxp.IMDATA
* `injection_id` - Index, Injection name, usually filename without the extension
* `broad_id` - Assigned biospeciemn label
* `program_id` - Biospecimen label as received (inherited from Sample Metadata)
* `injection_type` - Type of injection ("sample", "prefa", "prefb", "blank", "other-", "not_used-")
* `comments` - Comments about the injection
* `column_number` - Column number, in multi-column studies
* `injection_order` - Injection number, not skipping blanks or non-samples
* `batches` - Denotes batches ('batch start' or 'batch end')

## generated by blueshift
* `QCRole` - Role in drift correction ("QC-drift_correction", "QC-pooled_ref", "QC-not_used", "sample")

## Sample Metadata - bmxp.SMDATA
* `broad_id` - Assigned biospecimen label
* Arbitrary Metadata Columns - Any column label except labels in Injection Metadata