# Formation

## Quick use
1. Install with `pip install bmxp` (requires Python 3)
2. Get an LCMS dataset with "RT" and sample abundances
3. Have fun!

```python
from bmxp import formation
pdf_bytes = formation.report('dataset.csv')
# returns BytesIO of PDF
```

This will open the dataset, generate a report, write the pdf as "dataset.csv-formation-report.pdf", and return the BytesIO object of the pdf.

## Overview
This module will eventually be a part of assembling (and de-assembling) Feature Metadata, Sample Metadata, Injection Metadata, and Sample Abundances into a final format, and generating a QC report. Currently, only the QC report exists.

### QC Report
This takes a formatted dataset and writes a pdf report. Params are:
* `dataset` - Pandas dataframe, .csv filepath or .xlsx filepath
* `write_pdf` - `Bool` - Optionally disables writing the pdf
* `out_filepath` - `String` - Specify the name of the output

The formatted Dataset is our standard combine feature and injection metadata pivot table. It might look like:

|     |     |  |  |  |  |  |  |
|-----|-----|-------|-------|-------|-------|-------|-------|
|             |     |          |                   | Date_Injected | 12-01-01 | 12-01-01 | 12-01-01  |
|             |     |          |                   | injection_type| PREFA    | PREFB    | Sample    |
|             |     |          |                   | Column        | 1        | 1        | 2         |
|             |     |          |                   | raw_file_name | one.raw  | two.raw  | three.raw |
| Compound_ID | RT  | MZ       | HMDB_ID           | Metabolite    | One      | Two      | Three     |
| cmp_001     | 1.5 | 100.0001 | internal_standard | Compound232   | 43.2     | 54.2     | 100.9     |
| cmp_002     | 1.6 | 200.0001 |                   | Compound10    | 38.1     | 74.2     | 345.7     |
| ...         | ... | ...      |                   | ...           | ...      | ...      | ...       |

There must be at least 2 feature metadata columns (such that the top left cell is empty), and one of them must be "Compound_ID". The Sample Metadata must have "raw_file_name" and an unnamed column ("reporting_name") with the sample names (e.g. "One", "Two" and "Three" in the example).

The report includes the following charts:
* PCA of Samples labaled by Name
* Loadings Plot labeled with Metabolite (knowns labeled)
* Loadings Plot labeled with Compound_ID (all labeled)
* PCAs colored by metadata
* Internal Standard abundance plots
* Sample Median abundance plot
* Box and whisker plots of feature abundance
* Feature CVs


