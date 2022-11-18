# Gravity

## Quick use
1. Install with `pip install bmxp` (requires Python 3)
2. Get an LCMS dataset with "RT" and sample abundances
3. Have fun!

```python
from bmxp.gravity import cluster
df = pd.read_csv("your-file.csv", header=0)
results = cluster(df)
print(results)
```

This will run the data on default settings, `rt = 0.025`, `corr = 0.9`, and `method='Spearman'`.


## Overview
Gravity uses Retention Time and feature correlation across samples to identify redundant features. In the future, this might include peak profile.

### Acceptable Formats
Gravity must identify the retention time column and the sample columns. To do so, retention time must be `RT`, and the following columns are ignored:
`['RT', 'MZ', 'Intensity', "Non_Quant", "Compound_ID", "Adduct", "Annotation_ID", "Metabolite"]`
Everything else is assumed to be a sample column.

### Network Building Phase
Gravity works by building a network of features, based on the correlation and retention time being above a threshold. Gravity can handle large datasets and has been tested with 68,000 features by 130 samples, taking about 5 minutes to run.

### Deconvolution Phase
Once the network is built, the largest clique is identified and recorded as a cluster. These features are then removed from the graph, and the next largest clique is recorded. To break ties between cliques, the total correlation (stored as an edge attribute) is calculated and the higher one takes priority.
