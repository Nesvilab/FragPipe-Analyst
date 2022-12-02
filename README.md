# FragPipe-Analyst Forum

FragPipe-Analyst is the downstream analysis tool desgined for [FragPipe](https://fragpipe.nesvilab.org/). It's currently under beta test.
This repo is intended to be used as a user forum for questions/suggestions/bug reports. For bug report, please make sure you provide enough details for reproducing the errors. If you feel more safe to share data via email, please email it to [yihsiao@umich.edu](yihsiao@umich.edu). 

- [Latest Documentation/Tutorials](https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs)
- [Questions/Suggestions](https://github.com/Nesvilab/FragPipe-Analyst/discussions)
- [Bug Reports](https://github.com/Nesvilab/FragPipe-Analyst/issues)

## Public Servers

There are two server instances
- Latest dev (unstable) server is hosted at [http://fragpipe-analyst.nesvilab.org/](http://fragpipe-analyst.nesvilab.org/).
- Production (stable) server is coming soon.

If you are interested in the source code of FragPipe-Analyst, you can find it [here](https://github.com/MonashProteomics/FragPipe-Analyst).

## Features

- Differential expression analysis
- Enrichment analysis (GO/Pathways)
- Imputation (optional)
- Data visualization
  1. PCA
  2. Sample correlation
  3. Heatmaps
  4. Missing value inspection
  5. Sample coverage
  6. Protein intensity plots for slected protein(s)
  7. Imputation effect evaluation

- Report and multiple levels of exportable tables for further analysis
  - Table options
    - DE results
    - Unimputed data matrix: Original protein intensities before imputation
    - Imputed data matrix: Protein intensities after performing selected imputation method
