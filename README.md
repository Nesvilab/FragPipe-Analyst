# FragPipe-Analyst

FragPipe-Analyst is an easy-to-use, interactive web application developed to perform differential expression analysis with “one click” and to visualize quantitative proteomic datasets analyzed using [FragPipe] (https://fragpipe.nesvilab.org/) computational platform. It is compatible with the LFQ-MBR, TMT, and DIA quantification workflows in FragPipe. 

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


## Public Servers

There are two server instances
- Production (stable) server is available at [http://fragpipe-analyst.org] (http://fragpipe-analyst.org/).
- Development server is hosted at [http://fragpipe-analyst.nesvilab.org/](http://fragpipe-analyst.nesvilab.org/). It provides the most up to date version, including bug fixes, but may not be stable 

## Documentation, Questions, and Technical Support 
This repository is intended to provide documentation and a user forum for questions, suggestions, and bug reports. For bug report, please make sure you provide enough details for reproducing the errors. If you feel safer sharing data via email, please email it to [yihsiao@umich.edu](yihsiao@umich.edu). 

- [Latest Documentation/Tutorials](https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs)
- [Questions/Suggestions](https://github.com/Nesvilab/FragPipe-Analyst/discussions)
- [Bug Reports](https://github.com/Nesvilab/FragPipe-Analyst/issues)

## Source Code
The source code of FragPipe-Analyst can be found [here](https://github.com/MonashProteomics/FragPipe-Analyst).

## Credits

FragPipe-Analyst is based on the original LFQ-Analyst code developed by the Monash Proteomics and Metabolomics Facility & Monash Bioinformatics Platform, Monash University. FragPipe is being developed by the Nesvizhskii Lab (University of Michigan) in collaboration with the Monash University team (PI Ralf Schittenhelm).

