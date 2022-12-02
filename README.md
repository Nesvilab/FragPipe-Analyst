# FragPipe-Analyst

FragPipe-Analyst is an easy-to-use, interactive web application developed to perform various computational analyses and to visualize quantitative mass spectrometry-based proteomic datasets processed using [FragPipe](https://fragpipe.nesvilab.org/) computational pipeline. It is compatible with the LFQ-MBR, spectral counts, TMT, and DIA quantification workflows in FragPipe. 

FragPipe-Analyst provides multiple useful analysis features, including differential expression (DE) analysis using Limma, enrichment analysis (GO/Pathways) using Enrichr, and multiple options for missing value imputation. It provides rich data visualization capabilities, including PCA and volcano plot, sample correlation plots, heatmaps, plots for missing value and sample coverage inspection, protein abundance plots for selected protein(s), imputation effect evaluation, and more. The user can also export results in table format for further analysis, such as DE results (FC and p-values), imputed data matrix (protein abundances after performing selected imputation method), and GO/pathway enrichment results.

## Public Servers
There are two server instances
- Production (stable) server is available at [http://fragpipe-analyst.org](http://fragpipe-analyst.org/).
- Development server is hosted at [http://fragpipe-analyst.nesvilab.org/](http://fragpipe-analyst.nesvilab.org/). It provides the most up to date version, including bug fixes, but may not be stable 

## Documentation, Questions, and Technical Support 
This repository is intended to provide documentation and a user forum for questions, suggestions, and bug reports. For bug report, please make sure you provide enough details for reproducing the errors. If you feel safer sharing data via email, please email it to [yihsiao@umich.edu](yihsiao@umich.edu). 

- [Latest Documentation/Tutorials](https://github.com/MonashProteomics/FragPipe-Analyst/tree/main/docs)
- [Questions/Suggestions](https://github.com/Nesvilab/FragPipe-Analyst/discussions)
- [Bug Reports](https://github.com/Nesvilab/FragPipe-Analyst/issues)

## Source Code
The source code of FragPipe-Analyst can be found [here](https://github.com/MonashProteomics/FragPipe-Analyst).
The code is available under GNU General Public License v3.0

## Credits
FragPipe-Analyst is based on the original LFQ-Analyst code developed by the Monash Proteomics and Metabolomics Facility & Monash Bioinformatics Platform, Monash University. FragPipe is being developed by the Nesvizhskii Lab (University of Michigan) in collaboration with the Monash University team (PI Ralf Schittenhelm).

