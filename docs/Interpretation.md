# Interpret the result

## Principal Component Analysis (PCA)

![PCA Plot](Images/PCA_plot.png)

PCA plot is a good way to inspect your result regarding batch effect and reproducibility of samples from same condition or eplicates. One limitation of PCA is that only full observed proteins (may by imputation) are included in the PCA.  

## Missing Value Heatmap

![Missing Value Heatmap](Images/missing_heatmap.png)

Missing value heatmap shows the missing value pattern of samples. Note that only protein with missing values are shown in the heatmap.

## Correlation Heatmap

![Correlation Plot](Images/correlation_plot.png)

Correlation heatmap shows the correlation matrix between samples. Typically, we expect samples in the same condition should be clustered closer together.

## Coefficient of Variation Plot

![CV plot](Images/CV_plot.png)

Coefficient of variation plot shows distribution of protein level coefficient of variation for each condition. Each plot also contains a vertical line representing median CVs percentage within that condition. Low coefficient of variation means the spread of data values is low relative to the mean.

## Volcano Plot

![Volcano Plot](Images/volcano_plot.png)

Volcano plot is the most common used visualization for differentiatl expression analysis.

## Enrichment Result

![Enrichment Plot](TMT-tutorial/GO_term.png)

Enrichment analysis gives you a rough idea of changes between two different conditions. Our enrichment analysis is based on [Enrichr](https://maayanlab.cloud/Enrichr/) with background correction through hypergeometric test.
