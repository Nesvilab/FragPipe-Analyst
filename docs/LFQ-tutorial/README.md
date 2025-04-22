## LFQ data analysis using FragPipe-Analyst

After processing data in FragPipe as mentioned in the [last part of the tutorial](LFQ.md), here we show how to further analyze the data using FragPipe-Analyst for downstream analysis and visualization (PCA, differential expression, etc.). 

Just as a recap, a subset of a published gliomas dataset described by: [J. Bader et al.](https://pubmed.ncbi.nlm.nih.gov/36584682/) was used when running FragPipe.

J. Bader et al. “Proteomics separates adult-type diffuse high-grade gliomas in metabolic subgroups independent of 1p/19q codeletion and across IDH mutational status”, Cell Rep Med 2023 4(1):100877. doi: 10.1016/j.xcrm.2022.100877. 

In the original study, authors studied high-grade adult-type diffuse gliomas are malignant neuroepithelial tumors with poor survival rates in combined chemoradiotherapy. They used MS1-based label-free quantification (LFQ) mass spectrometry to characterize 42 formalin-fixed, paraffin-embedded (FFPE) samples from IDH-wild-type (IDHwt) gliomas, IDH-mutant (IDHmut) gliomas, and non-neoplastic controls. Here, we just used 6 samples, 3 IDHmut and 3 IDHwt earlier when preparing the data and now we are about to further perform downstream analysis. 

For the following steps, you can use the generated `combined_protein.tsv` and `experiment_annotation.tsv` or download them from [here](https://www.dropbox.com/sh/38rgczrcmdrcnm9/AADCiLYKhfouSFk-LD4gR-D0a?dl=1).

- Open FragPipe-Analyst (there is a link in the Run tab) or visit [http://fragpipe-analyst.nesvilab.org/](http://fragpipe-analyst.nesvilab.org/)
- Click `Analysis` on the left-hand side menu
- Choose `LFQ` in the dropdown menu for the data type  
- Upload `combined_protein.tsv` and `experiment_annotation.tsv`. That's a subset of a published dataset described in the following publication: https://pubmed.ncbi.nlm.nih.gov/36584682/from J. Bader et al. “Proteomics separates adult-type diffuse high-grade gliomas in metabolic subgroups independent of 1p/19q codeletion and across IDH mutational status”, Cell Rep Med 2023 4(1):100877. doi: 10.1016/j.xcrm.2022.100877.
- In the `Intensity Type`, change it to MaxLFQ Intensity.
- Click the `Run` button

When we work on such project, one would typically start with explorative analyses such as Principal Component Analysis (PCA) to see if the protein data exhibit tumor/normal difference. Following that, one would look for known and potential diagnostic markers for various tumor subtypes with differential expression (DE) analysis, comparing the expression of each protein in tumor samples of one type to that of other types (or with that of normal samples). If for a protein A there is a significant difference between the expression of the two groups, A is seen as a potential marker. Since the comparison is done for many genes, multiple test adjustment is implemented to control the overall false discovery rate for differential expression. 

And here are some findings reported in the paper by the authors related to IDHmut vs IDHwt comparisons:

The IDHwt gliomas in our study were most distinct from CNS ctrl and also segregated from IDHmut gliomas (Figures 1B and 1C). The IDHwt proteome was enriched with proteins linked to inflammation, MCM complex DNA polymerases, an integrin-, collagen-, and laminin-rich ‘‘basement membrane-like’’ extracellular matrix (ECM) profile, low in hyaluronic acid, which is associated with increased malignancy in gliomas (Figures 1C–1E). IDHwt/IDHmut differences aligned well with the CPTAC data and were largely unaffected by 1p/19q codeletion status in our data (Figures 1E, S1C, and S1D). In line with amore ‘‘aggressive’’ phenotype of IDHwt gliomas, many outlier proteins with high abundance in IDHwt are cancer drivers, several of them linked to invasion. Outlier proteins associated with IDHmut included tumor suppressors downregulated in IDHwt (Figures 1E andS1D; Table 2). Notably, these tumor suppressors include the histone proteins H1F0 and H2AFY2, which both maintain an epigenetic profile of differentiation and inhibit a return to a proliferative stem cell state through distinct mechanisms. Known proteome alterations driving progression of the IDHmut were apparent, such as the strong epigenetic downregulation of RBP1, as well as novel ones, such as AKR1C3 overexpression selectively in IDHmut (Figures 1E and S1E; Table 2). 

### Questions

Explore the data in FragPipe-Analyst to answer the following questions:

- Does the proteome data exhibit differences between the IDHwt and IDHmut samples? Inspect the PCA plot and the heat map.
Inspect ‘Protein Numbers’ and ‘Missing Values- heatmap’ plots to see if there are any issues with the data (e.g., too few proteins identified in one of the plexes).

- Check some proteins they highlighted in the paper. For example, RBP1 and AKR1C3 mentioned above. See them on the volcano plot. Check Limma computed adjusted p-values. Check their expression levels across the individual samples using ‘Protein Plot’.

- With differential analysis conducted, we get lists of up- and down-regulated genes. It is often difficult to make sense of individual genes, especially when there are many. Enrichment analysis enables us to aggregate the evidence to biological pathway level to gain a higher-level insight of tumor features.

- Conduct an Enrichment Analysis. Focus on the most significantly enriched pathways. Compare with what is discussed in the paper (see above, especially Figure 1E; see Appendix).

### Notes
- For LFQ, a Perseus-like imputation is used by default prior to Limma differential expression analysis.
- FragPipe-Analyst allows using the following databases for enrichment analysis: KEGG, Reactome, HALLMARK. For this tutorial, we recommend HALLMARK.
- You can use Reactome, or DAVID, or GSEA, or any tool you like directly. To use an external tool, download the results of the differential expression from FragPipe-Analyst by clicking on Download Report. Or you can also use the enrichment analysis tool that is available as part of FragPipe-Analyst. Choose ‘Pathways analysis’, select ‘KEGG’, and run. Then download the full table of enriched pathways by clicking on Download Table, since the plot you will see shows just a few selected pathways.

