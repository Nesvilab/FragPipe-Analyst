# Missing value imputation options

Multiple imputation options were provided in FragPipe-Analyst

-   **Perseus-type:** This method is based on popular missing value
    imputation procedure implemented in *Perseus* software by MaxQuant
    team. The missing values are replaced by random numbers drawn from a
    normal distribution of *1.8* standard deviation down shift and with a
    width of *0.3* of each sample.
-   **bpca:** Bayesian missing value imputation
-   **knn:** Missing values replace by nearest neighbor averaging
    technique
-   **min:** Replaces the missing values by the smallest non-missing
    value in the data.
-   **zero:** Replaces the missing values by **0**.
