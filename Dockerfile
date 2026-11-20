# syntax=docker/dockerfile:1
FROM rocker/shiny-verse:4.4
RUN apt-get update && apt-get install -yq \
    libbz2-dev \
    libhdf5-dev \
    libnetcdf-dev \
    build-essential \
    libgd-dev \
    libudunits2-dev \
    libproj-dev \
    libgdal-dev  \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra

RUN Rscript -e 'install.packages(c("devtools", "BiocManager", "tidyverse", "ggrepel", "httr", "rjson", "mvtnorm", "tmvtnorm", \
"plotly", "DT", "testthat", "RColorBrewer", "shiny","shinyalert","shinydashboard", \
"shinyjs", "svglite", "rhandsontable", "shinyBS", "shinyWidgets", "ggVennDiagram", "conflicted", "png", "vegan", "assertthat", "jsonlite", \
"shinycssloaders","data.table", "factoextra", "UpSetR","fastcluster","fdrtool"), dependencies=TRUE)'
### #FROM bioconductor/bioconductor_docker:RELEASE_3_19
# Install from the source archive rather than install_github(), which calls the
# GitHub API and fails with HTTP 403 when CI runners hit its rate limit
RUN Rscript -e 'remotes::install_url("https://github.com/Appsilon/shiny.info/archive/f909a3f3142418bf53a02ca06ecd2d3d28438b7c.tar.gz")'
RUN Rscript -e 'BiocManager::install(pkgs=c("ensembldb", "EnsDb.Hsapiens.v86", "SummarizedExperiment", "limma", "ComplexHeatmap", "MSnbase", "missForest", "impute", "pcaMethods", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "ReactomePA", "msigdbr", "imputeLCMD"), ask=F, dependencies=TRUE)'

### RUN Rscript -e 'install.packages("renv")'
COPY ./ /srv/shiny-server/fragpipe-analyst
COPY shiny-server.conf.prod /etc/shiny-server/shiny-server.conf
### WORKDIR /srv/shiny-server/fragpipe-analyst
### RUN Rscript -e 'renv::init()' # This is already run and pushed to github.
### RUN Rscript -e 'renv::restore()'
### RUN Rscript -e 'renv::isolate()'
RUN chmod -R +r /srv/shiny-server/fragpipe-analyst
