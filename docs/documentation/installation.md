## Installation

Two options are available right now for FragPipe-Analyst: local installation or run through Docker.

### Local Installation

#### Prerequisite
- R >= 4.4
- PDFlatex

Once all the prerequisites are installed, follow steps below to build and run the server locally.

``` sh
# Clone the repository
git clone https://github.com/MonashProteomics/FragPipe-Analys.git

# Move to the folder
cd FragPipe-Analyst

# Inside R console or R studio
> install.packages("renv")
> renv::init(bioconductor = T)

# Execute
> shiny::runApp()
```

### Installation through Docker:

``` sh
# Clone the repository
git clone https://github.com/MonashProteomics/FragPipe-Analyst.git

# Move to the folder
cd FragPipe-Analyst

# Build FragPipe-Analyst (Any name after -t)
docker buildx build -f Dockerfile.local -t fragpipe-analyst  --output=type=docker --platform=linux/amd64 .

# Run FragPipe-Analyst
docker run -it --platform=linux/amd64 -d -p 3838:3838 fragpipe-analyst

# Open local interface
http://localhost:3838/fragpipe-analyst/
```
