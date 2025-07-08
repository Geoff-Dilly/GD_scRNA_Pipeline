# Install_Packages.R
# Install all required packages for the scRNA-seq pipeline

# CRAN packages
cran_packages <- c(
  "BiocManager",
  "here",
  "dplyr",
  "patchwork",
  "ggplot2",
  "stringr",
  "data.table",
  "knitr",
  "devtools"
)

# Bioconductor packages
bioconductor_packages <- c(
  "Seurat",
  "SoupX",
  "DropletUtils",
  "DESeq2",
  "glmGamPoi"
)

# Install CRAN packages
for (package in cran_packages){
  install.packages(package)
}

# Install Bioconductor packages
for (package in bioconductor_packages){
  BiocManager::install(package, ask = FALSE, update = FALSE)
}

# Install DoubletFinder from GitHub
if (!"DoubletFinder" %in% rownames(installed.packages())) {
  devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
}

# Install Presto from GitHub
if (!"presto" %in% rownames(installed.packages())) {
  devtools::install_github("immunogenomics/presto")
}