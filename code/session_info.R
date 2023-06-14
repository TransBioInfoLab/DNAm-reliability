# -----------------------------------------------------------------------------------------------------------
# For reproducible research, please install the following R packages 
# and make sure the R and BiocManager versions are correct
# Session Info ----------------------------------------------------------------------------------------------
# setting  value 
# version  R version 4.2.1 (2022-06-23)
# os       macOS Ventura 13.3.1
# system   x86_64, darwin17.0
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2023-06-07
# rstudio  2023.03.1+446 Cherry Blossom (desktop)
# pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# -----------------------------------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.15")

list.of.packages <- c(
  "coMethDMR",
  "data.table",
  "devtools",
  "doParallel",
  "dorothea",
  "ExperimentHub",                                
  "GEOquery",                                     
  "ggpubr",    
  "gmodels",
  "GWASTools",  
  "hgu219.db",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "irr",
  "lme4",
  "lmerTest",
  "matrixStats",
  "MethReg",
  "plyr",                                         
  "readxl", 
  "RnBeads",
  "RnBeads.hg19",
  "sesame",
  "sesameData",
  "stats",                                        
  "SummarizedExperiment",                         
  "tidyverse",                                        
  "writexl"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  for(new in new.packages){
    if(new %in% available.packages()[,1]){
      install.packages(new)
    } else BiocManager::install(new)
  }
} 

