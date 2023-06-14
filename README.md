# Critical evaluation of the reliability of DNA methylation probes on the Illumina MethylationEPIC BeadChip microarrays for dementia research

Wei Zhang, Juan I. Young, Lissette Gomez, Michael A. Schmidt, David Lukacsovich, Achintya Varma, X. Steven Chen, Brian Kunkle, Eden R. Martin, Lily Wang

## Discription

This github repository includes scripts used for the analyses in the above manuscript.

DNA methylation (DNAm) plays a crucial role in many complex diseases, including dementia. Array-based technologies offer a cost-effective and comprehensive approach to measure DNA methylation profiles on a genome-wide scale. However, the accuracy of DNAm measurements obtained using Illumina arrays can vary across different probes. Existing research has primarily focused on assessing the reliability of DNAm in younger subjects, leaving a knowledge gap regarding the reliability of CpG probes in older subjects for dementia research. Additionally, most previous studies have compared duplicate samples between the 450k-450k or 450k-EPIC platforms, with limited investigations on EPIC-EPIC comparisons.

We conducted a comprehensive assessment of probe reliability on the Illumina EPIC arrays using 138 duplicated blood DNAm samples from the Alzheimer's Disease Neuroimaging Initiative (ADNI) study. We carefully selected DNA methylation samples from independent older subjects aged 65 years and above, and included duplicate samples that were distributed across different methylation plates within the ADNI dataset. To assess reliability of each probe, we computed intraclass correlations (ICCs) for each probe based on methylation beta values, which were measured on duplicates of blood samples collected from the same subject and at the same visit. Both the magnitude and patterns of reliability observed in the EPIC-EPIC comparison were compared with the findings of previous studies. Furthermore, we also investigated the impact of probe reliability on the analyses of epigenome-wide association studies (EWAS).



### 1. Preprocessing of DNA methylation data and estimation of probe reliability

We pre-processed the samples using the SeSAMe 2 pipeline described in Welsh et al. (2023), which was found to perform the best and produced the largest percentage of reliable CpG probes in a recent comparison of various pre-processing and normalization pipelines [1]. To estimate the reliability of the CpG probes, we computed intraclass correlations (ICCs) for each probe based on methylation beta values, which were measured on duplicates of blood samples collected from the same subject and at the same visit. 

| File                                 | Link                                                         |
| ------------------------------------ | ------------------------------------------------------------ |
| preprocessing/ADNI_preprocessing.Rmd | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/preprocessing/ADNI_preprocessing.Rmd) |



### 2. Distribution of the ICCs in EPIC-EPIC comparison 

To study the ICCs for probes in our EPIC-EPIC comparison, we compared the distributions of the ICCs in our EPIC-EPIC comparison with results from Sugden et al. (2020) [2], Logue et al. (2017) [3], and Bose et al. (2014) [4]. 

| File                                        | Link                                                         |
| ------------------------------------------- | ------------------------------------------------------------ |
| ADNI_reliability/ADNI_reliable_property.Rmd | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/ADNI_reliability/ADNI_reliable_property.Rmd) |



### 3. Half-width of the 95% confidence Limits of Agreement (HoLA) and modified ICC

To explicitly account for the amount of disagreement, we proposed to assess reliability using the modified ICC, which is defined as ICC – Half-width of 95% confidence Limits of Agreement (or HoLA). HoLA is calculated as $1.96\sigma_d$, where $\sigma_d$ is the standard deviation of the differences between the two duplicate measures. Moreover, we evaluated the impact of modified ICC on mQTL analysis, dementia studies, DNAm-to-gene expression correlations, and surrogate variables 

| File                                                       | Link                                                         |
| ---------------------------------------------------------- | ------------------------------------------------------------ |
| ADNI_reliability/ADNI_mICC_CpGs_property.Rmd               | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/ADNI_reliability/ADNI_mICC_CpGs_property.Rmd) |
| ADNI_reliability/ADNI_reliable_mICC_CpGs_vs_ADNI_study.Rmd | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/ADNI_reliability/ADNI_reliable_mICC_CpGs_vs_ADNI_study.Rmd) |
| DNAm_vs_RNA/DNAm_vs_RNA.Rmd                                | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/DNAm_vs_RNA/DNAm_vs_RNA.Rmd) |
| cell_type/cell_type_ICC.Rmd                                | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/cell_type/cell_type_ICC.Rmd) |



### 4. Comparison of reliability with different characteristics

To compare the reliability of probes with different characteristics (e.g., type I probes vs. type II probes), we performed mixed effects model analyses.

| File                          | Link                                                         |
| ----------------------------- | ------------------------------------------------------------ |
| compare_test/compare_test.Rmd | [Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/compare_test/compare_test.Rmd) |



### 5. Results for probe reliability, HoLA and modified ICC

| File                                               | Link                                                         |
| -------------------------------------------------- | ------------------------------------------------------------ |
| results/ADNI_ICC_annotation_with_HoLA_grouped.xlsx | [Link to the result](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/results/ADNI_ICC_annotation_with_HoLA_grouped.xlsx) |



## For reproducible research

The following R packages are required:

```r
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
```

The script for installing the packages can be found at: code/session_info.R ([Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/session_info.R))

The platform information is:

```r
version  R version 4.2.1 (2022-06-23)
os       macOS Ventura 13.3.1
system   x86_64, darwin17.0
ui       RStudio
language (EN)
collate  en_US.UTF-8
ctype    en_US.UTF-8
tz       America/New_York
date     2023-06-14
rstudio  2023.03.1+446 Cherry Blossom (desktop)
pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
```



## Acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf



## Reference

1. Welsh H, Batalha C, Li W, Mpye KL, Souza-Pinto NC, Naslavsky MS, Parra EJ: **A systematic evaluation of normalization methods and probe replicability using infinium EPIC methylation data.** *Clin Epigenetics* 2023, **15:** 41.
2. Sugden K, Hannon EJ, Arseneault L, Belsky DW, Corcoran DL, Fisher HL, Houts RM, Kandaswamy R, Moffitt TE, Poulton R, et al: **Patterns of Reliability: Assessing the Reproducibility and Integrity of DNA Methylation Measurement.** *Patterns (N Y)* 2020, **1**.
3. Logue MW, Smith AK, Wolf EJ, Maniates H, Stone A, Schichman SA, McGlinchey RE, Milberg W, Miller MW: **The correlation of methylation levels measured using Illumina 450K and EPIC BeadChips in blood samples.** *Epigenomics* 2017, **9:** 1363-1371.
4. Bose M, Wu C, Pankow JS, Demerath EW, Bressler J, Fornage M, Grove ML, Mosley TH, Hicks C, North K, et al: **Evaluation of microarray-based DNA methylation measurement using technical replicates: the Atherosclerosis Risk In Communities (ARIC) Study.** *BMC Bioinformatics* 2014, **15:** 312.
