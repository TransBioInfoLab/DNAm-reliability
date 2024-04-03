# **Critical evaluation of the reliability of DNA methylation probes on the Illumina MethylationEPIC v1.0 BeadChip microarrays** 
Wei Zhang, Juan I. Young, Lissette Gomez, Michael A. Schmidt, David Lukacsovich, Achintya Varma, X. Steven Chen, Brian Kunkle, Eden R. Martin, Lily Wang

#### Cite this article
Zhang, W., Young, J. I., Gomez, L., Schmidt, M. A., Lukacsovich, D., Varma, A., Chen, X. S., Kunkle, B., Martin, E. R., & Wang, L. (2024). Critical evaluation of the reliability of DNA methylation probes on the Illumina MethylationEPIC v1.0 BeadChip microarrays. *Epigenetics*, 19(1). https://doi.org/10.1080/15592294.2024.2333660

## Discription

This github repository includes scripts used for the analyses in the above manuscript.

DNA methylation (DNAm) plays a crucial role in a number of complex diseases. However, the reliability of DNAm levels measured using Illumina arrays varies across different probes. Previous research primarily assessed probe reliability by comparing duplicate samples between the 450k-450k or 450k-EPIC platforms, with limited investigations on Illumina EPIC v1.0 arrays. We conducted a comprehensive assessment of the EPIC v1.0 array probe reliability using 69 blood DNA samples, each measured twice, generated by the Alzheimer's Disease Neuroimaging Initiative study. We observed higher reliability in probes with average methylation beta values of 0.2 to 0.8, and lower reliability in type I probes or those within the promoter and CpG island regions. Importantly, we found that probe reliability has significant implications in the analyses of Epigenome-wide Association Studies (EWAS). Higher reliability is associated with more consistent effect sizes in different studies, the identification of differentially methylated regions (DMRs) and methylation quantitative trait locus (mQTLs), and significant correlations with downstream gene expression. Moreover, blood DNAm measurements obtained from probes with higher reliability are more likely to show concordance with brain DNAm measurements. Our findings, which provide crucial reliable information for probes on the EPIC v1.0 array, will serve as a valuable resource for future DNAm studies. 

### 1. Preprocessing of DNA methylation data

We pre-processed the samples using the SeSAMe 2 pipeline described in Welsh et al. (2023), which was found to perform the best and produced the largest percentage of reliable CpG probes in a recent comparison of various pre-processing and normalization pipelines [1]. 

| File and folder                                              | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [code/preprocessing](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/preprocessing) | Code for preprocessing of ADNI DNAm duplicated blood samples |
### 2. **Estimation of probe reliability and surrogate variables** 

The reliability is calculated by intraclass correlation (ICC). The ICC is defined as $\frac{\sigma_{b}^2}{\sigma_{b}^2 + \sigma_{w}^2}$, where $\sigma_{b}^2$ is the between-subject variance and $\sigma_{w}^2$ is the within-subject variance. As recommended by Koo and Li (2016) [2], ICC values were computed using a two-way random effect, absolute agreement, and single-rating model, as implemented in the irr R package. The cell-type proportions are computed using EpiDISH R package [3]. The coefficients for DNA methylation-based surrogate variables for BMI, smoking, alcohol use, total cholesterol, HDL cholesterol, LDL cholesterol, and total-to-HDL ratio were obtained from Additional file 1 in McCartney et al. (2018) [4].

| File and folder                                              | Description                                  |
| ------------------------------------------------------------ | -------------------------------------------- |
| [code/icc_calculation](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/icc_calculation) | Code for probe ICC calculation               |
| [code/surrogate_variables](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/surrogate_variables) | Code for surrogate variables ICC calculation |
| [results/ADNI_ICC_annotation.xlsx](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/results/ADNI_ICC_annotation.xlsx) | Results for probe reliability and annotation |

### 3. Distribution of the ICCs in EPIC-EPIC comparison 

To study the ICCs for probes in our EPIC-EPIC comparison, we compared the distributions of the ICCs in our EPIC-EPIC comparison with studies from Sugden et al. (2020) [5], Logue et al. (2017) [6], and Bose et al. (2014) [7] and several characteristics of probes (i.e. probe type). Furthermore, we evaluated the ICC distribution with our previous DNAm-AD association studies [8-9], DNAm-to-mRNA correlations, and blood-brain DNAm correlations.

| File and folder                                              | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [code/ADNI_reliability/ADNI_reliable_property.Rmd](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/ADNI_reliability/ADNI_reliable_property.Rmd) | Code for ICC comparisons with other studies and probe characteristics |
| [code/ADNI_reliability/ADNI_reliable_ICC_CpGs_vs_AD_study.Rmd](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/ADNI_reliability/ADNI_reliable_ICC_CpGs_vs_AD_study.Rmd) | Code for ICC comparisons with DNAm-AD association studies, and blood-brain DNAm correlations. |
| [code/DNAm_vs_RNA/DNAm_vs_RNA.Rmd](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/DNAm_vs_RNA/DNAm_vs_RNA.Rmd) | Code for ICC comparisons with DNAm-to-mRNA correlations      |
### 4. Comparison of reliability with different characteristics

To compare the reliability of probes with different characteristics (e.g., type I probes vs. type II probes), we performed mixed effects model analyses.

| File and folder                                                       | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [code/compare_test](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/compare_test) | Mixed-effects model analyses for reliability of probes and different characteristics |

## For reproducible research

To perform the analysis, begin by installing the packages found in `session_info.R` ([Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/session_info.R)). Then, load the auxiliary functions from `Utility.R` ([Link to the script](https://github.com/TransBioInfoLab/DNAm-reliability/blob/main/code/Utility.R)). Follow the sequence provided in the Description to conduct the analysis.

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

1. Welsh H, Batalha C, Li W, Mpye KL, Souza-Pinto NC, Naslavsky MS, Parra EJ: **A systematic evaluation of normalization methods and probe replicability using infinium EPIC methylation data.** *Clin Epigenetics* 2023, **15**:41.
2. Koo TK, Li MY: **A Guideline of Selecting and Reporting Intraclass Correlation Coefficients for Reliability Research.** *J Chiropr Med* 2016, **15**:155-163.
3. Teschendorff AE, Breeze CE, Zheng SC, Beck S: **A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies.** *BMC Bioinformatics* 2017, **18**:105.
4. McCartney DL, Hillary RF, Stevenson AJ, Ritchie SJ, Walker RM, Zhang Q, Morris SW, Bermingham ML, Campbell A, Murray AD, et al: **Epigenetic prediction of complex traits and death.** *Genome Biol* 2018, **19**:136.
5. Sugden K, Hannon EJ, Arseneault L, Belsky DW, Corcoran DL, Fisher HL, Houts RM, Kandaswamy R, Moffitt TE, Poulton R, et al: **Patterns of Reliability: Assessing the Reproducibility and Integrity of DNA Methylation Measurement.** *Patterns (N Y)* 2020,**1**.
6. Logue MW, Smith AK, Wolf EJ, Maniates H, Stone A, Schichman SA, McGlinchey RE, Milberg W, Miller MW: **The correlation of methylation levels measured using Illumina 450K and EPIC BeadChips in blood samples.** *Epigenomics* 2017, **9**:1363-1371.
7. Bose M, Wu C, Pankow JS, Demerath EW, Bressler J, Fornage M, Grove ML, Mosley TH, Hicks C, North K, et al: **Evaluation of microarray-based DNA methylation measurement using technical replicates: the Atherosclerosis Risk In Communities (ARIC) Study.** *BMC Bioinformatics* 2014, **15**:312.
8. T CS, Young JI, Zhang L, Gomez L, Schmidt MA, Varma A, Chen XS, Martin ER, Wang L: **Cross-tissue analysis of blood and brain epigenome-wide association studies in Alzheimer's disease.** *Nat Commun* 2022, **13**:4852.
9. T CS, Zhang W, Young JI, Gomez L, Schmidt MA, Varma A, Chen XS, Martin ER, Wang L: **Distinct sex-specific DNA methylation differences in Alzheimer's disease.** *Alzheimers Res Ther* 2022, **14**:133.
