# Utility functions
## ----------------------------------------------------------------------------------------------------------------------------------------
## Create ICC value table
## ----------------------------------------------------------------------------------------------------------------------------------------
create_rating <- function(samples.id, beta, model = "twoway", type = "agreement", unit = "average", random = F, scale = F, id = "cpg"){
  
  barcodes <- samples.id[,c("barcodes1", "barcodes2")]
  if(random){
    sam <- sample.int(size = nrow(barcodes), n = 2, replace = T)
    barcodes <- sapply(1:nrow(barcodes),  function(i) c(barcodes[i,sam[i]],barcodes[i,-sam[i]])) %>% t() 
    rownames(barcodes) <- samples.id$RID_Phase_Edata
  }
  
  
  doParallel::registerDoParallel(20)
  rating_matrix <- plyr::adply(
    beta,
    .margins = 1,
    .fun = function(probe){
      rating <- plyr::aaply(
        barcodes,
        .margins = 1, 
        .fun = function(bar){
          probe[bar %>% unlist()]
        },.expand = F
      )
      if(scale) rating <- scale(rating)
      ic <- icc(rating, model = model, type = type, unit = unit)
      data.frame(value = ic$value, r0 = ic$r0, Fvalue = ic$Fvalue,
                 df1 = ic$df1, df2 = ic$df2, p.value = ic$p.value, 
                 lbound = ic$lbound, ubound = ic$ubound)
    },.parallel = T, .id = id
  )
  
  return(rating_matrix)
}
## ----------------------------------------------------------------------------------------------------------------------------------------
## Mixed-model test
## ----------------------------------------------------------------------------------------------------------------------------------------
model.fn <- function(data, Type = "ICC_Beta", test_var = NULL, model = "GEE", method = "del.eff"){
  
  data[[test_var]] <- factor(data[[test_var]])
  
  if(model == "GEE"){
    
    formu <- formula(paste(Type, "~ -1 + ", test_var))
    results <- geeglm(
      formu, 
      data = data,
      id = cluster,
      corstr = "exchangeable"
    )
    
  }
  
  if(model == "LMER"){
    
    formu <- formula(paste(Type, "~ -1 +", test_var, "+ (1 | cluster) + (1 | Gene) + (1 | seqnames)"))
    results <- lmer(
      formu, 
      data = data
    )
    
    
  }
  
  e <- emmeans(results, specs = test_var)
  c_df <- data.frame(Type = Type, contrast(e, type = "average", method = method))   # method = "del.eff"
  c_df$contrast <- gsub(" effect","",c_df$contrast)
  c_df$df <- nrow(data)
  c_df$p.t <- 2 * (1 - pt(abs(c_df$z.ratio), c_df$df))
  
  return(c_df)
  
}
## ----------------------------------------------------------------------------------------------------------------------------------------
## Add annotation
## ----------------------------------------------------------------------------------------------------------------------------------------
annotate_results <- function(result, dir.extra.aux){
  
  load(file.path(dir.extra.aux,"great_EPIC_array_annotation.rda"))
  load(file.path(dir.extra.aux,"E073_15_coreMarks_segments.rda"))
  data <- readr::read_tsv(
    file.path(dir.extra,"AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
  )
  CellType.selected <- readxl::read_xlsx(
    file.path(dir.extra, "Nassser study selected biosamples.xlsx"),col_names = FALSE
  ) %>% dplyr::pull(1)

  result$Islands.UCSC.Relation_to_Island <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC[result$cpg,"Relation_to_Island"]
  result$UCSC_RefGene_Name <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg,"UCSC_RefGene_Name"]       
  result$UCSC_RefGene_Group <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg,"UCSC_RefGene_Group"]     
  result$UCSC_RefGene_Group[result$UCSC_RefGene_Group == ""] <- "Intergenic"
  
  result <- dplyr::left_join(result, great, by = "cpg")

  data.filtered <- data %>% dplyr::filter(CellType %in% CellType.selected) %>% 
    dplyr::filter(!isSelfPromoter)  %>% 
    dplyr::filter(class != "promoter")
  
  nasser.enhancer.gr <- data.filtered %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "chr",
    keep.extra.columns = TRUE
  )
  result.gr <- result %>% makeGRangesFromDataFrame(
      start.field = "start",
      end.field = "end",
      seqnames.field = "seqnames"
    )
  hits <- findOverlaps(result.gr, nasser.enhancer.gr) %>% as.data.frame()
  result$enhancer <- 0
  result$enhancer[unique(hits$queryHits)] <- 1
  
  annotType <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Manifest
  
  annotTypeCompleteCol1 <- annotType[result$cpg, ]
  result$Probe_Type <- annotTypeCompleteCol1$Type
  
  return(result)
}
# ---------------------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------------------
summ.table <- function(df, Group, type = "ICC_Beta"){
  
  d <- df %>% group_by(get(Group)) %>% 
    dplyr::summarise(
      num_of_probes = n(),
      Mean = mean(get(type), na.rm = T),
      SD = sd(get(type), na.rm = T),
      Min = min(get(type), na.rm = T),
      Median = median(get(type), na.rm = T),
      Max = max(get(type), na.rm = T)
    ) 
  colnames(d)[1] <- Group
  
  return(d)  
}
# ---------------------------------------------------------------------------------------
# Calculate blood biomarker scores 
# ---------------------------------------------------------------------------------------
get_scores <- function(beta, coef_df){
  cpg_selected <- coef_df$CpG
  weights <- coef_df$Beta
  names(weights) <- cpg_selected
  
  common_cpg <- intersect(cpg_selected, rownames(beta))
  beta <- beta[common_cpg,]
  
  weights <- weights[common_cpg]
  
  scores <- as.numeric(weights %*% as.matrix(beta))
  names(scores) <- colnames(beta)
  
  scores
}

get_replicate_scores <- function(beta_list, coef_df){
  icc_df <- plyr::ldply(
    beta_list,
    .fun = function(ls){
      s <- get_scores(beta = ls, coef_df = coef_df)
      data.frame(barcodes = names(s), scores = s)
    }
  ) %>% column_to_rownames("barcodes") %>%
    t()
  
  colnames(icc_df) <- gsub("^X", "", colnames(icc_df))
  
  icc_df
}
