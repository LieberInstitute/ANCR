###############################################################################
## This code provides results for Figure 2, Figure 3AB and Figure 4. 
## 
## Code to compute network-wide risk convergence onto a module (z_kme),
## for all modules(n_genes>=30) in a set of  networks. Median MAMGMA is also 
## computed, considering different proportions of top module MAGMA genes. z_kme is
## computed removing different proportions of top MAGMA genes in the network,
## defined by the "magma_quant" variable. 
## 
## Then in compute_ANCR_step2.R the across module association between module median_MAGMA and z_kme, 
## accounting for module size using a robust linear model...... from rlm(median_MAGMA ~ z_kme + module_size) 
## This across module association is used as an index of ANCR
##  



####################
###   PREPPING   ###
####################
# LOAD LIBRARIES
library(limma)
library(dplyr)
library(magrittr)
library(tidyr)
library(purrr)
library(furrr)
library(qs)
library(MASS)
library(geomtextpath)
library(ggplot2)

rm(list = ls())
setwd(here::here())
subdir =  "ANCR_RemNetInNet"
dir.create(paste0("test_sets/",subdir), showWarnings = FALSE)
dir.create(paste0("test_sets/",subdir, "/results"), showWarnings = FALSE)
# dir.create(paste0("test_sets/",subdir, "/plots"), showWarnings = FALSE)
options(stringsAsFactors = FALSE)
seed <- 1997
set.seed(seed)

####### Function modulewise across-network risk convergence ----
## This function runs for each module in the network, 
## and computes the risk convergence strength onto that module based on the 
## across-gene association between MAGMA and kme to said module
modulewise.fn = function(currentModule , allModules, 
                         df_small_quantcrop,df_small_nocrop, 
                         validConfounders_quantcrop, validConfounders_nocrop, 
                         confoundMod_quantcrop,confoundMod_nocrop,  
                         Mcol = Mcol, KMEcols = KMEcols, KTOTcol = KTOTcol , KWITcol = KWITcol,  
                         zVars, simVars = NULL, networkName = networkName){
  
  
  #subTMP_og is data for genes in module, NOT removing genes based on network wide MAGMA rank (not considering magma_quant cutoff)
  subTMP_og   <- df_small_nocrop %>% filter(.data[[Mcol]] %in% currentModule) %>%
    {if (length(KWITcol)>0) {
      filter(., is.finite(.data[[KWITcol]]))
    } else {
      .
    }
    }
  subTMP_og = subTMP_og[complete.cases(subTMP_og),]
  if(nrow(subTMP_og) < 10) return(NULL)
  
  
  #The KME for current module
  ### otherTMP is data for genes outside module, removing genes based on network wide MAGMA rank ( considering magma_quant cutoff)
  ### Filter out genes with NA/NaN for current module KME column
  currentKMEcol = paste0(networkName,".KME",currentModule) %>% .[. %in% colnames(df_small_quantcrop)]
  otherTMP <- df_small_quantcrop %>% filter(!.data[[Mcol]] %in% currentModule) %>%
    {if (length(currentKMEcol)>0) {
      filter(., is.finite(.data[[currentKMEcol]]))
    } else {
      .
    }
    }
  
  ### Compute convergence of risk onto this module for each type of MAGMA
  allResults = map_dfr(c(zVars),~{
    currentVar = .x
    
    ## For Figure 2 and Figure 3AB, z_kme is computed considering the previous network wide gene removal based on magma_quant
    ## thus otherTMP is used
    otherTMP_var = otherTMP[!is.na(otherTMP[currentVar]),]
    
    ## Figure 3CD ## For Figure 3CD, z_kme is computed NOT considering the previous network wide gene removal based on magma_quant, thus otherTMP_og is used
    ## Figure 3CD # otherTMP_var = otherTMP_og[!is.na(otherTMP_og[currentVar]),]
    
    #If any confounder does not have enough levels in otherTMP, remove 
    #If any confounder has NA values in otherTMP, remove 
    moduleconfounders = validConfounders_quantcrop[validConfounders_quantcrop %>% map_lgl(~{
      (length(unique(as.character(otherTMP_var[,.x]))) > 1) & !any(is.na(as.character(otherTMP_var[,.x])))
    })]
    
    ## The lm association between MAGMA and module kme is computed considering present confounders
    if (length(currentKMEcol)>0) {
      otherKmeMod      = lm(reformulate(response = currentKMEcol, termlabels = c(moduleconfounders,currentVar)), data = otherTMP_var)
      otherKmeCoeff    = summary(otherKmeMod)$coefficients
      #t_Kme            = otherKmeCoeff[grep("ZSTAT",rownames(otherKmeCoeff)), "t value"]
      t_Kme            = otherKmeCoeff[currentVar, "t value"]
      residual_Kme     = otherKmeMod$df.residual
      z_Kme            = zscoreT(t_Kme    , residual_Kme    , approx=FALSE, method = "hill")
    } else {
      otherKmeMod      = NA
      otherKmeCoeff    = NA
      t_Kme            = NA
      residual_Kme     = NA
      z_Kme            = NA
    }
    
    ## Here the median MAGMA is computed for said module considering different proportions of top module MAGMA genes.
    ## These different median MAGMA values will give an idea of what proportion
    ## of module MAGMA genes contribute to the prediction of network wide risk convergence onto that module
    ## This is depicted in Figure 2, "% of top risk score genes considered in module to compute its risk enrichment" on the x axis ,
    ## and resulting ANCR score on the y axis (final across module ANCR computation is shown in "compute_ANCR_step2_xModule.R"). 
    allQuant <- c(1000L, 100L, 20L, 10L, 4L, 2L, 1L)
    median_H <- allQuant %>% set_names(paste0("median_H_Q",.)) %>% map_dfc(~{
      whichQuant = .x
      
      
      ## For Figure 2, 3AB and 4, this value is computed NOT considering the previous 
      ## network wide gene removal based on magma_quant (considering all initial netowrk genes)
      residuals(confoundMod_nocrop)[rownames(subTMP_og),currentVar]%>%
        sort(decreasing = TRUE)                %>%
        head(ceiling(nrow(subTMP_og)/whichQuant)) %>%
        median()
      
      
    })
    
    tmp.df = data.frame(simVar = currentVar, t_Kme, residual_Kme, z_Kme, median_H)
  })
  
  allResults = allResults %>% mutate(modules = currentModule, modSize = nrow(subTMP_og))
  
  #out.df = data.frame(modules = currentModule, modSize = nrow(subTMP), t_kWithin, t_Kme, median_H, residual_kWithin, residual_Kme, z_kWithin, z_Kme) %>% tibble::rownames_to_column("ZSTAT")
  
}
#############################



#Main----

#Read dataframe containing all gene network information
df <- readRDS("../data/wideforms_WGCNA/gene_net_info.rds")
df_extra <- readRDS("../data/wideforms_WGCNA/gene_net_info_extra.rds")
df = merge(df,df_extra)
colnames(df)[colnames(df)=="Werling2020|silver"] = "Werling2020|grey60"
rownames(df) <- df$ensemblID
colnames(df) <- gsub(" ", "_", colnames(df))
df[ ,sapply(df, function(x) all(is.na(x)))] <- list(NULL) #Remove a column if all values in column are NA

#choose initial confounders
confounders = c("seqnames","start","strand","gene_type","NumTx","Length","GC_content")

#Define networks to evaluate, here we evaluate only bulk tissue networks
allNet_og =  gsub(".modules", "",grep("modules",names(df), value = T))
allNet = allNet_og
# allNet = allNet_og[!grepl("Ruz_|Bat_",allNet_og)]

#Define MAGMA traits to evaluate
zVars   <- grep("35.10", colnames(df), value = T)
zVars   <- zVars[!grepl("PGC2|.eur.|CUD|BARI|LIBD|SA\\.model", zVars)]
zVars_select = zVars[grepl("SCZ.PGC3.kbp35.10.ZSTAT|AD.|UC.|CD.|ALS.|BIP.|MDD.|ASD.|ADHD.|PD.|ASD.|RA", zVars)] 

plan(multisession,workers = 6)

setwd(paste0("test_sets/" , 
             subdir, "/results"))



#####################
###   INITIALIZE   ###
######################

## Here the df will be cropped based on the a network-wide MAGMA cut off.
## magma_quant defines the percentage of lower MAGMA genes that will be kept in 
## the network to compute the z_kme (this code is for Figure 2 and Figure 3AB)
## For Figure 3CD magma_quant is used to compute module median_MAGMA, and z_kme is computed always considering all genes, 
## changes to this end are commentented in appropriate positions, identified by a "## Figure 3CD #" margin .
for(magma_quant in c("100%", "90%","80%", "70%", "60%", "50%","40%", "30%", "20%")) {
  discard.var = map(allNet,~{ #Iterate for each network
    networkName = .x
    map(as.list(zVars_select), ~ {
      zvar = .x
      Mcol    <- grep(paste0(networkName,".modules" ), colnames(df), value = T)
      KMEcols <- grep(paste0(networkName,".KME"     ), colnames(df), value = T)
      KTOTcol <- grep(paste0(networkName,".kTotal"  ), colnames(df), value = T)
      KWITcol <- grep(paste0(networkName,".kWithin" ), colnames(df), value = T)
      
      
      df_small = df[,c("ensemblID",confounders, 
                       zVars,    #Outcome variable (pathology zscores)
                       Mcol, KMEcols, KTOTcol, KWITcol)] %>% set_rownames(.$ensemblID)
      df_small = df_small %>% filter(!is.na(.data[[zvar]]))
      df_small = df_small %>% filter(!is.na(.data[[Mcol]]))
      df_small = df_small %>% filter(gene_type %in% names(table(gene_type)[table(gene_type,useNA = "ifany")>40]))
      df_small = df_small %>% mutate(across(where(is.character), factor))
      df_small[,KTOTcol] <- log(df_small[,KTOTcol]+.0001)
      df_small[,KWITcol] <- log(df_small[,KWITcol]+.01)
      df_small[[Mcol]]  = unlist(df_small[[Mcol]] )
      df_small[[Mcol]] == NULL
      
      rownames(df_small) = df_small$ensemblID
      
      confounders = confounders[map_lgl(confounders,~{
        (length(unique(as.character(df_small[,.x]))) > 1) & !any(is.na(as.character(df_small[,.x])))})]
      responseVars = paste0("cbind(",paste0(c(zVars),collapse = ","),")")
      confoundMod <- lm(reformulate(response = responseVars, termlabels = confounders),data = df_small)
      residuals_zvar = residuals(confoundMod)[,zvar]
      
      df_small = df_small[rownames(df_small) %in% names(residuals_zvar),]
      
      n_quant = 100
      quants = quantile(residuals_zvar, probs = seq(0, 1, 1/n_quant), 
                        na.rm = T,
                        names = TRUE)
      df_small_quantcrop = df_small[residuals_zvar <= quants[magma_quant],]
      confounders_quantcrop = confounders[map_lgl(confounders,~{
        (length(unique(as.character(df_small[,.x]))) > 1) & !any(is.na(as.character(df_small[,.x])))})]
      confoundMod_quantcrop <- lm(reformulate(response = responseVars, termlabels = confounders_quantcrop),data = df_small_quantcrop)
      
      df_small_nocrop = df_small
      confoundMod_nocrop <- lm(reformulate(response = responseVars, termlabels = confounders),data = df_small_nocrop)
      
      min_size <- 30
      #Find all modules for current network and filter out grey and small modules
      modules = df_small %>% dplyr::count(.data[[Mcol]]) %>% filter(n>min_size & .data[[Mcol]]!= "grey") %>% pull(.data[[Mcol]]) %>% as.character()
      
      ANCR_DF <- purrr::map_dfr(modules,~ {modulewise.fn(currentModule = .x, allModules = modules, 
                                                          df_small_quantcrop = df_small_quantcrop,df_small_nocrop = df_small_nocrop, 
                                                          validConfounders_quantcrop= confounders_quantcrop, validConfounders_nocrop = confounders, 
                                                          confoundMod_quantcrop = confoundMod_quantcrop, confoundMod_nocrop = confoundMod_nocrop, 
                                                          Mcol = Mcol, KMEcols = KMEcols, KTOTcol = KTOTcol , KWITcol = KWITcol,
                                                          zVars = zvar, simVars = NULL, networkName = networkName)})
      ANCR_DF = ANCR_DF %>% mutate(tissue_age = networkName, enrichment = ifelse(grepl("ZSTAT",simVar),strsplit2(simVar,"\\.")[,1],"SCZ"), list.type = ifelse(enrichment == "SCZ","PGC3",enrichment),.before =1)
      final = ANCR_DF
      qsave(final,paste0(networkName,zvar,"_quant",gsub("%", "" ,magma_quant),".qs"))
      
    })
    
    gc()
  }, .progress = magma_quant)  
}
