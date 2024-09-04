###############################################################################
## Code to compute the across module association between module median_MAGMA and z_kme, 
## accounting for module size using a robust linear model. This is used as an index of ANCR.


####################
###   PREPPING   ###
####################
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
library(sfsmisc)

rm(list = ls())
setwd(here::here())
subdir =  "ANCR_NetRemInMod"
basedir = "/set/my/dir/" 

dir.create(paste0(basedir , subdir, "/plots_new"), showWarnings = FALSE)
dir.create(paste0(basedir , subdir, "/results"), showWarnings = FALSE)
resultsdir = "/results"
setwd(paste0(basedir,subdir, resultsdir))
magma_quant_list = rev(as.list(as.character(seq(20,100,10))) %>% setNames(.,.))



#####Compute per-network ANCR index as the across-module rlm() of median MAGMA vs zKme accounting for module size
#### 
magma_quant_list_qq = imap(magma_quant_list, ~{
  magma_quant = ..2
  print(..2)
  # if ( !file.exists(paste0("allVars.results.processed",magma_quant,".qs")) ){
    ll = list.files(pattern = paste0(magma_quant,".qs")) %>% .[!grepl("processed",.)]
    ll = ll[!grepl("Li2018|Fromer",ll)]
    qq = map_dfr(ll,qread)
    qq1 = qq %>% group_nest(tissue_age, enrichment, list.type,simVar) %>% mutate(id = glue::glue("{tissue_age}..{enrichment}..{list.type}..{simVar}"))
    names(qq1$data) = qq1$id
    
    
    qq1$processed.data  =  qq1$data %>% imap(~{
      # print(.y)
      qq_small = unique(.x)
      quants = grep("median_H" , colnames(qq_small),value = T)
      ANCR = map_dfr(quants, ~{
        z_Kme = qq_small[["z_Kme"]]
        modSize = qq_small[["modSize"]]
        medianMAGMA = qq_small[[.x]]
        
        s2 = broom::tidy(rlm(medianMAGMA ~ z_Kme + modSize)  )  %>% .[grep("z_Kme"    ,.$term),] %>% pull(statistic)
        s22 <- f.robftest(rlm(medianMAGMA ~ z_Kme + modSize) , var = "z_Kme")[["p.value"]]
        
        return(c(s2,s22) %>%
                 set_names(c(
                             "rlm.median_H__z_kme",
                             "rlm.median_H__z_kme..pval")))
      }) %>% set_rownames(quants)
      
      return(list(median_H_stats = ANCR      ))
      
    }, .progress = T) 
    qq1
})

#####Transorm info long dataframe, assign trait types
#### 
ANCR_df = imap_dfr(magma_quant_list_qq, ~{
  qq1 = ..1
  test_median_H = future_map_dfr(qq1$processed.data, .id = "id",~{
    .x$median_H_stats %>% tibble::rownames_to_column("quant") %>% mutate(quant = as.numeric(gsub("median_H_Q","",quant)))
  }) %>% full_join(qq1[,!colnames(qq1) %in% c("data","processed.data")])
  net.PercentTopRiskRemoved = 100 - as.numeric(rep(..2, nrow(test_median_H)))
  out = cbind(test_median_H,net.PercentTopRiskRemoved)
})
ANCR_df$net.PercentTopRiskRemoved <- factor(ANCR_df$net.PercentTopRiskRemoved,
                                          levels = paste(sort(as.integer(unique(ANCR_df$net.PercentTopRiskRemoved)), decreasing = F)))

ANCR_df$Type = "None"
ANCR_df$Type[grepl("AD|PD|ALS",ANCR_df$simVar)] = "Neurological"
ANCR_df$Type[grepl("ASD|ADHD",ANCR_df$simVar)] = "Child Psychiatric"
ANCR_df$Type[grepl("RA|CD|UC",ANCR_df$simVar)] = "Immune"
ANCR_df$Type[grepl("OCD|SCZ|BIP|MDD",ANCR_df$simVar)] = "Adult Psychiatric"
ANCR_df$Type = factor(ANCR_df$Type, levels = c("Adult Psychiatric", "Child Psychiatric", 
                                               "Neurological", "Immune", "Environmental"))

ANCR_df$enrichment = gsub("BIP","BPD",ANCR_df$enrichment)
ANCR_df$mod.PercentTopRiskKept = (1/ANCR_df$quant)*100

#### Removing invalid networks and networks with fewer than 10 modules. Cleaning of superfluous columns
magma_quant_0rem = magma_quant_list_qq[["0"]]
qq1_nmods = map_int(unique(magma_quant_0rem$tissue_age)  %>%setNames(.,.), function(net){
  out = magma_quant_0rem %>% 
    filter(tissue_age == net ) %>% 
    filter(grepl("PGC3",list.type) )
  nrow(out[["data"]][[1]])
  # out = nrow(out[["data"]])
}) 
nets_rem = names(qq1_nmods)[qq1_nmods <10]
names(ANCR_df)[names(ANCR_df) == "tissue_age"] = "Network"
ANCR_df = ANCR_df %>% filter(!Network %in% nets_rem)
ANCR_df = ANCR_df %>% filter(!Network %in% c("DG","DG_QSVA", "HP_noQSVA"))
ANCR_df = ANCR_df[!grepl("sim",names(ANCR_df))]
ANCR_df = ANCR_df[!names(ANCR_df) == "quant"] 
ANCR_df = ANCR_df[!names(ANCR_df) == "list.type"] 
saveRDS(ANCR_df,"ANCR_NetRemInMod_df.rds")
