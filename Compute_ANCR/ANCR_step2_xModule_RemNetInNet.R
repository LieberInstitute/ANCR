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
subdir =  "ANCR_RemNetInNet"
basedir = "/set/my/dir/" 
dir.create(paste0(basedir , subdir, "/plots_new"), showWarnings = FALSE)
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


ct_annots = list(
  Ex_L2 = c("Ex.L2"  ,    "Ex.L23", "_L2_" ),
  # L3 = c("Ex.L2"  ,    "Ex.L23", "Ex.L4", "L3"),
  Ex_L3 = c("Ex.L3", "_L3_"),
  Ex_L4 = c("Ex.L45" ,    "Ex.L4" , "_L4_"),
  # L5 = c("Ex.L45" , "Ex.L56"  ,   "Ex.L5b" , "L5"),
  Ex_L5 = c( "Ex.L56"  ,   "Ex.L5b" , "_L5_"),
  Ex_L6 = c("Ex.L6","Ex.L6b", "_L6_" ),
  # L56 = c(c( "Ex.L56"  ,   "Ex.L5b" , "_L5_"), c("Ex.L6","Ex.L6b", "_L6_" )),
  In_PV = c("In.PV","PVALB"),
  In_SST = c("In.SST","SST"),
  In_VIP = c("In.VIP","VIP"),
  In_ID2 = c("ID2"),
  In_INH_OTHER = c("_Inhib_1","_Inhib_5"),
  In_REELIN = c("In.Reelin"),
  In_ROSEHIP = c("In.Rosehip"),
  GLIAL = c("Oli" ,"Ast" , "Glia", "Oligos"),
  OPC = c("OPC" ),
  OTHER = "Other"
)

ct_superannots = list(
  Excit = "^L",
  Inhib = "PV|SST|VIP|ID2|REELIN|ROSEHIP",
  Glia = "GLIAL",
  OPC = c("OPC" ),
  OTHER = "OTHER"
)


ANCR_df$Type = "None"
ANCR_df$Type[grepl("AD|PD|ALS",ANCR_df$simVar)] = "Neurological"
ANCR_df$Type[grepl("ASD|ADHD",ANCR_df$simVar)] = "Child Psychiatric"
ANCR_df$Type[grepl("SA.|CUD|PTSD",ANCR_df$simVar)] = "Environmental"
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

ct_annots = list(
  Ex_L2 = c("Ex.L2"  ,    "Ex.L23", "_L2_" ),
  # L3 = c("Ex.L2"  ,    "Ex.L23", "Ex.L4", "L3"),
  Ex_L3 = c("Ex.L3", "_L3_"),
  Ex_L4 = c("Ex.L45" ,    "Ex.L4" , "_L4_"),
  # L5 = c("Ex.L45" , "Ex.L56"  ,   "Ex.L5b" , "L5"),
  Ex_L5 = c( "Ex.L56"  ,   "Ex.L5b" , "_L5_"),
  Ex_L6 = c("Ex.L6","Ex.L6b", "_L6_" ),
  # L56 = c(c( "Ex.L56"  ,   "Ex.L5b" , "_L5_"), c("Ex.L6","Ex.L6b", "_L6_" )),
  In_PV = c("In.PV","PVALB"),
  In_SST = c("In.SST","SST"),
  In_VIP = c("In.VIP","VIP"),
  In_ID2 = c("ID2"),
  In_INH_OTHER = c("_Inhib_1","_Inhib_5"),
  In_REELIN = c("In.Reelin"),
  In_ROSEHIP = c("In.Rosehip"),
  GLIAL = c("Oli" ,"Ast" , "Glia", "Oligos"),
  OPC = c("OPC" ),
  OTHER = "Other"
)

ct_superannots = list(
  Excit = "^L",
  Inhib = "PV|SST|VIP|ID2|REELIN|ROSEHIP",
  Glia = "GLIAL",
  OPC = "OPC" ,
  OTHER = "OTHER"
)


ANCR_df$data_source = stringr::str_split_fixed(ANCR_df$Network, pattern = "__", n = 2)[,1]
ANCR_df$ct = stringr::str_split_fixed(ANCR_df$Network, pattern = "__", n = 2)[,2]
ANCR_df$ct_prefix = stringr::str_split_fixed(ANCR_df$ct, pattern = "_", n = 2)[,1]
ANCR_df$ct_group =NA
for(ct_name in names(ct_annots)) {
  ct_v =  ct_annots[[ct_name]]
  ANCR_df$ct_group[grepl(paste(ct_v, collapse = "|"),ANCR_df$Network)] = ct_name
}
ANCR_df$ct_supergroup = NA
for(ct_name in names(ct_superannots)) {
  ct_v =  ct_superannots[[ct_name]]
  ANCR_df$ct_supergroup[grepl(ct_v,ANCR_df$ct_group)] = ct_name
}
ANCR_df = ANCR_df[!ANCR_df$net %in% nets_rem,]
#### celltypes removed if not in both ruz and bat after removing nets with fewer than 10 modules
ct_rem = c("Ex_L6","In_PV","In_ID2","In_REELIN", "In_ROSEHIP", "OPC", "OTHER")
ANCR_df = ANCR_df %>% filter(!ct_group %in% ct_rem)


ANCR_df_BULKsumm = ANCR_df %>% 
  filter(!grepl("Bat_|Ruz_", Network)) %>% 
  group_by(enrichment,net.PercentTopRiskRemoved,mod.PercentTopRiskKept) %>% 
  summarise(
    median_ANCRz = median(rlm.median_H__z_kme),
    median_ANCRp = median(rlm.median_H__z_kme..pval)
  )


ANCR_df_scRNAsumm = ANCR_df %>% 
  filter(grepl("Bat_|Ruz_", Network)) %>% 
  group_by(enrichment,net.PercentTopRiskRemoved,mod.PercentTopRiskKept,ct_group) %>% 
  summarise(
    median_ANCRz = median(rlm.median_H__z_kme),
    median_ANCRp = median(rlm.median_H__z_kme..pval)
  )
