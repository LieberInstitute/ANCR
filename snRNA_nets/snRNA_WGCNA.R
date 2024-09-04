#### Compute networks per celltype in a snRNAseq dataset. ####
## First identify the minimum soft power beta (B) that allows for matched median connectivity across all celltype networks 
## Use said B power to construct co-expression networks.
## Example for Batiuk. Run for Ruzicka by uncommenting lines 10-12. 

#### Pseudobulked expression matrices for Batiuk (constructed in Batiuk_pseudobulk.R)
mat_list = readRDS("data/Batiuk_aggctMats_prepro_sva.rds")
fp_recas = "WGCNA/Batiuk/"

# #### Pseudobulked expression matrices for Ruzicka (constructed in Ruzicka_pseudobulk.R)
# mat_list = readRDS("data/Ruzicka_aggctMats_prepro_sva.rds")
# fp_recas = "WGCNA/Ruzicka/"

#### Load dependencies ####
suppressPackageStartupMessages({
  library(dynamicTreeCut)
  library(fastcluster)
  library(WGCNA)
  library(SummarizedExperiment)
  library(qs)
  library(dendextend)
  library(limma)
  library(pheatmap)
  library(gtools)       #Optional. Just to sort path of files read
  library(furrr)
  library(doParallel)
  library(WGCNA)
  library(tidyverse)
})

set.seed(123)
CPU = detectCores() - 2
registerDoParallel(cores=CPU)
getDoParWorkers()
options(mc.cores = CPU, future.globals.maxSize= 3565158400)
# plan(multicore, workers = CPU )

fp_sft.pwrs= paste0(fp_recas,"results/sft.pwrs/") 
dir.create(fp_sft.pwrs)
setwd(fp_sft.pwrs)


### Get B pwrs per celltype ####
get_sft_parameters = function(ct_mats, mod, net.type, assay.name, nname){
  cat(paste0(nname,"...sft.",mod,".",net.type), "- ")
  input_matrix = t(ct_mats$ranknorm) #Genes in columns, subjects in rows
  
  if(mod == "pearson"){
    li_more = list(
      corOpt = list(nThreads = CPU)
      #Not specifying corFn defaults to pearson
    )
  }
  
  WGCNA::allowWGCNAThreads(nThreads = CPU)
  WGCNAnThreads()
  
  li.args = list(data = input_matrix, powerVector = 1:30, RsquaredCut = 0.8,verbose = 5, blockSize = 20000, moreNetworkConcepts = TRUE)
  
  ret     = do.call(WGCNA::pickSoftThreshold, c(li.args, li_more, networkType = net.type))
  saveRDS(ret, paste0(nname,"...sft.",mod,".",net.type,".rds"))
  gc()
  return(ret)
}
assay.name = "ranknorm"
pwr.type   = "corr"
mod        = "pearson"
net.type = "signed hybrid"
sft.pearson.shybrid   =  imap(mat_list,
  ~ get_sft_parameters(.x, nname = .y, mod = mod , 
    assay.name = assay.name, net.type = "signed hybrid"))
saveRDS(sft.pearson.shybrid, "sft.pearson.shybrid.rds"))



#### selected the parameter B such that connectivity was matched across all networks (by_all)  ####
dir      = getwd()
ls       = gtools::mixedsort(list.files(dir, pattern = ".rds"))
obj.name = substring(text = ls,first = 1,last = nchar(ls)-4)
fp       = file.path(dir,ls)
# fp       = fp[grepl("sACC|Amygdala",fp)]
li = list()
for (f in seq_along(fp)){
  li[[obj.name[f]]] = readRDS(fp[f])
}
li = li[!grepl("soft-thresholding-pwrs",names(li))]

li_fitIndices                = map(li,"fitIndices")           ##List of sft dataframe 

metadf = tibble(full.name = names(li))
metadf = metadf %>% mutate(
  Celltype     = strsplit2(full.name,"\\.\\.\\.")[,1]                                ,
  correlation    = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,2]    ,
  network        = strsplit2(strsplit2(full.name,"\\.\\.\\.")[,2],"\\.")[,3]    ,
  # samples        = n.samples[Type]                                        ,
  powerEstimate  = unlist(map(li,"powerEstimate"))[.$full.name]                 ,
  class          = ifelse(grepl("gr",full.name),"parsed","non-parsed")          ,
  nl = map(li_fitIndices, ~{
    scalefreeIndex = -sign(.x$slope) * .x$SFT.R.sq                                   
    
    if(any(scalefreeIndex > 0.8)){
      allowed.pwrs   = .x$Power[scalefreeIndex > 0.8]                   
      allowed.conn   = .x$median.k.[scalefreeIndex > 0.8]
      
    } else {                              #No pwr makes network scalefree
      allowed.pwrs   = .x$Power[which.min(abs(scalefreeIndex - 0.8))]                   
      allowed.conn   = .x$median.k.[which.min(abs(scalefreeIndex - 0.8))]
    }
    
    return(dplyr::lst(allowed.pwrs, allowed.conn))
  })) %>% unnest_wider(nl) %>% arrange(correlation, network)
# 
mm = metadf
mm = mm %>% mutate(
  by_individual = pmap(list(allowed.pwrs, allowed.conn) ,~{
    pwr = min(.x)
    conn = round(.y[which(.x == pwr)],3)
    return(dplyr::lst(pwr,conn))
  }),
  notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}) %>% unnest_wider(by_individual,names_sep = ".")
mm = mm [!is.na(mm$powerEstimate),]

mm = mm %>% group_by(correlation,network) %>% 
  mutate(by_connLIBD = { 
    th.conn = 1.123
    # min.conn = LIBD_allowedconn
    index = map(allowed.conn,~ which.min(abs(. - th.conn)))
    
    pmap(list(allowed.pwrs,allowed.conn,index, th.conn),~{
      pwr = ..1[..3]
      conn = round(..2[..3],3)
      return(dplyr::lst(pwr,conn, matched.conn = ..4))
    })},
    notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
  )  %>% unnest_wider(by_connLIBD,names_sep = ".") %>% ungroup()

mm = mm %>% group_by(correlation,network) %>% 
  mutate(by_conn10 = { 
    th.conn = 10
    # min.conn = LIBD_allowedconn
    index = map(allowed.conn,~ which.min(abs(. - th.conn)))
    
    pmap(list(allowed.pwrs,allowed.conn,index, th.conn),~{
      pwr = ..1[..3]
      conn = round(..2[..3],3)
      return(dplyr::lst(pwr,conn, matched.conn = ..4))
    })},
    notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
  )  %>% unnest_wider(by_conn10,names_sep = ".") %>% ungroup()

mm = mm %>% group_by(correlation,network) %>% 
  mutate(by_conn30 = { 
    th.conn = 30
    # min.conn = LIBD_allowedconn
    index = map(allowed.conn,~ which.min(abs(. - th.conn)))
    
    pmap(list(allowed.pwrs,allowed.conn,index, th.conn),~{
      pwr = ..1[..3]
      conn = round(..2[..3],3)
      return(dplyr::lst(pwr,conn, matched.conn = ..4))
    })},
    notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
  )  %>% unnest_wider(by_conn30,names_sep = ".") %>% ungroup()

mm = mm %>% group_by(correlation,network) %>% 
  mutate(by_all = { 
    min.conn = min(map_dbl(allowed.conn,max))
    index = map(allowed.conn,~ which.min(abs(. - min.conn)))
    
    pmap(list(allowed.pwrs,allowed.conn,index, min.conn),~{
      pwr = ..1[..3]
      conn = round(..2[..3],3)
      return(dplyr::lst(pwr,conn, matched.conn = ..4))
    })},
    notScalefree = {ifelse(anyNA(powerEstimate),"yes","no")}
  )  %>% unnest_wider(by_all,names_sep = ".") %>% ungroup()

mm = mm  %>% arrange(correlation, network)        

saveRDS(mm, "soft-thresholding-pwrs.rds")



#### Compute net per celltype using "by_all" B power ####
readRDS("soft-thresholding-pwrs.rds")
get_networks = function(ct_mats,mod, net.type, assay.name, beta.pwr,pwr.type, nname){
  if(mod == "pearson"){
    li_more = list(nThreads = CPU, corType = "pearson")
    #Not specifying corFn defaults to pearson
  }
  
  WGCNA::allowWGCNAThreads(nThreads = CPU)
  WGCNAnThreads()
  
  pwr = unlist(beta.pwr[beta.pwr$Celltype == nname & beta.pwr$correlation == mod & beta.pwr$network == net.type, paste0(pwr.type,".pwr")])
  
  cat(paste0(nname,"...sft.",mod,".",net.type), "- ")
  input_matrix = t(ct_mats$ranknorm) #Genes in columns, subjects in rows

  li.args = list(datExpr = input_matrix, randomSeed = 123, verbose = 5, maxBlockSize = 30000,
                 TOMType = "signed", saveTOMs = FALSE, saveTOMFileBase = (nname), power = pwr,
                 minModuleSize = 50,
                 pamStage = TRUE, pamRespectsDendro = TRUE,  deepSplit = 4,                   
                 reassignThreshold = 1e-06,                                                   
                 mergeCutHeight    = 0.1,                                                    
                 numericLabels     = TRUE)
  
  print(nname)
  print(pwr)
  cor_default    = cor
  cor            = WGCNA::cor      #Overwrite default cor method as WGCNA::cor instead of stats::cor
  net            = do.call("blockwiseModules", c(li.args, li_more, networkType = net.type))
  
  li.args.IMC    = list(datExpr =  input_matrix, colors = net$colors, power = pwr, getWholeNetworkConnectivity = TRUE)
  IMC            = do.call("intramodularConnectivity.fromExpr", c(li.args.IMC, networkType = net.type)) 
  IMC$modules    = net$colors 
  rownames(IMC)  = names(net$colors)  
  net$IMC        = IMC
  
  net$KME        = WGCNA::signedKME(input_matrix, net$MEs)
  cor            = cor_default     #Reassign default cor method
  
  
  print(paste0(nname,"...sft.",mod,".",net.type,"|pwr:",pwr,"|pwr.type: ",pwr.type))
  saveRDS(net, paste0(nname,"...sft.",mod,".",net.type," & pwr-",pwr," & pwr.type-",pwr.type, ".rds"))
  #gc()

  adj = adjacency(datExpr = input_matrix, type = net.type, power = pwr)
  qsave(adj,  paste0(fp_recas,"results/adj","/adj_",nname, ".qs") ,nthreads = 50)
  rm(adj)

  return(net)
}

pwr.df <- readRDS("soft-thresholding-pwrs.rds")
assay.name = "ranknorm"
pwr.type   = "by_all"
mod        = "pearson"
sdir = paste("net",assay.name,pwr.type,mod,sep= "__")

fp_nets = paste0(fp_recas,"results/nets/")
dir.create(fp_nets)
dir.create(fp_nets %>% paste0(.,sdir))
setwd(fp_nets %>% paste0(.,sdir))
#Run "WGCNA::blockwiseModules" to get network
#pearson correlation on the ranknorm assay, signed hybrid network, "by_all" power
net.pearson.shybrid   =  imap(mat_list,~ get_networks(
  .x, nname = .y, mod = mod , 
  assay.name = assay.name, beta.pwr = pwr.df, 
  pwr.type = pwr.type, net.type = "signed hybrid"))

