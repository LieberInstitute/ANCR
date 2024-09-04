#### Preprocess pre-pseudobulked snRNAseq data from Ruzicka et al. ####

#### Load dependencies ####
require(SummarizedExperiment)
require(tidyverse)

#### Mean counts data were downloaded https://www.synapse.org/#!Synapse:syn25922167 (pseudobulk_mean_logcounts_filtered.rds)  of Ruzicka et al ####
pseudobulk_mean_logcounts_filtered <- readRDS("../data/Ruzicka/pseudobulk_mean_logcounts_filtered.rds")

#### Get colData #### 
coldata =  as.data.frame(colData(pseudobulk_mean_logcounts_filtered))
coldata = cbind(Identifier = rownames(coldata),coldata)

#### Extract counts ####  
agg_counts_ctMats_logMean = pseudobulk_mean_logcounts_filtered@assays@data@listData

#### Same preprocessing as usual pre WGCNA #### 
per0_th = .2
agg_ctMats_logMean_prepro = imap(agg_counts_ctMats_logMean, function(mat,name) {
  subjects = colnames(mat)
  genes = rownames(mat)
  mat = as.matrix(mat)
  median.gene.exprs = matrixStats::rowMedians(mat)

  cat("\n",name, median(median.gene.exprs), 
      mean(median.gene.exprs), min(median.gene.exprs), max(median.gene.exprs),"\n")
  print(quantile(median.gene.exprs, probs = seq(0, 1, 0.1)))
  
  names(median.gene.exprs) = rownames(mat)
  gene_filter0 = as.vector(rowSums(mat ==0) <= ncol(mat)*per0_th)                #Genes to keep (not more than 20% samples are zero)
  gene_filter1 = as.vector(median.gene.exprs   >= 0)                              #Genes to keep (median exp > 0.1)
  
  median.gene.exprs_filtered = median.gene.exprs[gene_filter0 & gene_filter1]
  gene_filter2 = as.vector(abs(scale(median.gene.exprs_filtered)) <  3)            #Genes to keep (z-scoreed_median_exp < 3)
  genes_to_keep  = names(median.gene.exprs_filtered)[gene_filter2]
  
  mat = mat[rownames(mat) %in% genes_to_keep, ]
  mat = mat[!grepl("^MT-",rownames(mat)), ]
  mat = as.data.frame(mat)
  list(logMean = mat)
}, .progress = T)


n_genes = map_int(agg_ctMats_logMean_prepro, function(mat) {
  nrow(mat[["logMean"]]) +1
})
jpeg(paste0("plots/preprocess/n_genes_Ruzicka","_per0th_", as.character(per0_th),".jpeg"), height = 1000, width = 1000)
barplot(n_genes ,main="Number of genes after prepro",ylab="Freqency",cex.names = .5,las=2)
dev.off()



gene_th = 5000
agg_ctMats_logMean_prepro = agg_ctMats_logMean_prepro[map_lgl(agg_ctMats_logMean_prepro, ~{
  nrow(..1$logMean) > gene_th
  })]


#### Same preprocessing as usual pre WGCNA #### 
#### remove MB8  (technical replicate of MB8-2)#### 	
library(sva)
library(jaffelab)
set.seed(123)
library(limma)

facts = names(coldata)
scree_fp = "plots/scree_PBct_ruzicka/"
dir.create(scree_fp)
PCheat_fp = "plots/PCheat_PBct_ruzicka/"
dir.create(PCheat_fp)

mat.name = "Ex-L3"
mat_list = agg_ctMats_logMean_prepro[[mat.name]]
protect = 1
lets = paste(paste(letters, collapse ="|"),paste(LETTERS, collapse ="|"), sep ="|")
conf_possible = names(coldata)
conf = grep("Age|PMI|Ancestry|Gender|mito_perc|Batch|Phenotype",conf_possible,value = T)
get.cleaned.obj = function(mat_list, mat.name, protect){
  mat = mat_list$logMean
  coldata_cr = coldata[coldata$Identifier %in% colnames(mat) ,]
  coldata_cr = coldata_cr[match(colnames(mat),coldata_cr$Identifier),]
  if(!identical(coldata_cr$Identifier,colnames(mat))) {
    stop("Mat and coldata not identical")
  }
  # cat("Mat and coldata identical?", identical(coldata_cr$Identifier,colnames(mat)))
  pcobj = prcomp(t(mat),center = T, scale. = T)
  if(identical(coldata_cr$Identifier,rownames(pcobj$x))) {
    coldata_cr = cbind.data.frame(coldata_cr,pcobj$x)
  } else {stop("pcobj$x and coldata not identical")}
  var_explained = pcobj$sdev^2 / sum(pcobj$sdev^2)
  jpeg(paste0(scree_fp,mat.name,".jpeg"))
  # plot(1:length(var_explained), var_explained)
  print(qplot(1:length(var_explained), var_explained) + 
          geom_line() + 
          xlab("Principal Component") + 
          ylab("Variance Explained") +
          ggtitle(mat.name) +
          ylim(0, 1))
  dev.off()
  
  coldata_cr_scale.factor = imap_dfc(coldata_cr, function(c,name) {
    if(is.numeric(c)) {
      # cat(name)
      # c = gsub("\\,","\\.",c)
      c = scale((c))[,1]
    } else {
      c = as.factor(c)
    }
  })
  
  check = map_dfc(conf %>% setNames(.,.), function(x) {
    round(cor(pcobj$x,as.numeric(coldata_cr_scale.factor[[x]]))*100,1)
  })
  
  check = cbind(PCimportance = t(summary(pcobj)$importance["Proportion of Variance",,drop=F])*100,
                 check)
  
  jpeg(paste0(PCheat_fp,mat.name,".jpeg"), height = 1200, width = 350)
  pheatmap::pheatmap(check,cluster_rows = F, cluster_cols = F, cellwidth = 20)
  dev.off()
  
  form = reformulate(c(conf, "PC1"))
  mod = model.matrix(~Age + Phenotype + Batch + Gender + PMI + 
                       EUR_Ancestry +
                       EAS_Ancestry +
                       AMR_Ancestry +
                       SAS_Ancestry + AFR_Ancestry +
                       mito_perc +
                       PC1, data = coldata_cr_scale.factor)
  
  print(mat.name)
  
  #Check of multi-colinearity between known covariates and PCs
  test = alias(lm(t(mat)~mod))$Complete
  cat(sum(alias(lm(t(mat)~mod))$Complete))
  
  cleaned_matrix = cleaningY(mat, mod = mod, P = protect)
  
  all.equal(cleaned_matrix,
            mat)
  
  check.shapiro = apply(cleaned_matrix,1, function(r) shapiro.test(r)$p.value)
  print(paste0("Non-normal genes (p <0.05) before ranknorm: ", sum(check.shapiro <0.05), " / ",length(check.shapiro)))
  
  exp_ranknorm           = t(apply(cleaned_matrix, 1, RNOmni::RankNorm))  #rankNorm each gene in the matrix
  # assays(rse)$ranknorm   = exp_ranknorm
  
  mat_list$svacleaned = cleaned_matrix
  mat_list$ranknorm = as.data.frame(exp_ranknorm)
  
  return(mat_list)
}

agg_ctMats_logMean_prepro_sva = imap(agg_ctMats_logMean_prepro, function(mat_list,name) {
  get.cleaned.obj(mat_list = mat_list, mat.name = name, protect = 1)
}, .progress = T)
saveRDS(agg_ctMats_logMean_prepro_sva,file =  "results/agg_mats_Ruzicka/Ruzicka_aggctMats_CPM_prepro_sva.rds")

