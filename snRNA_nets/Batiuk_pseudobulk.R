#### Pseudobulk and preprocess snRNAseq data from Batiuk et al., 2022 ####

#### Load dependencies ####
require(readxl)
require(tidyverse) 

### Split original snRNA data per subject to lighten #### 
snRNA_seq_raw_countmatrices <- readRDS("../data/Batuik/snRNA-seq_and_spatial_transcriptomics/snRNA-seq_raw_countmatrices.RDS")
fp = "../data/Batuik/snRNA_seq_raw_countmatrices/"
split_snRNA = function() {
  dir.create(fp)
  iwalk(snRNA_seq_raw_countmatrices, function(subject_data, subject_name) {
    saveRDS(subject_data, paste0(fp,subject_name, ".rds" ))
  }, .progress = T)
}
split_snRNA()


#### Choose ct annotation file and get ct names #### 
annotations_final <- readRDS("../data/Batuik/snRNA-seq_and_spatial_transcriptomics/annotations_final.RDS")
annot = annotations_final$high
cts = names(table(annot)) %>% setNames(.,.)
coldata =  read_excel("../data/Batuik/sciadv.abn8367_tables_s1_to_s6.xlsx", sheet = "Table 1 ", skip = 1)
names(coldata)


#### Function to get ct pseudobulk means per subject mat #### 
get_ctAgg_mats = function(mat) {
  nuclei = colnames(mat)
  genes = rownames(mat)
  out = imap(cts, function(ct, ctname) {
    out = mat[,nuclei %in% names(annot)[annot == ct]]
    list(
      n_nucei = ncol(out),
      genemeans = rowSums(as.matrix(out))
    )
  })
  return(out)
}

#### Iterate function through mats saved in fp #### 
fp = "../data/Batuik/snRNA_seq_raw_countmatrices/" # fp from first step split
ll = list.files(fp) %>% setNames(.,gsub(".rds","",.))
agg_counts_ct = imap(ll, function(fn,name) {
  get_ctAgg_mats(readRDS(paste0(fp,fn)))
}, .progress = T)

#### Create pseudobulk matrices per ct with th of number of nuclei #### 
agg_counts_ctMats = imap(cts[], function(ct,name) {
  genes = names(agg_counts_ct[["MB10"]][[ct]][["genemeans"]])
  out = imap_dfr(agg_counts_ct, function(subject_data,name, n_nuc_th = 10) {
    n_nuc = subject_data[[ct]]$n_nucei
    if(length(n_nuc) == 0) {n_nuc = 0}
    genemeans = subject_data[[ct]]$genemeans 
    if(n_nuc >= n_nuc_th) {
      genemeans
    } else {
      rep(NA, times = length(genemeans)) %>% setNames(genes)
    }
  }, .progress = F) 
  out = as.data.frame(out)
  rownames(out) = names(agg_counts_ct)
  out = out[rowSums(is.na(out))<ncol(out),]
  out
}, .progress = T)

#### Normalize per person to the number of total reads multiple by 10^6 (per Batuik) #### 
#### Dont normalize by length wiht snRNA, cols are subjects #### 
agg_counts_ctMats_CPM = imap(agg_counts_ctMats, function(mat,name) {
  subjects = rownames(mat)
  mat = apply(mat,1,function(subject_row) {
    (subject_row/(sum(subject_row))) * 10^6
  }) 
  mat = as.data.frame((mat))
  colnames(mat) = subjects
  mat
}, .progress = T)

per0_th = .2
#### Same preprocessing as usual pre WGCNA #### 
agg_counts_ctMats_CPM_prepro = imap(agg_counts_ctMats_CPM, function(mat,name) {
  subjects = colnames(mat)
  genes = rownames(mat)
  mat_CPM = (mat)
  mat_logCPM  = as.matrix(log2(mat+1))
  mat = mat_logCPM
  median.gene.exprs = matrixStats::rowMedians(mat)
  
  cat("\n",name, median(median.gene.exprs), 
      mean(median.gene.exprs), min(median.gene.exprs), max(median.gene.exprs),"\n")
  
  print(quantile(median.gene.exprs, probs = seq(0, 1, 0.1)))
  
  names(median.gene.exprs) = rownames(mat)
  gene_filter0 = as.vector(rowSums(mat ==0) <= ncol(mat)*per0_th)                #Genes to keep (not more than 20% samples are zero)
  gene_filter1 = as.vector(median.gene.exprs   >= .1)                              #Genes to keep (median exp > 0.1)

  median.gene.exprs_filtered = median.gene.exprs[gene_filter0 & gene_filter1]
  gene_filter2 = as.vector(abs(scale(median.gene.exprs_filtered)) <  3)            #Genes to keep (z-scoreed_median_exp < 3)
  genes_to_keep  = names(median.gene.exprs_filtered)[gene_filter2]

  mat_logCPM = mat_logCPM[rownames(mat_logCPM) %in% genes_to_keep, ]
  mat_logCPM = mat_logCPM[!grepl("^MT-",rownames(mat_logCPM)), ]
  mat_logCPM = as.data.frame(mat_logCPM)
  
  mat_CPM = mat_CPM[rownames(mat_CPM) %in% genes_to_keep, ]
  mat_CPM = mat_CPM[!grepl("^MT-",rownames(mat_CPM)), ]
  mat_CPM = as.data.frame(mat_CPM)
  
  # mat_logCPM  = log2(mat+1)
  list(CPM = mat_CPM,
       logCPM = mat_logCPM)
}, .progress = F)

n_genes = map_int(agg_counts_ctMats_CPM_prepro, function(mat) {
  nrow(mat[["logCPM"]]) 
})
jpeg(paste0("plots/preprocess/n_genes_Batuik","_per0th_", as.character(per0_th),".jpeg"), height = 1000, width = 1000)
barplot(n_genes ,main="Number of genes after prepro",ylab="Freqency",cex.names = .5,las=2)
dev.off()



#### Same preprocessing as usual pre WGCNA #### 
#### remove MB8  (technical replicate of MB8-2)#### 	
library(sva)
library(jaffelab)
set.seed(123)
library(limma)

facts = names(coldata)
scree_fp = "plots/scree_PBct_batuik/"
dir.create(scree_fp)
PCheat_fp = "plots/PCheat_PBct_batuik/"
dir.create(PCheat_fp)

mat.name = "Glia"
mat_list = agg_counts_ctMats_CPM_prepro[["VIP_CRH"]]
lets = paste(paste(letters, collapse ="|"),paste(LETTERS, collapse ="|"), sep ="|")
get.cleaned.obj = function(mat_list, mat.name, protect){
  mat = mat_list$logCPM
  mat = mat[,colnames(mat) != "MB8"]
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
    if(sum(grepl(lets,c)) == 0 && !grepl("_date",name)) {
      # cat(name)
      c = gsub("\\,","\\.",c)
      c = scale(as.numeric(c))
    } else {
      c = as.factor(c)
    }
  })
  
  check = cbind.data.frame(
    PCimportance          = t(summary(pcobj)$importance["Proportion of Variance",,drop=F])*100 ,
    PC_Age_cor            = round(cor(pcobj$x,(coldata_cr_scale.factor$Age))*100,1) , 
    PC_PMI_hrs_cor            = round(cor(pcobj$x,(coldata_cr_scale.factor$PMI_hrs))*100,1)   ,
    PC_MeanRead_hrs_cor            = round(cor(pcobj$x,(coldata_cr_scale.factor$Mean_reads_per_nucleus_cellranger))*100,1) ,                     
    PC_Dx_cor             = round(cor(pcobj$x,as.numeric(factor(coldata_cr_scale.factor$Diagnosis)))*100,1)               )
  
  
  jpeg(paste0(PCheat_fp,mat.name,".jpeg"))
  pheatmap::pheatmap(check,cluster_rows = F, cluster_cols = F, cellwidth = 20)
  dev.off()
  
  mod = model.matrix(~ Age + PMI_hrs + Mean_reads_per_nucleus_cellranger +
                     Diagnosis +
                     Gender + PC1,
                     data = coldata_cr_scale.factor)

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

agg_counts_ctMats_CPM_prepro_sva = imap(agg_counts_ctMats_CPM_prepro, function(mat_list,name) {
  get.cleaned.obj(mat_list = mat_list, mat.name = name, protect = 1)
}, .progress = T)
saveRDS(agg_counts_ctMats_CPM_prepro_sva,file =  "results/agg_mats_Batuik/Batuik_aggctMats_CPM_prepro_sva.rds")
