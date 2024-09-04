#### Code to compute ranked connectivity to PGC3 SCZ risk genes, taking into account connectivity to background ####
#### Load dependencies ####
library(tidyverse)
library(limma)
library(MASS)
library(psych)
library(limma)
library(qs)
library(magrittr)
library(Microsoft365R)
library(Matrix)

#### Functions ####
#### Function to compute connectivity to PGC3 and connectivity to all network genes ####
GTEXPGCconktotalreg = function(newdir , dir,adj_dir, genes_list ) {
  setwd(dir)
  dir.create(newdir, showWarnings = FALSE)
  dir.create(paste0(newdir,"/k_total") , showWarnings = FALSE)
  dir.create(paste0(newdir,"/GS.newlists") , showWarnings = FALSE)
  
  ll = list.files(adj_dir)
  print(ll)
  
  walk(ll,~{
    l = .x
    nn = substr(basename(l), 1, nchar(basename(l))-3)
    print(nn)
    if(!file.exists(paste0(nn,"...GS.top.loci.newlists.qs"))| !file.exists(paste0("k_total/",nn,"...ktotal.qs"))){
      
      mat = qread(paste0(adj_dir,l),nthreads = 2)

      colnames(mat) = strsplit2(colnames(mat),"\\.")[,1]
      rownames(mat) = strsplit2(rownames(mat),"\\.")[,1]
      
      diag(mat) = 0
      
      mat = round(mat, 10)
      
      bins = names(genes_list) %>% set_names(.,.)
      bin_genes = imap_dfc(bins,~{
        rr = (colnames(mat) %in% genes_list[[.x]])
      }) %>% set_rownames(rownames(mat))
      
      out = imap_dfc(bins,~{
        rr = data.frame(Matrix::rowSums(mat[,bin_genes[[.x]]])) %>% set_colnames(.y)
      }, .progress = l) %>% set_rownames(rownames(mat))
      qsave(out, paste0(newdir,"/","GS.newlists/", nn,"...GS.newlists.qs"))
      out = data.frame(Matrix::rowSums(mat)) %>% set_rownames(rownames(mat))%>% set_colnames(nn)
      qsave(out, paste0(newdir,"/","k_total/",nn,"...ktotal.qs"))
      
      rm(mat)
      gc()
    }
  },.progress = T)
}
#### Function to regress vector of connectivity to geneset to vector of connectivity to all network genes ####
GTEXPGC_nnresiduals = function(newdir , dir  ) {
  options(future.globals.maxSize= 8912896000)
  setwd(dir)
  setwd(newdir)
  setwd("GS.newlists/")
  lll = list.files(pattern = "...GS.newlists.qs")
  for(i in 1:length(lll)) {
    # i=2
    setwd(dir)
    setwd(newdir)
    setwd("GS.newlists/")
    ll = list.files(pattern = "...GS.newlists.qs")[i]
    names(ll) = ll
    
    net_name = gsub("\\.\\.\\..*","",names(ll))
    
    #Connectivity of each gene to the PGC/GWAS gene loci
    mm_GScon = imap(ll,~ {
      data.frame(qread(.x)) %>% 
        set_names(paste0(strsplit2(.y, "\\.\\.\\.")[,1],"__",colnames(.))) %>% 
        tibble::rownames_to_column("ENSEMBL")
    }) %>% set_names(ll)
    
    setwd(dir)
    setwd(newdir)
    nn_GScon = purrr:::reduce(mm_GScon,full_join, by = "ENSEMBL") %>% tibble::column_to_rownames("ENSEMBL")
    nn_GScon = nn_GScon[sort(rownames(nn_GScon)),]
    names(nn_GScon)=  gsub(" ","_",names(nn_GScon))
    names(nn_GScon)=  gsub("DLPFC_preds2","D_preds2",names(nn_GScon))
    names(nn_GScon)=  gsub("CN__","CN_full__",names(nn_GScon))
    names(nn_GScon)=  gsub("DG__","DG_full__",names(nn_GScon))
    names(nn_GScon)=  gsub("DLPFC__","DLPFC_full__",names(nn_GScon))
    names(nn_GScon)=  gsub("HP__","HP_full__",names(nn_GScon))
    names(nn_GScon)=  gsub("sACC__","sACC_full__",names(nn_GScon))
    names(nn_GScon)=  gsub("AMY__","AMY_full__",names(nn_GScon))
    
    dir.create("nn_GScon")
    qsave(nn_GScon,paste0("nn_GScon/", net_name,"__nn_GScon.qs"))
    
    
    setwd("k_total/")
    ll = list.files(pattern = "...ktotal.qs")[i]
    names(ll) = ll
    
    #Ktotal, Connectivity to all network genes
    mm_ktotal = imap(ll,~ {data.frame(qread(.x)) %>% tibble::rownames_to_column("ENSEMBL")})
    setwd(dir)
    setwd(newdir)
    nn_ktotal = purrr:::reduce(mm_ktotal,full_join, by = "ENSEMBL") %>% tibble::column_to_rownames("ENSEMBL")
    names(nn_ktotal)=  gsub(".*_","",strsplit2(names(ll), "\\.\\.\\.")[,1],"__",names(ll))
    
    ## Compute residuals
    xx = as.matrix(nn_GScon)
    yy = as.numeric(nn_ktotal[, 1]); names(yy) = rownames(nn_ktotal)
    yy = yy[match(rownames(xx),names(yy))]
    stopifnot(all.equal(rownames(xx), names(yy)))
    res = round(lm(xx ~ yy)$residual,8) %>% .[order(rownames(.)),]
    nn_residual_wide = res 
    nn_residual_wide = nn_residual_wide[ order(row.names(nn_residual_wide)), ]
    qsave(nn_residual_wide,paste(net_name,"__nn_residual_wide.qs"))
    
    ## Compute rank residual
    genes_rank = nn_residual_wide %>% as.data.frame() %>%
      mutate(across(everything(),~ rank(-.x,ties.method = "first", na.last = "keep")))
    qsave(genes_rank, paste(net_name,"__genes_rank.qs", nthreads = 2))
    
    setwd(dir)
  }
}

#### Run #### 
library(tidyverse); library(readxl); library(biomaRt)

#### Compute connectivity to PGC3 genelists and connectivity to all network genes ####
dir = getwd()
genes_list = readRDS("data/PGC3_genelistGSconn.rds") ### list of PGC3 genesets, prioritized genes (PGC), genes mapping to 200kbp windows around PGC3 sig SNPs (kbp_200), and genes mapping to 200kbp windows around  to least significant PGC3 SNPs (kbp_200_negative)
newdir = "PGC3_real"
setwd(dir)
adj_dir = "../dir/to/adj_matrices" ##Adjacency matrices available on request
GTEXPGCconktotalreg(newdir , dir,adj_dir, genes_list,is.SYM  )

#### Compute connectivity to PGC3 null sets matched to PGC3 prioritizedgenelist ####
dir = getwd()
genes_list = readRDS("data/PGC3_genelistGSconn.rds") ### list of PGC3 genesets, PGC3 GWAS prioritized genes (PGC), genes mapping to 200kbp windows around PGC3 GWAS sig SNPs (kbp_200), and genes mapping to 200kbp windows around  to least significant PGC3 SNPs (kbp_200_negative)
DLPFC_permutation.output <- readRDS("../data/DLPFC_permutation.output.rds") 
names(DLPFC_permutation.output) = c( "PGC"   ,  "kbp_0"  , "kbp_20"  ,"kbp_50" , "kbp_100",
                                     "kbp_150" ,"kbp_200" ,"kbp_250", "kbp_500")
gs_PGC = genes_list[[c("PGC")]] ### PGC3 GWAS prioritized genes (PGC)
gs_kbp200 = genes_list[[c("kbp_200")]] ### Genes mapping to 200kbp windows around PGC3 GWAS sig SNPs (kbp_200)
gs_kbp500 = readRDS("data/PGC3_kbp500.rds") ### Genes mapping to 500kbp windows around PGC3 GWAS sig SNPs (kbp_500)
PGC_supp_3 <- read_excel("data/PGC3_Supp/Supplementary Table 3.xls")
gs_PGC_supp_3 = map(PGC_supp_3$genes_all, function(st) {str_split_1(string = st ,pattern = ",")}) %>% unlist() %>% unique() ### Genes mapping to PGC3 extended GWAS significant loci

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
biomaRt::listFilters(mart)
mapping_gene <- getBM(filters= "external_gene_name", attributes= c("ensembl_gene_id",
                                                                "external_gene_name"),
                      values=gs_PGC_supp_3,mart= mart)
gsEns_PGC_supp_3 = mapping_gene$ensembl_gene_id
DLPFC_nulls_og = DLPFC_permutation.output[c("PGC")] ### Null sets matched to PGC3 prioritized genes based on median expression in the LIBD DLPFC, GC content, and gene length
og_sizes = c(PGC= length(DLPFC_nulls_og$PGC[[1]]))

#### Remove GWAS associated genes from null sets ####
DLPFC_nulls = map_depth(DLPFC_nulls_og,.depth = 2, function(gs) {
  gs[!gs %in% gs_PGC]
  gs[!gs %in% gsEns_PGC_supp_3]
  gs[!gs %in% gs_kbp200]
  gs[!gs %in% gs_kbp500]
}, .progress = T)

#### Run through permutations, compute connectivities ####
n_perm = 1000
DLPFC_nulls_new = imap(DLPFC_nulls, function(gs_list, name) {
  uni = unique(unlist(gs_list))
  map(1:n_perm, function(i) {
    set.seed(i)
    sample(x = uni, size = og_sizes[name])
  })
}, .progress = T)
PGC_null = DLPFC_nulls_new[["PGC"]][1:n_perm]
names(PGC_null) = paste0("PGC___",as.character(1:n_perm))
genes_list = c(PGC_null)
newdir = "PGC3_null"
setwd(dir)
adj_dir = "../dir/to/adj_matrices" ##Adjacency matrices available on request
GTEXPGCconktotalreg(newdir , dir,adj_dir, genes_list,is.SYM  )


#### Regress vector of connectivity to geneset to vector of connectivity to all network genes ####
newdir = "PGC3_real"
setwd(dir)
GTEXPGC_nnresiduals(newdir , dir )

dir = "~/Library/CloudStorage/OneDrive-JohnsHopkins/R Projects/Omni_conntoGSs"
newdir = "PGC3_null"
setwd(dir)
GTEXPGC_nnresiduals(newdir , dir )
