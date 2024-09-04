##### Compute association between PGC3 connectivity (for real and null sets) and Zscore of response to CRISPRa activation of PGC3 #####
### First, connectivity to PGC3 vectors are ranked then scaled for each network.
### They are ranked first to remove outlier bias, then scaled so that networks with more genes 
### are not over-represented when taking the median value per gene. Then the median value per gene is taken 
### for each gene across the four DLPFC networks chosen to prioritize consensus genes (LIBD Adult, LIBD Older Adult, GTEx Frontal Cortex, CMC)
### We then computed a linear regression between the median scaled ranked PGC3 conn values and the CRISPa response Zscore.
### This procedure was run for connectivity to PGC3 and for connectivity to 1000 null sets matched to PGC3 (connectivity values were computed in PGC3conn_Fig5.R)

#### Load dependencies ####
library(qs)
library(tidyverse)
library(reshape2)
library(data.table)
library(furrr)
library(sjPlot)
library(readxl)

#### Gene network information ####
setwd(here::here())
gi <- readRDS("../data/wideforms_WGCNA/wf3.6_AmyACC_Hartl_allMAGMA.rds")

#### CRISPRa data from Townsley et al. 2021 ####
CRISPr_screen_meta <- read_excel("../data/Townsley2022/media-1.xlsx",sheet = "META_DEGs")
CRISPr_screen_meta$Symbol = CRISPr_screen_meta$MarkerName
CRISPr_screen_meta = merge(CRISPr_screen_meta,gi[,c("ensemblID","Symbol")])

#### Get scaled ranked connectivity values of connectivity to PGC3 for each network. For real set.####
### Ranked first to remove outlier bias, then scaled so that networks with more genes 
### are not over-represented when taking the median value per gene

nnres_PGC = qread("PGC3_real/nn_residual_wide.qs")  %>% as.data.frame() %>%
  mutate(across(everything(),~ rank(-.x,ties.method = "first", na.last = "keep")))%>%
  scale() %>%
  as.data.frame()
nnres_PGC = t(as.matrix(nnres_PGC)) %>% as.data.frame()
nnres_PGC$Info = rownames(nnres_PGC)
nnres_PGC$Network = gsub("adj_","",stringr::str_split_fixed(rownames(nnres_PGC), n = 2, pattern = "__")[,1])
nnres_PGC$RiskSet = stringr::str_split_fixed(rownames(nnres_PGC), n = 2, pattern = "__")[,2]
nnres_PGC$Iteration = "Original"
nnres_PGC = nnres_PGC[grepl("DLPFC|MIR|Frontal",nnres_PGC$Network),]

nnres_PGC_lg = melt(as.data.table(nnres_PGC), measure.vars = grep("ENSG000", names(nnres_PGC), value = T),
                    value.name = "Rank", variable.name = "Gene")
nnres_PGC_lg$ensemblID = nnres_PGC_lg$Gene
nnres_PGC_lg = merge(CRISPr_screen_meta[c("ensemblID","Zscore")],nnres_PGC_lg, by = "ensemblID")
nnres_PGC_lg = merge(nnres_PGC_lg, gi[,c("ensemblID","Symbol","Length","GC_content")])
nnres_PGC_lg$Network = gsub("_"," ",nnres_PGC_lg$Network)
nnres_PGC_lg$PGC_conn = nnres_PGC_lg$Rank

#### Median value per gene across 4 DLPFC networks
### the four DLPFC networks chosen to prioritize consensus genes (LIBD Adult, LIBD Older Adult, GTEx Frontal Cortex, CMC)
nnres_PGC_lg = nnres_PGC_lg %>% 
  filter(RiskSet == "PGC") %>%
  filter(grepl("DLPFC O|DLPFC Ad|MIR|Frontal",Network)) %>%
  group_by(Gene, RiskSet,Iteration) %>% 
  summarize(n_nets =     sum(!is.na(Rank)),
            Rank_product = prod(Rank),
            Rank_median = median(Rank, na.rm = T),
            Zscore =      median(Zscore, na.rm =T),
            Length =      median(Length, na.rm = T),
            GC_content = median(GC_content, na.rm =T)) %>%
  filter(n_nets == 4)

#### Compute lm
nnres_PGC_lm = nnres_PGC_lg%>%
  group_by(RiskSet,Iteration) %>%
  do(broom::tidy(lm(Zscore~Rank_median + Length + GC_content, .)))


#### Get scaled ranked connectivity values of connectivity to null PGC3 sets for each network.####
### Ranked first to remove outlier bias, then scaled so that networks with more genes 
### are not over-represented when taking the median value per gene

ll_nnres = list.files("PGC3_null", full.names = T) %>% setNames(.,list.files("PGC3_nullv4", full.names = F)) %>%
  grep("residual_wide",., value = T) %>%
  grep("DLPFC Adult|DLPFC Older|MIR|Pergola|DRD2|Frontal",., value = T)%>%
  setNames(.,gsub(" __.*|adj_","",names(.)))

nn_residual = map(ll_nnres, function(fp) {
  qread(fp)
}) 

nn_df  = imap_dfr(nn_residual, function(net_df, net_name) {
  out = net_df[rownames(net_df),]%>% as.data.frame() %>%
    mutate(across(everything(),~ rank(-.x,ties.method = "first", na.last = "keep")))%>%
    scale() 
  out = t(out) %>% as.data.frame()
  out$Network = net_name
  out
}, .progress = T)
nn_df$Info = stringr::str_split_fixed(rownames(nn_df), n = 2, pattern = "__")[,2]
nn_df$RiskSet = stringr::str_split_fixed(nn_df$Info, n = 2, pattern = "___")[,1]
nn_df$Iteration = stringr::str_split_fixed(nn_df$Info, n = 2, pattern = "___")[,2]
nn_df = nn_df[grepl("__PGC__", rownames(nn_df)),]

#### Median value per gene across 4 DLPFC networks. For nulls, per iteration.
options(dplyr.summarise.inform = FALSE)
genes = grep("ENSG000", names(nn_df), value = T) %>% setNames(.,.) ; gene = genes[1]
nn_df_lg = map_dfr(genes, function(gene) {
  nn_dff = nn_df[c("Info","RiskSet","Iteration","Network",gene)]
  nn_dff = data.table::melt(as.data.table(nn_dff), measure.vars = gene,
                  value.name = "Rank", variable.name = "Gene")
  test = nn_dff  %>% 
    group_by(Gene, RiskSet,Iteration) %>% 
    summarize(n_nets =     sum(!is.na(Rank)),
              Rank_product = prod(Rank),
              Rank_median =       median(Rank, na.rm = T)) 
}, .progress = T )
nn_df_lg = nn_df_lg %>%
    filter(n_nets == 4)
nn_df_lg$ensemblID = nn_df_lg$Gene
nn_df_lg = merge(CRISPr_screen_meta[c("ensemblID","Zscore")],nn_df_lg, by = "ensemblID")
nn_df_lg = merge(nn_df_lg, gi[,c("ensemblID","Symbol","Length","GC_content")])

#### Compute lm. For nulls, per iteration.
nn_df_lm = nn_df_lg %>%
  group_by(RiskSet,Iteration) %>%
  do(broom::tidy(lm(Zscore~Rank_median + Length + GC_content, .)))

#### plot histogram
ppp = ggplot(nn_df_lm %>% filter(term == "Rank_median"), 
             aes(x = statistic*-1,  color = 'lightgrey' , fill = 'lightgrey')) + 
  geom_histogram(binwidth=.5, alpha = 1,  color = 'lightgrey' , fill = 'lightgrey') +
  geom_vline(xintercept = nnres_PGC_lm[nnres_PGC_lm$term== "Rank_median",
                                       colnames(nnres_PGC_lm) == "statistic"][[1]]*-1,
             color = "darkmagenta") +
  xlab("t value") +
  theme_bw(base_size = 15)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")
ppp
# save_plot(fig = ppp, filename = "plots/PGCvNullCRISPr_hist.svg",width = 10, height = 14)

