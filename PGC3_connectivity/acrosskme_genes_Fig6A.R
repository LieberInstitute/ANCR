#### Code for Figure 6A, per-gene connectivity to SCZ risk based on across-module association between kME and median magma ####

#### Load dependencies ####
library(qs)
library(limma)
library(furrr)
library(corrplot)
library(matrixStats)
library(MASS)
library(tidyverse)
library("org.Hs.eg.db")
ENT2SYM = function(IDs) {
  return(mapIds(org.Hs.eg.db,
                keys= as.character(IDs) , 
                column="SYMBOL",
                keytype="ENTREZID",
                multiVals="first"))
}
prioritize(dplyr)

##### Read gene_network_info df ####
gi <-readRDS(  "/Users/christopher.boscuk/Library/CloudStorage/OneDrive-LieberInstituteforBrainDevelopment/Neuron Omni/240729 Submission/Data and code/gene_network_info.rds")
gi = gi[!is.na(gi$SCZ.PGC3.Adult_brain.ZSTAT), ]

## Specifically evaluate published adult DLPFC networks
pub_DLPFC = c(
  "Gandal2018"    ,
  "Gandal2018PE" ,
  "Gandal2018PE_cs" ,
  "Pergola2017" ,
  "Pergola2019"  ,
  "Pergola2020"  ,
  "Radulescu2020"
) %>% paste0(., ".modules")
rem_nets = grep(measure, names(gi), value = T) %>% .[!. %in% pub_DLPFC]
gi = gi %>% dplyr::select(!rem_nets)

#### Extract modules list ####
measure = ".modules"
gm_AB = map(gi %>% dplyr::select(pub_DLPFC) %>% setNames(., gsub(".modules", "", names(.))), ~
              {
                genes = gi$ensemblID
                mod_v = ..1
                mod_names = unique(na.omit(mod_v)) %>% setNames(., .)
                # print(mod_names)
                out = imap(mod_names, ~ {
                  genes[which(mod_v == ..1)]
                })
              })

#### Compute mean magma per module ####
st = "SCZ.PGC3.Adult.brain.ZSTAT"
magmaModStats = imap(gi %>% dplyr::select(dplyr::contains(measure)) %>% setNames(., gsub(measure, "", names(.))),
            ~ {
              print(paste0(..2, ".KME"))
              out = cbind(gi %>% dplyr::select(dplyr::contains(paste0(..2, ".KME")))) %>% t %>% as.data.frame
              colnames(out) = gi$ensemblID
              out = out[, colSums(is.na(out)) == 0]
              out = cbind(module = gsub(paste0(..2, ".KME"), "", rownames(out)), out)
              MAGMA_info = imap_dfr(gm_AB[[..2]], ~ {
                genes_iter = ..1
                v_MAGMA = gi[gi$ensemblID %in% genes_iter, grepl(c(st), colnames(gi))]
                out1 = c(
                  module = ..2,
                  size = length(genes_iter),
                  mean = mean(v_MAGMA),
                  median = median(v_MAGMA)
                )
                names(out1) = c(names(out1)[1:2], paste0(names(out1)[-c(1, 2)], st))
                out1
              })
              out = merge(MAGMA_info, out)
            })
magmaModStats = magmaModStats[map_lgl(magmaModStats, ~ {
  nrow(..1) > 0
})]

#### Compute per-gene ranked across-module statistic of kME vs mean_MAGMA ####
rlm_magmavskme = imap(magmaModStats, ~ {
  print(..2)
  net = ..1
  net_name = ..2
  genes = grep("ENSG", names(net), value = T) %>% setNames(., .)
  net = net[!net$module == "grey", ]
  pb <- progress_estimated(length(genes))
  magmaModStats_rlm = imap_dfr(genes[],  ~ {
    pb$tick()$print()
    mean_MAGMA = as.numeric(net[[3]])
    gene_kme = net[[..1]]
    mod_size = as.numeric(net[["size"]])
    out = broom::tidy(rlm(mean_MAGMA ~ gene_kme +  mod_size))   %>% filter(term == "gene_kme")
    c(ensemblID = ..2, module = gi[gi$ensemblID %in% ..2, c(which(names(gi) == paste0(net_name, ".modules")))], out)
  })
  gene_info = gi[gi$ensemblID %in% magmaModStats_rlm$ensemblID, c(1, 8:40)]
  out = merge(magmaModStats_rlm, gene_info)
  out = cbind(rank = 1:nrow(out), out[order(out$statistic, decreasing = T), ])
})
# saveRDS(rlm_magmavskme,file = "rlm_magmavskme.rds")
# rlm_magmavskme = readRDS("rlm_magmavskme.rds")
rlm_magmavskme_lf = imap_dfr(rlm_magmavskme, ~ {
  out = ..1
  out$module = as.character(out$module)
  cbind(Network = ..2, out)
})

#### Take the sd of median value of the ranked connectivity to high risk modules ####
## remove genes in fewer than 4 networks
rlm_magmavskme_lf$Region = rlm_magmavskme_lf$Network
sum_genes = rlm_magmavskme_lf %>%
  group_by(Symbol) %>%
  summarise(
    mean = mean(rank),
    median = median(rank),
    sd = sd(rank),
    n = length(Symbol)
  ) %>%
  filter(n > 3) %>% arrange((median), .by_group = TRUE) %>% cbind(rank = 1:nrow(.), .)


#### Plot figure ####
## Cerebrum wide SCZ-connected genes
con_genes = c(
  "GRIN1",
  "BSN",
  "SCN8A",
  "SLC12A5",
  "UNC13A",
  "PIP5K1C",
  "UNC5A",
  "GRIN2B",
  "KCNB1",
  "DNM1",
  "WNK2",
  "ADGRL1",
  "CACNA1E",
  "GABRB3"
)
SCHEMA = readRDS(file = "data/SCHEMA_genes_list.rds") # SCHEMA genes from Singh et al., 2022
WGCNAAging_conn = readRDS(file = "data/WGCNAcon_genes_list.rds") # Consensus molecular genes from Pergola et al., 2023 Sci Adv
PGC3_genes = readRDS("data/PGC3_genelistGSconn.rds")
PGC3_genes_ind = which(sum_genes$Symbol %in% ENS2SYM(PGC3_genes$PGC))
SCHEMA_10_ind = which(sum_genes$Symbol %in% SCHEMA$SCHEMA_genes_10)
WGCNAconn_10_ind = which(sum_genes$Symbol %in% WGCNAAging_conn$WGCNAcon_genes_adult)
con_genes_ind = which(sum_genes$Symbol %in% con_genes)

library(RColorBrewer)
cols = RColorBrewer::brewer.pal(name = "Dark2", n = 6)
col_point_size = 1.5
plot(
  sum_genes$median,
  sum_genes$sd,
  xlab = "Median rank of connecitiy to highly risk-enriched modules",
  ylab = "SD rank of connecitiy to highly risk-enriched modules",
  xaxt = 'n',
  yaxt = 'n',
  pch = 16,
  cex = .9,
  col = alpha("lightgrey", .3)
)
axis(side = 1, at = seq(0, 20000, by = 2000))
axis(side = 2, at = seq(0, 10000, by = 2000))
points(
  sum_genes$median[PGC3_genes_ind],
  sum_genes$sd[PGC3_genes_ind],
  col = cols[2],
  pch = 15,
  cex = col_point_size
)
points(
  sum_genes$median[WGCNAconn_10_ind],
  sum_genes$sd[WGCNAconn_10_ind],
  col = cols[4],
  pch = 18,
  cex = col_point_size + .2
)
points(
  sum_genes$median[SCHEMA_10_ind],
  sum_genes$sd[SCHEMA_10_ind],
  col = cols[3],
  pch = 17,
  cex = col_point_size
)
points(
  sum_genes$median[con_genes_ind],
  sum_genes$sd[con_genes_ind],
  col = cols[1],
  pch = 16,
  cex = col_point_size
)

