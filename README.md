## Scripts for the manuscript: 
## *Network-wide risk convergence in gene co-expression identifies reproducible genetic hubs of schizophrenia risk*

### Data repository
* [Zenodo](https://doi.org/10.5281/zenodo.13381309): ANCR scores, connectivity to PGC3 ranks, and other data
* Data pertaining to certain scripts may also be found in their respective sub-directories

### Computing ANCR (Compute_ANCR, Figures 2,3,4):
* ANCR_step1_xGene_RemNetInMod.R : Compute network-wide risk convergence onto a module ("z_kme", gene risk vs gene kme), for all modules in a set of networks. (Fig 2B, 3B, 4). 
* ANCR_step1_xGene_RemNetInNet.R : Compute network-wide risk convergence onto a module ("z_kme", gene risk vs gene kme), for all modules in a set of networks. (Fig 3D). 
* ANCR_step2_xModule_RemNetInMod.R : Compute across module association between module median_MAGMA and z_kme, accounting for module size. This across module association is used as an index of ANCR.(Fig 2B, 3B, 4). 
* ANCR_step2_xModule_RemNetInNet.R : Compute across module association between module median_MAGMA and z_kme, accounting for module size. This across module association is used as an index of ANCR. (Fig 3D). 

### Pseudobulk and pre-process published snRNA-seq data  (snRNA_nets, Figure 4):
* Batiuk_pseudobulk.R] : Pseudobulk and preprocess snRNAseq data from Batiuk et al. (Fig 4). 
* [Batiuk_data](https://zenodo.org/records/6921620) : snRNA-seq_and_spatial_transcriptomics.zip. Paper :  https://doi.org/10.1126/sciadv.abn8367
* Ruzicka_pseudobulk.R] : Preprocess pre-pseudobulked snRNAseq data from Ruzicka et al. (Fig 4). 
* [Ruzicka_data](https://www.synapse.org/#!Synapse:syn25922167) : pseudobulk_mean_logcounts_filtered.rds. Paper :  https://doi.org/10.1101/2022.08.31.22279406
* snRNA_WGCNA.R : Compute networks per celltype in a snRNAseq dataset, using the pre-processed expression matrices constructed in the above scripts. (Fig 4). 

### Computing and evaluating gene connectivity to PGC3 risk loci genes (PGC3_connectivity, Figures 5,6):
* PGC3conn_Fig5.R : Compute ranked connectivity to PGC3 SCZ risk genes (and null sets), taking into account connectivity to background. (Fig 5). 
* PGC3conn_CRISPRa_null_Fig5.R : Evaluate association between PGC3 connectivity (for real and null sets) and Zscore of response to CRISPRa activation of PGC3 eGenes (from Townsley et al). (Fig 5).
* [Townsley.et.al_CRISPRaPGC3_data](https://www.biorxiv.org/content/10.1101/2022.03.29.486286v2.supplementary-material) : Published data of Zscore of response to CRISPRa activation of PGC3 eGenes : see their Supplemental Data 1.
* acrosskme_genes_Fig6A.R : Compute per-gene connectivity to SCZ risk based on across-module association between kME and median magma. (Fig 6).

For any data or code inquiries please contact Giulio Pergola: [Giulio.Pergola@libd.org] [https://www.libd.org/team/giulio-pergola-phd/]
