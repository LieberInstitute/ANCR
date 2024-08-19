## Scripts for the manuscript: 
## *Network-wide risk convergence in gene co-expression identifies reproducible genetic hubs of schizophrenia risk*
"add citation when published""

### Data repository
* [figshare](link to be added at pub date): ANCR scores, connectivity to PGC3 ranks, and other data
* Data pertaining to certain scripts may also be found in their respective sub-directories

### Scripts for computing ANCR (Across Network Convergence of Risk, Figures 2,3,4):
* [compute_ANCR_step1_xGene.R] : Compute network-wide risk convergence onto a module (gene risk vs gene z_kme), for all modules(n_genes>=30) in a set of networks.
* [compute_ANCR_step2_xModule.R] : Compute across module association between module median_MAGMA and z_kme, accounting for module size using a robust linear model. This across module association is used as an index of ANCR.

### Pseudobulk and pre-process published snRNA-seq data  (Figure 4):
* [Batiuk_pseudobulk.R] : Pseudobulk and preprocess snRNAseq data from Batiuk et al.
* [Batiuk_data](https://zenodo.org/records/6921620) : snRNA-seq_and_spatial_transcriptomics.zip. Paper :  https://doi.org/10.1126/sciadv.abn8367
* [Ruzicka_pseudobulk.R] : Preprocess pre-pseudobulked snRNAseq data from Ruzicka et al.
* [Ruzicka_data](https://www.synapse.org/#!Synapse:syn25922167) : pseudobulk_mean_logcounts_filtered.rds. Paper :  https://doi.org/10.1101/2022.08.31.22279406
* [snRNA_WGCNA.R] : Compute networks per celltype in a snRNAseq dataset, using the pre-processed expression matrices constructed in the above scripts.

### Scripts for computing and evaluating gene connectivity to PGC3 risk loci genes (Figures 5,6):
* [PGC3conn_Fig5.R] : Compute ranked connectivity to PGC3 SCZ risk genes (and null sets), taking into account connectivity to background.
* [PGC3conn_CRISPRa_null_Fig5.R] : Evaluate association between PGC3 connectivity (for real and null sets) and Zscore of response to CRISPRa activation of PGC3 eGenes (from Townsley et al).
* [Townsley.et.al_CRISPRaPGC3_data](https://www.biorxiv.org/content/10.1101/2022.03.29.486286v2.supplementary-material) : Published data of Zscore of response to CRISPRa activation of PGC3 eGenes : see their Supplemental Data 1.
* [acrosskme_genes_Fig6A.R] : Compute per-gene connectivity to SCZ risk based on across-module association between kME and median magma.

For any data or code inquiries please contact Giulio Pergola: [Giulio.Pergola@libd.org] [https://www.libd.org/team/giulio-pergola-phd/]
