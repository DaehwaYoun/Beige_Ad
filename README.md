# Cross-talks between Metabolic and Translational Controls during Beige Adipocyte Differentiation
This repository provides codes to review or reproduce results of "Cross-talks between Metabolic and Translational Controls during Beige Adipocyte Differentiation".

Updated : 2024. 10. 11.


## Dependencies
Analyses in this repository were performed with following packages.
* python = 3.9.16
* pandas = 1.5.2
* numpy = 1.23.4
* scipy = 1.9.3
* seaborn = 0.13.2
* matplotlib = 3.8.4
* matplotlib-venn = 0.11.7

Analyses with edgeR were performed with following packages
* R = 4.3.2
* edgeR = 4.0.2

Analyses of single-cell transcriptome were performed with following packages
* R = 4.2.0
* Seurat = 4.3.0.1


## Code location for figures

In the below table, 
* S means Extended Data Figure.
* L means the left subpanel.
* R means the rght subpanel.
* T means the upper (top) subpanel.
* B means the lower (bottom) subpanel.

Figure|Path
---|---
1A|(Experimental result & Manually drawn)
1B|/Results/Plots/Metagene/Abundant/v20240128/adi_otherP_comp_abndt-metagene_heatmap.ipynb
1C_L|/Results/Plots/RPF-PT-corr_better-than-RNA/scatter-cumul/v20240416/RD-up-down-genes_volcano-cumul.ipynb
1C_R|/Results/Plots/RPF-PT-corr_better-than-RNA/scatter-cumul/v20240416/RD-up-down-genes_volcano-cumul.ipynb
1D|/Results/Plots/Mitochondria/GSEA/GSEA_NES_scatter/v20240225/GSEA_NES_scatter.ipynb
1E|/Results/Plots/Mitochondria/PPI/Hierarchical_edge_bundling/RD_GSEA_GOCC_lowest3terms_lead/v20240309/circleHeat.ipynb
1F|/Results/Plots/Mitochondria/violin_box_kde/v20240425/OXPHOS-others_logRD-logFCz_box.ipynb
1G|/Results/Plots/Mitochondria/violin_box_kde/v20240425/OXPHOS-others_logRD-logFCz_box.ipynb
S1A_T|(Experimental result)
S1A_B|/Results/Plots/Marker_Exp/v20240426/adi_Marker_Exp_barplot.ipynb
S1B|(Experimental result)
S1C|/Results/Plots/3ntP/v20240501/3ntP-swarm_CDS-UTR-bar.ipynb
S1D|/Results/Plots/PCA/RNA-RPF-PT/v20240426/adi_RNA-RPF-PT_PCA_plot.ipynb
S1E|/Results/Plots/RPF-PT-corr_better-than-RNA/correlation_heatmap/v20240426/logFC_correlation_heatmap.ipynb
S1F|/Results/Plots/RPF-PT-corr_better-than-RNA/scatter-cumul/v20240220/RD-up-down-genes_scatter.ipynb
S1G|/Results/Plots/RPF-PT-corr_better-than-RNA/scatter/v20240918/RD-DEG_density-scatter_lowess-reg.ipynb
S1H|/Results/Plots/Mitochondria/GSEA/GSEA_NES_scatter/v20240224/GSEA_NES_scatter.ipynb
S1I|/Results/Plots/Mitochondria/GSEA/GSEA_NES_scatter/v20240502/GSEA_NES_scatter.ipynb
S1J|/Results/Plots/Mitochondria/violin_box_kde/v20240509/mito_RNAlogNc_box.ipynb
2A|/Results/Plots/Mitochondria/OMM_proxy/v20240502/APEX-seq_RDlogFC_scatter.ipynb
2B|/Results/Plots/Mitochondria/EnrichR/v20240316/APEXseq_G3_EnrichR_dot_plots.ipynb
2C|/Results/Plots/Mitochondria/OMM_proxy/v20240502/APEX-seq_RDlogFC_scatter.ipynb
2D|/Results/Plots/Mitochondria/violin_box_kde/OXPHOS_complex/v20240811/OXPHOS_RDlogFC_strip-point.ipynb
2E|/Results/Plots/Mitochondria/violin_box_kde/OXPHOS_complex/v20240811/OXPHOS_PTlogFC_strip-point.ipynb
2F|/Results/Plots/Mitochondria/bar/v20240417/OXPHOS_PT_Exp-ratio_barplot.ipynb
2G|/Results/Plots/Mitochondria/scatter/v20240811/mito-OXPHOS_RNA-RD_logFC_scatter-strip.ipynb
2H|/Results/Plots/Mitochondria/violin_box_kde/v20240419/mitoTranslation_PTlogFC_box.ipynb
2I|(Manually drawn)
2J_L|(Experimental result)
2J_R|/Results/Experiments/OXPHOS_complex_western_chloramphenicol/v20241005/OXPHOS_band_PT-level_barplot.ipynb
S2A|(Manually drawn)
S2B|/Results/Plots/Mitochondria/OMM_proxy/v20240502/APEX-seq_RDlogFC_scatter.ipynb
S2C|/Results/Plots/Mitochondria/heatmap/v20240503/mtRP_RDlogFC_heatmap.ipynb
S2D|/Results/Plots/Mitochondria/violin_box_kde/OXPHOS_complex/v20240811/OXPHOS_RDlogFC_strip-point.ipynb
S2E|/Results/Plots/Mitochondria/violin_box_kde/v20240527/mitoTranslation_R-P-logFC_box.ipynb
3A_T|/Results/Plots/Mitochondria/heatmap/v20240419/comp_mtDNA_RDlogFC_heatmap.ipynb
3A_B|/Results/Plots/Mitochondria/heatmap/v20240419/comp_TCA_RDlogFC_heatmap.ipynb
3B|/Results/Plots/Mitochondria/violin_box_kde/OXPHOS_complex/v20240419/comp_OXPHOS_R-P-relFC_strip-point.ipynb
3C_L|(Experimental result)
3C_R|/Results/Experiments/OXPHOS_complex_western_white-beige/Received_20240416/v20240417/OXPHOS_band_beige-white_PT-level_barplot.ipynb
3D|/Results/Plots/Mitochondria/bar/v20240419/Martinez_OXPHOS_PT_Exp-ratio_barplot.ipynb
3E|/Results/Plots/Mitochondria/bar/v20240827/Martinez_TCA-mito_PT_Exp-ratio_barplot.ipynb
3F|/Results/Experiments/OXPHOS_complex_activity/v20240509/CII-CI-activity.ipynb
3G|/Results/Experiments/OXPHOS_complex_activity/v20240509/CII-CI-activity.ipynb
3H|/Results/Experiments/OXPHOS_complex_activity/v20240509/CII-CI-activity.ipynb
3I|(Manually drawn)
S3A_T|/Results/Plots/Mitochondria/heatmap/v20240319/comp_mtDNA_RNAlogFC_heatmap.ipynb
S3A_B|/Results/Plots/Mitochondria/heatmap/v20240319/comp_mtDNA_RPFlogFC_heatmap.ipynb
S3B_T|/Results/Plots/Mitochondria/heatmap/v20240319/comp_TCA_RNAlogFC_heatmap.ipynb
S3B_B|/Results/Plots/Mitochondria/heatmap/v20240319/comp_TCA_RPFlogFC_heatmap.ipynb
S3C|/Results/Plots/Mitochondria/heatmap/v20240518/comp_TCA_PTlogFC_heatmap.ipynb
S3D|(Experimental result)
S3E|/Results/Plots/Mitochondria/violin_box_kde/v20240418/Martinez_mito-PTexp_box.ipynb
S3F|/Results/Plots/Mitochondria/violin_box_kde/v20240827/Martinez-mito_TCA-PTexp_strip-point.ipynb
4A_L|/Results/Plots/Ribosome_stalling/v20240411/adi_stalling-score_cumul.ipynb
4A_R|/Results/Plots/Ribosome_stalling/v20240422_metagene/adi_top-stalling-score_metagene.ipynb
4B|/Results/Plots/Ribosome_stalling/v20240325/adi_top5pct_stalling-score_codon-count-z_heatmap.ipynb
4C|/Results/Plots/Ribosome_stalling/v20240422/adi_top5pct_stalling-score_codon-count_swarm.ipynb
4D|/Results/Plots/Ribosome_stalling/v20240422_SeqLogo/adi_top-stalling-score_SeqLogo.ipynb
4E|/Results/Plots/Ribosome_stalling/v20240423/stalling_change_volcano-pie.ipynb
4F|/Results/Plots/Ribosome_stalling/v20240422_SeqLogo/adi_SigUpStalling_SeqLogo.ipynb
4G|/Results/Plots/Codon_analysis/metagene_around_codon/v20241010_metaAA/adi_metaAA_plot.ipynb
4H|/Results/Plots/Gene_coverage/v20230725/Tef_gene_coverage_plot.ipynb
4I|/Results/Plots/Ribosome_stalling/v20240423/comp_top5pct-stalling-score_AA-count_swarm.ipynb
S4A|/Results/Plots/Codon_analysis/codon_usage_change_barplot/classic/v20240325/AAcodon_usage_dPCT_barplot.ipynb
S4B|/Results/Plots/Codon_analysis/codon_usage_change_barplot/classic/v20240325/codon_usage_dPCT_barplot.ipynb
S4C|/Results/Plots/Codon_analysis/metagene_around_codon/v20241010/adi_metacodon_plot.ipynb
S4D|/Results/Plots/Codon_analysis/metagene_around_codon/v20241010/adi_metacodon_plot.ipynb
S4E|/Results/Plots/Gene_coverage/v20230725/Tef_gene_coverage_plot.ipynb
S4F|/Results/Plots/Ribosome_stalling/v20240327/Martinez-Xie_top5pct-stalling-score_codon-count-z_heatmap.ipynb
5A|/Results/Plots/Glu_starvation/Glu-related_genes/Ours/v20240928/Glu-related-genes_logFC_heatmap.ipynb
5B|/Results/Experiments/Gln-Glu-assay/Received_20230808/v20240413/Glu-Gln-assay_barplot.ipynb
5C|/Results/Experiments/tRNA charging assay/Received_20231128/v20240413/tRNA-charging-assay_barplot.ipynb
5D|(Experimental result)
5E_L|(Experimental result & Manually drawn)
5E_R|/Results/Experiments/Reporter assay/Received_20240807/v20240922/Glu-Ala_Tet-DOX-BPTES_reporter.ipynb
5F|/Results/Plots/Glu_starvation/Glu-related_genes/Reid_2017/v20241009/Reid_Glu-related_logFC_barplot.ipynb
5G|/Results/Plots/Glu_starvation/scRNAseq/scHeatmap/v20240423/Emont_scHeatmap.ipynb
5H|/Results/Experiments/Marker-gene-Exp_adipogenesis-w-drugs/v20240427/MarkerExp-MSO_barplot.ipynb
S5A|/Results/Plots/Glu_starvation/Glu-related_genes/Ours/v20241006/Glu-related-genes_logFC_barplot.ipynb
S5B|/Results/Plots/Glu_starvation/Glu-related_genes/Reid_2017/v20241009/Reid_Glu-related_logFC_barplot2.ipynb
S5C|/Results/Experiments/Glu_Gln_Conc_MSO/v20241009/Glu-Gln-Conc-MSO_barplot.ipynb
S5D|/Results/Experiments/Glul_KD/v20241008/Glul-level-siGlul_barplot.ipynb
S5E|(Experimental result) 
S5F|/Results/Experiments/Glul_KD/v20241008/Glu-Gln-ratio-siGlul_barplot.ipynb
S5G|/Results/Experiments/Glul_KD/v20241008/Markers-level-siGlul_barplot.ipynb
S5H|/Results/Plots/Glu_starvation/scRNAseq/UMAP/v20240320/Emont_UMAP.ipynb
S5I|/Results/Plots/Glu_starvation/scRNAseq/Gls_Glul/v20231029/Gls_Glul_UMAP_violin.ipynb
S5J|/Results/Experiments/Marker-gene-Exp_adipogenesis-w-moreGlu/v20240925/Glu-Gln-ratio-moreGlu_barplot.ipynb
S5K|/Results/Experiments/Marker-gene-Exp_adipogenesis-w-moreGlu/v20240925/MarkerExp-moreGlu_barplot.ipynb
6A|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation2.ipynb
6B|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation2.ipynb
6C|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation2.ipynb
6D|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/delta_mRNA_stability_scatter.ipynb
6E|/Results/Plots/EnrichR/v20240415/Glu-rich-stalled_EnrichR_dotplot.ipynb
6F|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation.ipynb
6G|(Experimental result)
6H|(Experimental result)
6I|/Results/Experiments/Marker-gene-Exp_adipogenesis-w-CytoD/v20240516/Marker-CytoD_barplot.ipynb
S6A|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation2.ipynb
S6B|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation2.ipynb
S6C|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/Effect_of_stalling_on_translation2.ipynb
S6D|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/bias-corrected-stabl_scatters.ipynb
S6E|/Results/Plots/Ribosome_stalling_effect_on_protein/v20240422/miRNA-target_stabl-cumulative.ipynb
