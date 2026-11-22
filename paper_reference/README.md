# Reference values from the paper

Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557, https://doi.org/10.1038/s41380-025-03040-x

compare_to_paper.py and plot_paper_figures.py check our results against these files.

| File | From | Contents |
|---|---|---|
| paper_summary.tsv | Fig. 2B | Number of significant genes, down and up |
| paper_degs.tsv | Fig. 2D | Top 30 down and top 30 up genes, log2 fold change and adjusted p-value |
| paper_volcano_labels.tsv | Fig. 2A | The genes the paper labels on its volcano plot |
| paper_reactome_ora.tsv | Fig. 2E | Top 10 Reactome pathways |
| paper_gsea_go_mf.tsv | Fig. 2F | Top 10 GO Molecular Function terms in each direction |
| paper_fig1b_qpcr.tsv | Fig. 1B | AKAP11 exon 7 RT-qPCR, relative to WT (= 1). Read off the bar chart, so approximate (± 0.02) |

If you get the paper's Supplementary Table 1 (all significant genes), you can replace paper_degs.tsv with it, keeping the same three columns.
