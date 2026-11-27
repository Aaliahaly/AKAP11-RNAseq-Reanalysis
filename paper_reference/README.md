# Reference values from the paper

Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557, https://doi.org/10.1038/s41380-025-03040-x

compare_to_paper.py and plot_paper_figures.py check our results against these files.

**Why these values are in files and not in the code:** they are facts copied from the paper, not choices made by the analysis. Keeping them here means the Python code contains no paper-specific numbers, sample names or gene names, so the same code can be pointed at another paper by replacing these files. The file names themselves are set in `settings.py` (PART 9).

All files are tab-separated, with a header row.

| File | From | Contents |
|---|---|---|
| paper_info.tsv | whole paper | Citation, the gene measured in Fig. 1B and how, the sample that looks most different in the paper's heatmap (Fig. 2C), and the paper's headline GSEA finding (its keywords and direction) |
| paper_summary.tsv | Fig. 2B and others | Number of significant genes (total, down, up), the paper's two cutoffs, and two values read off figures: the longest bar in Fig. 2E (−log10 FDR) and the top of the Fig. 2F color scale (−log10 padj) |
| paper_methods.tsv | Methods | Each step of the paper's pipeline next to ours. Text in `{curly brackets}` in the `this_project` column is filled in by the code from `settings.py` and the results, so this table can't go out of date. `affects_gene_list` says whether the step comes before the list of significant genes is made (`yes`) or after it (`no`, e.g. GSEA); only `yes` steps are given as reasons for differences in that list |
| paper_degs.tsv | Fig. 2D | Top 30 down and top 30 up genes, log2 fold change and adjusted p-value |
| paper_volcano_labels.tsv | Fig. 2A | The genes the paper labels on its volcano plot |
| paper_reactome_ora.tsv | Fig. 2E | Top 10 Reactome pathways |
| paper_gsea_go_mf.tsv | Fig. 2F | Top 10 GO Molecular Function terms in each direction. The direction labels must match `higher_label` / `lower_label` in `settings.py` |
| paper_fig1b_qpcr.tsv | Fig. 1B | AKAP11 exon 7 RT-qPCR per clone, relative to WT (= 1), in the paper's bar order, with the paper's bar labels. Read off the bar chart, so approximate (± 0.02) |
| paper_fig2c_columns.tsv | Fig. 2C | The paper's heatmap column order, with a short label for each sample |
| paper_not_reproducible.tsv | Figs. 1C, 1D, 3 | Panels that can't be made from RNA-seq data, and why |

Values read off a figure (rather than from a table in the paper) are approximate; the report says so where they are used.

The code counts the rows of these files instead of assuming their size. For example, our Fig. 2D shows as many top genes per direction as paper_degs.tsv has (30), and Fig. 2F shows as many terms per direction as paper_gsea_go_mf.tsv has (10).
