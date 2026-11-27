# Reproducing the paper's figures

Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557, https://doi.org/10.1038/s41380-025-03040-x

Every panel of the paper's figures, what this project made for it, and why any differences are expected.

| Panel | What it shows | Our figure | Status |
|---|---|---|---|
| Fig. 1B | AKAP11 per clone (RT-qPCR on exon 7) | `12_fig1B_target_gene_relative.png` | Reproduced with differences |
| Fig. 1C | Protein (Western blot) | — | Not reproducible from RNA-seq |
| Fig. 1D | Neuron images (immunostaining) | — | Not reproducible from RNA-seq |
| Fig. 2A | Volcano plot | `12_fig2A_volcano.png` | Reproduced with differences |
| Fig. 2B | Number of significant genes | `12_fig2B_deg_counts.png` | Reproduced with differences |
| Fig. 2C | Clustered heatmap of significant genes | `12_fig2C_heatmap.png` | Reproduced |
| Fig. 2D | Top down and up genes | `12_fig2D_top_genes.png` | Reproduced with differences |
| Fig. 2E | Reactome ORA | `10_reactome_ora.png` | Reproduced with differences |
| Fig. 2F | GSEA, GO Molecular Function | `12_fig2F_gsea_go_mf.png` | Reproduced with differences |
| Fig. 3A–E | DNA methylation (WGBS) | — | Not reproducible from RNA-seq |

## Fig. 1B: AKAP11 per clone (RT-qPCR on exon 7)

- Our RNA-seq levels relative to WT: KO clone 21 = 0.55, KO clone 20 = 0.77, KO clone 16 = 0.65. The paper (RT-qPCR on exon 7): KO clone 21 = 0.43, KO clone 20 = 0.64, KO clone 16 = 0.44.
- Matches: AKAP11 is lower than in WT in every knockout clone in both, and KO clone 20 keeps the most AKAP11 in both.
- Differs: the drop is smaller in RNA-seq (average 0.65 vs 0.50). The paper measured RT-qPCR on exon 7, the region targeted by the knockout edit, so it sees only intact copies. RNA-seq counts reads along the whole gene, including the edited copy, which is still transcribed. A heterozygous knockout also keeps one working copy, so about 0.5 is the expected level for the intact region.
- To check the targeted region directly, you could count reads on it in the STAR BAM files from Part 1.

## Fig. 1C: Protein (Western blot)

- Not reproducible from RNA-seq: it measures protein on a gel. The closest RNA-seq evidence is Fig. 1B (the knocked-out gene's RNA in every clone).

## Fig. 1D: Neuron images (immunostaining)

- Not reproducible from RNA-seq: these are microscope images. The paper's cell-type deconvolution of the RNA-seq data, which found neurons predominant in both groups, is the sequencing-based counterpart.

## Fig. 2A: Volcano plot

- Matches: the same layout, cutoffs and colors. Of the 19 genes the paper labels, 16 can be plotted here and 12 of those are significant; 16 of the 16 that are in the paper's gene table are on the same side as in the paper.
- Differs: our y-axis reaches 72 (top gene CYP26C1), the paper's about 7 (top gene NDN). Our model gives much smaller p-values for the strongest genes, so the genes the paper labels sit lower down in our plot, under genes such as CYP26C1, OTX1, HOXC10.
- Why: the paper's model (limma-voom, robust fit) shrinks extreme statistics; ours (PyDESeq2, a Wald test) does not. The two methods agree on which genes change and in which direction (55 of 55 of the paper's top genes go the same way, r = 0.98), but not on how small the smallest p-values are.
- Not in our results (renamed between Ensembl 104 (paper) and GENCODE v47 (ours), or removed by our low-count filter): C5orf17, ZBED6CL.
- In our results but not plotted, because DESeq2 gave them no adjusted p-value (it leaves out some low-count or outlier genes when correcting p-values): PCDHGC4.

## Fig. 2B: Number of significant genes

- Matches: 848 significant genes here, 849 in the paper, with the same cutoffs (adjusted p < 0.05, |log2FC| ≥ 1).
- Differs: the split is 475 down / 373 up here and 421 / 428 in the paper.
- Why: genes close to the cutoffs can tip either way. These steps differ from the paper's: trimming (Trim Galore (nf-core/rnaseq) vs Trimmomatic); counting (Salmon + tximport vs HTSeq); annotation (GENCODE v47 vs Ensembl 104); gene filter (>= 10 counts in >= 3 samples (26,288 genes) vs mean count >= 20 (23,590 genes)); statistics (PyDESeq2 vs limma-voom, robust fit). Each one moves a few genes across a cutoff.

## Fig. 2C: Clustered heatmap of significant genes

- 848 genes split into two blocks, high in WT and high in the knockouts, as in the paper's heatmap.
- Share of these genes that change in the same direction as the overall result: KO clone 16 100.0%, KO clone 21 99.8%, KO clone 20 96.6%.
- Similarity of the clones' changes: KO clone 16 vs KO clone 21 r = 0.96; KO clone 20 vs the others r = 0.83 and 0.81.
- KO clone 20 follows the shared pattern least closely, and the paper's heatmap shows the same: its column looks different from the other clones. It is also the most deeply sequenced knockout sample (DESeq2 size factor 1.19), which points to a technical difference.
- Differs: the row order and exact colors, because the paper clustered its own gene list, with its own normalization (limma-voom, robust fit) and clustering settings, which the methods don't state.

## Fig. 2D: Top down and up genes

- Overlap with the paper's lists: 7 of our top 30 down genes and 1 of our top 30 up genes are also in the paper's lists (shaded in the figure).
- Differs: the two lists are mostly different genes, even though 55 of the 55 paper genes we found change the same way here (r = 0.98). Both tables rank genes by adjusted p-value, and the two statistical methods order the strongest genes differently: ours puts very strong, consistent changes such as CYP26C1, OTX1, HOXC10 at the top.
- The gene-by-gene comparison of the paper's 60 genes is in 11_paper_gene_comparison.tsv and 11_paper_vs_our_fold_change.png.

## Fig. 2E: Reactome ORA

- Our version is 10_reactome_ora.png, drawn the same way: the top 10 pathways by FDR among the significant genes.
- 7 of our top 10 pathways are also in the paper's top 10. Our top 3: Signaling By GPCR, GPCR Ligand Binding, GPCR Downstream Signaling. 7 of the paper's 10 pathways are significant here. Not significant: Non-integrin membrane-ECM interactions, ECM proteoglycans, Complement cascade.
- Differs: our strongest pathway reaches -log10(FDR) = 8.9, the paper's about 3.2 (read off its figure). Our gene list is a little different, our ORA (hypergeometric test, reactome_human.gmt) tests against all genes we measured, and the paper's (ShinyGO 0.77, curated Reactome) used its own background. These change the p-values more than the order of the pathways.

## Fig. 2F: GSEA, GO Molecular Function

- 9 GO Molecular Function terms are significant here (9 down, 0 up). The largest term tested has 1,068 genes.
- Of the paper's 20 terms, 9 are in our results; 9 of those change in the same direction, and 2 are also significant here.
- Not in our results: Structural Constituent Of Ribosome, Rna Binding, Dna Binding Transcription Activator Activity, Dna Binding Transcription Factor Activity, Regulatory Region Nucleic Acid Binding, Transcription Regulator Activity, Structural Molecule Activity, Oxidoreductase Activity, Phospholipid Binding, Lipid Binding, Signaling Receptor Binding. GO term names and gene sets differ between the GO version the paper used and the GO_Molecular_Function_2025 library we used.
- Differs: the paper's colour scale reaches -log10(padj) of about 30; ours stops at 3. The paper's GSEA (fgsea, GO Molecular Function) estimates very small p-values precisely; GSEApy with 1,000 permutations can't report an FDR much below 1/1,000. The gene ratio and dot size are our closest equivalents to the paper's.

## Fig. 3A–E: DNA methylation (WGBS)

- Not reproducible from this data: Fig. 3 comes from whole-genome bisulfite sequencing, a different experiment with its own FASTQ files and pipeline (for example nf-core/methylseq). Fig. 3E also needs those methylation results, because it uses genes that change in both expression and methylation. It is outside this RNA-seq project, but would be the natural next step.

