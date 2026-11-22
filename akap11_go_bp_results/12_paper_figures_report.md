# Reproducing the paper's figures

Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557, https://doi.org/10.1038/s41380-025-03040-x

Every panel of the paper's Figures 1 to 3, what this project made for it, and why any differences are expected.

| Panel | What it shows | Our figure | Status |
|---|---|---|---|
| Fig. 1B | AKAP11 exon 7 RT-qPCR per clone | `12_fig1B_akap11_relative.png` | Reproduced with differences |
| Fig. 1C | AKAP11 protein (Western blot) | — | Not reproducible from RNA-seq |
| Fig. 1D | Neuron images (immunostaining) | — | Not reproducible from RNA-seq |
| Fig. 2A | Volcano plot | `12_fig2A_volcano.png` | Reproduced with differences |
| Fig. 2B | Number of significant genes | `12_fig2B_deg_counts.png` | Reproduced with differences |
| Fig. 2C | Clustered heatmap of significant genes | `12_fig2C_heatmap.png` | Reproduced |
| Fig. 2D | Top 30 down and up genes | `12_fig2D_top_genes.png` | Reproduced with differences |
| Fig. 2E | Reactome ORA | `10_reactome_ora.png` | Reproduced with differences |
| Fig. 2F | GSEA, GO Molecular Function | `12_fig2F_gsea_go_mf.png` | Reproduced with differences |
| Fig. 3A–E | DNA methylation (WGBS) | — | Not reproducible from RNA-seq |

## Fig. 1B: AKAP11 exon 7 RT-qPCR per clone

- Our RNA-seq levels relative to WT: clone 21 = 0.55, clone 20 = 0.77, clone 16 = 0.65. The paper's RT-qPCR: 0.43, 0.64, 0.44.
- Matches: AKAP11 is lower in every knockout clone, and clone 20 keeps the most AKAP11 in both.
- Differs: the drop is smaller in RNA-seq. The paper's qPCR measures only exon 7, the exon targeted by the CRISPR edit, so it sees only intact copies. RNA-seq counts reads along the whole gene, including the edited copy, which is still transcribed. A heterozygous knockout also keeps one working copy, so about 0.5 is the expected level for intact exon 7.
- To check exon 7 directly, you could count reads on exon 7 in the STAR BAM files on Expanse.

## Fig. 1C: AKAP11 protein (Western blot)

- Not reproducible from RNA-seq: it measures protein on a gel. The closest RNA-seq evidence is Fig. 1B (AKAP11 RNA is lower in every clone).

## Fig. 1D: Neuron images (immunostaining)

- Not reproducible from RNA-seq: these are microscope images. The paper's cell-type deconvolution of the RNA-seq data, which found neurons predominant in both groups, is the sequencing-based counterpart.

## Fig. 2A: Volcano plot

- Matches: the same layout, cutoffs and colors. Of the 19 genes the paper labels, 16 can be plotted here and 12 of those are significant, all on the same side as in the paper.
- Differs: our y-axis reaches 72 (top gene CYP26C1), the paper's about 7 (top gene NDN). Our model gives much smaller p-values for the strongest genes, so the genes the paper labels sit lower down in our plot, under genes like CYP26C1 and the HOX genes.
- Why: the paper used limma-voom with robust empirical Bayes, which shrinks extreme statistics; we used DESeq2's Wald test, which does not. The two methods agree on which genes change and in which direction (r = 0.98 for the paper's top genes), but not on how small the smallest p-values are.
- Not in our results, because they were renamed between Ensembl 104 (paper) and GENCODE v47 (ours): C5orf17, ZBED6CL.
- In our results but not plotted, because DESeq2 gave them no adjusted p-value (it leaves out some low-count or outlier genes when correcting p-values): PCDHGC4.

## Fig. 2B: Number of significant genes

- Matches: 848 significant genes here, 849 in the paper, with the same cutoffs (adjusted p < 0.05, |log2FC| ≥ 1).
- Differs: the split is 475 down / 373 up here and 421 / 428 in the paper.
- Why: genes close to the cutoffs can tip either way. Counting (Salmon vs HTSeq), annotation (GENCODE v47 vs Ensembl 104), gene filter and model (DESeq2 vs limma-voom) all differ, and each moves a few genes across a cutoff.

## Fig. 2C: Clustered heatmap of significant genes

- Matches: 848 genes split cleanly into two blocks, high in WT and high in the knockouts, as in the paper's heatmap.
- Share of these genes that change in the same direction as the overall result: clone16 100.0%, clone21 99.8%, clone20 96.6%.
- Similarity of the clones' changes: clone16 vs clone21 r = 0.96; clone20 vs the others r = 0.83 and 0.81.
- Clone 20 follows the shared pattern least closely, and has its own block of genes that are up mainly in clone 20. The paper's heatmap shows the same thing: its clone 20 column looks different from clones 16 and 21. This fits clone 20's technical differences (more intronic and duplicate reads) and its higher remaining AKAP11 (Fig. 1B).
- Differs: the row order and exact colors, because the paper clustered its own gene list, with its own normalization (limma-voom) and clustering settings, which the methods don't state.

## Fig. 2D: Top 30 down and up genes

- Overlap with the paper's lists: 7 of our top 30 down genes and 1 of our top 30 up genes are also in the paper's top 30 (shaded in the figure).
- Differs: the two top-30 lists are mostly different genes, even though the paper's top 60 genes all change the same way here (r = 0.98). Both tables rank genes by adjusted p-value, and DESeq2 and limma-voom order the strongest genes differently: ours puts very strong, consistent changes such as CYP26C1, the HOX genes and OTX1 at the top.
- The gene-by-gene comparison of the paper's 60 genes is in 11_paper_gene_comparison.tsv and 11_paper_vs_our_fold_change.png.

## Fig. 2E: Reactome ORA

- Our version is 10_reactome_ora.png, drawn the same way: top 10 pathways by FDR among the significant genes.
- Matches: GPCR pathways lead both lists (GPCR ligand binding, signalling by GPCR, G alpha (q), class A/1 rhodopsin-like receptors), followed by extracellular matrix organization. 7 of the paper's 10 pathways are significant here. Not significant: Non-integrin membrane-ECM interactions, ECM proteoglycans, Complement cascade.
- Differs: our strongest pathway reaches -log10(FDR) = 8.9, the paper's about 3.2. Our gene list is a little different, and we test against all genes we measured with an up-to-date Reactome file, while the paper used ShinyGO 0.77 with its own curated Reactome version and background. These change the p-values but not which pathways come first.

## Fig. 2F: GSEA, GO Molecular Function

- 9 GO Molecular Function terms are significant here (9 down, 0 up). The largest term tested has 1,068 genes.
- Of the paper's 20 terms, 9 are in our results; 9 of those change in the same direction, and 2 are also significant here.
- Not in our results: Structural Constituent Of Ribosome, Rna Binding, Dna Binding Transcription Activator Activity, Dna Binding Transcription Factor Activity, Regulatory Region Nucleic Acid Binding, Transcription Regulator Activity, Structural Molecule Activity, Oxidoreductase Activity, Phospholipid Binding, Lipid Binding, Signaling Receptor Binding. GO term names and gene sets differ between the GO version the paper used (its names look like MSigDB's) and Enrichr's 2025 library.
- Differs: the paper's colour scale reaches -log10(padj) of about 30; ours stops at 3. The paper used fgsea, which estimates very small p-values precisely; GSEApy with 1,000 permutations can't report an FDR much below 1/1,000. The gene ratio and dot size are our closest equivalents to the paper's.

## Fig. 3A–E: DNA methylation (WGBS)

- Not reproducible from this data: Fig. 3 comes from whole-genome bisulfite sequencing, a different experiment with its own FASTQ files and pipeline (for example nf-core/methylseq). Fig. 3E also needs those methylation results, because it uses genes that change in both expression and methylation.
- It is outside this RNA-seq project, but would be the natural next step.

