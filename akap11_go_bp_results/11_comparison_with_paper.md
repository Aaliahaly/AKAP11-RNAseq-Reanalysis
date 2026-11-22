# Comparison with Farhangdoost et al. (2025)

Paper: Molecular Psychiatry 30:4543-4557, https://doi.org/10.1038/s41380-025-03040-x

## Summary

| Check | Result | Verdict |
|---|---|---|
| Number of significant genes | paper 849, ours 848 (0.1% apart) | similar |
| Paper's top 60 genes | 55 of 55 found go the same way; 37 significant in ours; r = 0.98 | agree |
| Reactome ORA pathways | 7 of the paper's 10 top pathways are significant in ours | agree |
| GSEA: ribosome/protein-synthesis terms lower in knockout | 16 significant terms lower, 0 higher | agree |
| GSEA: the paper's 20 GO Molecular Function terms | 9 of 9 tested go the same way | agree |

## 1. Number of significant genes (paper Fig. 2B)

| | Paper | Ours |
|---|---|---|
| Total | 849 | 848 |
| Lower in comparison group | 421 | 475 |
| Higher in comparison group | 428 | 373 |

The totals differ by 0.1%.

## 2. The paper's top 60 genes (paper Fig. 2D)

- Found in our results: 55 of 60
- Same direction of change: 55 of 55 found (100%)
- Also significant in our results, same direction: 37 of 60
- Correlation of log2 fold changes: 0.98
- Not found (renamed in our annotation, or removed by our low-count filter): C5orf17, ZBED6CL, C9orf64, FRG1DP, LINC02302

Gene-by-gene details: `11_paper_gene_comparison.tsv`. Plot: `11_paper_vs_our_fold_change.png`.

## 3. Reactome pathways (paper Fig. 2E)

Of the paper's top 10 pathways, 10 were tested in our ORA and 7 have FDR < 0.05 in ours.

## 4. GSEA (paper Fig. 2F)

The paper's main GSEA finding is that ribosome and protein-synthesis terms go down in the knockout. In our GSEA, 16 such terms (ribosome, translation, rRNA, tRNA) are significantly lower and 0 are significantly higher.

Top three of them in our results:

- Cytoplasmic Translation (GO:0002181): NES -3.41, FDR 0 (estimated)
- Translation (GO:0006412): NES -3.11, FDR 0 (estimated)
- Ribosomal Small Subunit Biogenesis (GO:0042274): NES -2.38, FDR 0 (estimated)

Of the paper's 20 GO Molecular Function terms, 9 were tested in our GSEA and 9 change in the same direction.

## Why our numbers are not identical to the paper's

| Step | Paper | This project |
|---|---|---|
| Trimming | Trimmomatic | Trim Galore (nf-core/rnaseq) |
| Alignment | STAR | STAR |
| Counting | HTSeq, Ensembl 104 | Salmon + tximport, GENCODE v47 |
| Gene filter | mean count >= 20 (23,590 genes) | >= 10 counts in >= 3 samples |
| Statistics | limma-voom, robust fit | PyDESeq2 |
| GSEA | fgsea, GO Molecular Function | GSEApy, setting in settings.py |
| ORA | ShinyGO, curated Reactome | hypergeometric test, Reactome (.gmt) |

