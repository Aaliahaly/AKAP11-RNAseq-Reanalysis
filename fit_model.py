"""
===================================================================
fit_model.py
===================================================================

WHAT THIS FILE DOES
    The main statistics: differential expression. For every gene it
    asks:

        "Is this gene expressed at a different level in KO than in WT,
         more than we'd expect from normal sample-to-sample variation?"

    It uses PyDESeq2, a Python version of DESeq2, one of the most
    widely used tools for this question. It also makes two expression
    tables (normalized counts and VST) that later plots use.

WHERE IT FITS
    Step 2. run_analysis.py calls fit_model() right after
    read_counts(), and passes it the three tables read_counts() made.

INPUTS
    counts            filtered counts from read_counts.py (genes x samples)
    metadata          each sample's group (WT or KO)
    annotation        each gene's name and plot label
    reference_group   "WT"
    comparison_group  "KO"
    output            the results folder

OUTPUT FILES
    01_sample_qc.tsv                 one row per sample: group, QC
                                     numbers, DESeq2 size factor
    02_differential_expression.tsv   one row per gene: fold change,
                                     p-value, adjusted p-value, direction
    02_significant_genes.tsv         only the genes passing both cutoffs
    02_normalized_counts.tsv         normalized counts (genes x samples)
    02_vst_expression.tsv            VST expression (genes x samples)

HOW DESeq2 WORKS, IN SHORT
    1. Size factors: some samples were sequenced more deeply, so every
       gene has more reads in them. DESeq2 works out one scaling number
       per sample (its "size factor") to make samples comparable.
    2. Dispersion: it estimates how much each gene varies between
       replicates of the same group. With only 3 replicates, one gene's
       estimate is shaky, so DESeq2 borrows information from genes with
       similar expression levels.
    3. Test: for each gene it estimates the log2 fold change (KO vs WT)
       and runs a Wald test: is the change big compared with its
       uncertainty? The result is a p-value per gene.
    4. Correction: the p-values are adjusted for testing ~26,000 genes
       at once (Benjamini-Hochberg), giving "padj".

PYTHON IDEAS USED IN THIS FILE
    .T          "transpose": swap a table's rows and columns
    object.method()
                PyDESeq2 creates an object (the model) that we then
                ask to do things, like model.deseq2()
    elif        "else if": checked only when the if above was False
===================================================================
"""

import pandas as pd

# DeseqDataSet holds the counts and fits the model.
# DeseqStats runs the statistical test on a fitted model.
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from settings import fdr_cutoff            # 0.05
from settings import fold_change_cutoff    # 1.0


def fit_model(counts, metadata, annotation, reference_group, comparison_group, output):

    # ===============================================================
    # STEP 1: Fit the model and compare the two groups
    # ===============================================================

    # Our count table has genes in rows and samples in columns:
    #
    #                 WT_rep1  WT_rep2  ...  AKAP11_KO_clone21
    #   ENSG...0003   9120     8433     ...  9775
    #
    # PyDESeq2 wants the opposite: samples in rows, genes in columns.
    # .T ("transpose") swaps rows and columns.
    counts_by_sample = counts.T

    # The model needs each sample's group. Double brackets [[ ]] keep
    # this as a one-column table (PyDESeq2 needs a table, not a column).
    conditions = metadata[["condition"]]

    # Set up the model.
    #   counts   -> the read counts (samples x genes)
    #   metadata -> the group of each sample
    #   design   -> "~condition" is a formula meaning "expression
    #               depends on the condition (WT or KO)". It's the
    #               simplest possible design: one factor, two levels.
    model = DeseqDataSet(
        counts=counts_by_sample,
        metadata=conditions,
        design="~condition",
    )

    # Fit the model. This one line does steps 1 to 3 of
    # "HOW DESeq2 WORKS" above: size factors, dispersions, and the
    # fold change of every gene. It prints progress messages and takes
    # a minute or two.
    model.deseq2()

    # Say which comparison we want. The order is:
    #   [column name, comparison group, reference group]
    # so this means "KO compared with WT". As a result:
    #   positive log2 fold change -> higher in KO
    #   negative log2 fold change -> lower in KO (higher in WT)
    contrast = ["condition", comparison_group, reference_group]

    # Set up the test. alpha=fdr_cutoff tells PyDESeq2 which adjusted
    # p-value cutoff we will use, which it needs for one of its
    # internal steps (independent filtering, see below).
    test = DeseqStats(model, contrast=contrast, alpha=fdr_cutoff)

    # Run the test on every gene and print a short preview in Terminal.
    test.summary()

    # The test results are a table called results_df, one row per gene:
    #   baseMean       -> average normalized count across all samples
    #   log2FoldChange -> the change, KO vs WT, on a log2 scale
    #   lfcSE          -> the uncertainty of that change
    #   stat           -> the Wald statistic: change / uncertainty
    #   pvalue         -> the p-value of the test
    #   padj           -> the p-value adjusted for testing many genes
    #
    # A few genes get padj = NaN (empty). PyDESeq2 does this on purpose
    # for genes with very low counts ("independent filtering") or with
    # an extreme outlier sample, because their test isn't reliable.
    test_results = test.results_df

    # ===============================================================
    # STEP 2: Add gene names and decide which genes are significant
    # ===============================================================

    # Both tables use gene IDs as row labels, so we can join them side
    # by side. how="right" means: keep every gene that is in the test
    # results (the table on the RIGHT, inside the brackets), and add
    # its gene_name and label columns from the annotation table.
    results = annotation.join(test_results, how="right")

    # Now give every gene a "direction" label:
    #   "Higher in comparison" -> significant AND higher in KO
    #   "Lower in comparison"  -> significant AND lower in KO
    #   "Below thresholds"     -> did not pass both cutoffs
    #
    # First take the two columns we need as plain lists.
    all_padj = list(results["padj"])
    all_fold_changes = list(results["log2FoldChange"])

    directions = []

    # Go through the genes one at a time, by position (0, 1, 2, ...).
    for position in range(len(results)):
        padj = all_padj[position]
        fold_change = all_fold_changes[position]

        # Cutoff 1: the adjusted p-value must be 0.05 or less.
        # A missing padj (NaN) never counts as <= 0.05, so genes
        # without a padj automatically fail this cutoff.
        passes_fdr = padj <= fdr_cutoff

        # Cutoff 2: the fold change must be at least 1 in size, in
        # either direction. abs() removes the minus sign
        # (abs(-1.5) = 1.5), so +1.5 and -1.5 both pass.
        passes_fold_change = abs(fold_change) >= fold_change_cutoff

        # A gene is significant only if it passes BOTH cutoffs.
        if passes_fdr and passes_fold_change:
            if fold_change > 0:
                directions.append("Higher in comparison")
            elif fold_change < 0:
                directions.append("Lower in comparison")
            else:
                # A fold change of exactly 0 can't pass cutoff 2, so
                # this line never actually runs. It's here so every
                # possible case gives the gene a label.
                directions.append("Below thresholds")
        else:
            directions.append("Below thresholds")

    # Add the list as a new column. The list was built in the same
    # order as the table's rows, so each label lands on its own gene.
    results["direction"] = directions

    # Sort the genes so the smallest adjusted p-values come first.
    # If two genes have the same padj, the one with the smaller raw
    # p-value goes first. Genes with no padj go to the bottom.
    results = results.sort_values(["padj", "pvalue"])

    # Save the full results table as a tab-separated file.
    # "output" is a Path (a folder); "/" joins the folder and the file
    # name into one path, e.g. akap11_go_bp_results/02_differential...
    results.to_csv(output / "02_differential_expression.tsv", sep="\t")

    # ===============================================================
    # STEP 3: Make two expression tables for the plots
    # ===============================================================

    # The row and column labels we'll need below.
    sample_names = metadata.index
    gene_ids = counts.index

    # ---------------------------------------------------------------
    # Table A: normalized counts
    # ---------------------------------------------------------------

    # Normalized count = raw count / the sample's size factor.
    # Example: clone 20 was sequenced about twice as deeply, so its
    # size factor is larger and its counts are divided by more. After
    # this, a gene's numbers can be compared fairly across samples.
    #
    # PyDESeq2 stores them in model.layers["normed_counts"] as a plain
    # grid of numbers with no labels, so we wrap them in a DataFrame
    # and give it the sample names (rows) and gene IDs (columns).
    normalized = pd.DataFrame(
        model.layers["normed_counts"],
        index=sample_names,
        columns=gene_ids,
    )

    # Swap rows and columns so genes are in rows again.
    normalized = normalized.T

    # Save it. plot_paper_figures.py uses it for the paper's Fig. 1B.
    normalized.to_csv(output / "02_normalized_counts.tsv", sep="\t")

    # ---------------------------------------------------------------
    # Table B: variance-stabilized (VST) expression
    # ---------------------------------------------------------------

    # Raw and normalized counts have a problem for plots: genes with
    # big counts also vary by big amounts, so in a heatmap or a
    # correlation they would dominate everything. The VST
    # ("variance-stabilizing transformation") is a log-like transform
    # that evens this out, so every gene counts about equally.
    #
    # use_design=False means "don't use the WT/KO groups for this
    # step", so the transform can't make the groups look more
    # different than they are.
    model.vst(use_design=False)

    # Same wrapping as for Table A.
    transformed = pd.DataFrame(
        model.layers["vst_counts"],
        index=sample_names,
        columns=gene_ids,
    )

    # Swap rows and columns so genes are in rows again.
    transformed = transformed.T

    # Save it. plot_paper_figures.py uses it for the paper's Fig. 2C.
    transformed.to_csv(output / "02_vst_expression.tsv", sep="\t")

    # ===============================================================
    # STEP 4: Save the sample QC table and the significant genes
    # ===============================================================

    # Add each sample's size factor to the sample table, so we can see
    # how much each sample was scaled. A value near 1 means "average";
    # clone 20's is the largest (about 1.19) because it is the deepest.
    metadata["size_factor"] = model.obs["size_factors"]

    # Save the sample table. index_label="sample" puts the word
    # "sample" at the top of the first column in the file.
    metadata.to_csv(output / "01_sample_qc.tsv", sep="\t", index_label="sample")

    # Keep only the genes that passed both cutoffs: True/False per gene,
    # then keep the True rows.
    is_significant = results["direction"] != "Below thresholds"
    significant = results[is_significant]

    significant.to_csv(output / "02_significant_genes.tsv", sep="\t")

    # Print the number in Terminal, e.g.
    #   Genes passing both thresholds: 848
    # print() with commas prints the pieces separated by spaces.
    number_significant = len(significant)
    print("Genes passing both thresholds:", number_significant)

    # Send three tables back to run_analysis.py:
    #   results     -> test results for every gene in the model
    #   normalized  -> normalized counts (genes x samples)
    #   transformed -> VST expression (genes x samples)
    return results, normalized, transformed
