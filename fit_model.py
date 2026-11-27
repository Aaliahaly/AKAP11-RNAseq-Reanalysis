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
    zip(a, b)   go through two lists side by side, one pair at a time
    def         make our own small function, to avoid writing the
                same code twice
===================================================================
"""

import pandas as pd

# DeseqDataSet holds the counts and fits the model.
# DeseqStats runs the statistical test on a fitted model.
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Settings (see settings.py for why each value was chosen).
from settings import fdr_cutoff              # 0.05
from settings import fold_change_cutoff      # 1.0
from settings import condition_column        # "condition"
from settings import sample_column           # "sample"
from settings import higher_label            # "Higher in comparison"
from settings import lower_label             # "Lower in comparison"
from settings import not_significant_label   # "Below thresholds"

# The names of the five files this step writes.
from settings import sample_qc_file
from settings import differential_expression_file
from settings import significant_genes_file
from settings import normalized_counts_file
from settings import vst_expression_file


def layer_to_table(model, layer_name, sample_names, gene_ids):
    """Turn one of PyDESeq2's number grids into a labeled table.

    PyDESeq2 keeps its results as plain grids of numbers with no
    labels (samples in rows, genes in columns). This adds the labels,
    then swaps rows and columns (.T) so genes are in rows, like the
    rest of the project. We need this twice (STEP 3), so it is a
    function instead of the same code written twice.
    """
    table = pd.DataFrame(model.layers[layer_name], index=sample_names, columns=gene_ids)
    return table.T


def fit_model(counts, metadata, annotation, reference_group, comparison_group, output):

    # ===============================================================
    # STEP 1: Fit the model and compare the two groups
    # ===============================================================

    # Set up the model.
    #   counts   -> the read counts. Ours has genes in rows, but
    #               PyDESeq2 wants samples in rows, so .T swaps them.
    #   metadata -> each sample's group. Double brackets [[ ]] keep it
    #               as a one-column table, which PyDESeq2 needs.
    #   design   -> "~condition" means "expression depends on the
    #               group (WT or KO)": one factor, two levels. We build
    #               it from the column name in settings.py, so it always
    #               matches the sample sheet: "~" + "condition".
    model = DeseqDataSet(
        counts=counts.T,
        metadata=metadata[[condition_column]],
        design="~" + condition_column,
    )

    # Fit the model: size factors, dispersions, and every gene's fold
    # change (steps 1 to 3 of "HOW DESeq2 WORKS" above). It prints
    # progress messages and takes a minute or two.
    model.deseq2()

    # Run the test, comparing KO with WT. The contrast is
    #   [group column, comparison group, reference group]
    # so a positive log2 fold change means higher in KO, and a
    # negative one means lower in KO (higher in WT).
    # alpha=fdr_cutoff tells PyDESeq2 our adjusted p-value cutoff, which
    # it uses to decide which genes are worth testing ("independent
    # filtering": genes with very low counts get padj = NaN, empty).
    test = DeseqStats(
        model,
        contrast=[condition_column, comparison_group, reference_group],
        alpha=fdr_cutoff,
    )

    # Test every gene, and print a short preview in Terminal.
    test.summary()

    # ===============================================================
    # STEP 2: Add gene names and decide which genes are significant
    # ===============================================================

    # test.results_df has one row per gene:
    #   baseMean       -> average normalized count across all samples
    #   log2FoldChange -> the change, KO vs WT, on a log2 scale
    #   lfcSE          -> the uncertainty of that change
    #   stat           -> the Wald statistic: change / uncertainty
    #   pvalue, padj   -> the p-value, and the p-value adjusted for
    #                     testing many genes
    #
    # Both tables use gene IDs as row labels, so .join() puts them side
    # by side. how="right" keeps every gene in the test results (the
    # table in the brackets) and adds its name and label.
    results = annotation.join(test.results_df, how="right")

    # Give every gene one of three labels (set in settings.py):
    #   higher_label          -> significant AND higher in KO
    #   lower_label           -> significant AND lower in KO
    #   not_significant_label -> did not pass both cutoffs
    directions = []

    # zip() walks through the two columns side by side, so each round
    # of the loop gets one gene's padj and its fold change.
    for padj, fold_change in zip(results["padj"], results["log2FoldChange"]):

        # Cutoff 1: padj of fdr_cutoff (0.05) or less. A missing padj
        # (NaN) never counts as <= 0.05, so it fails automatically.
        passes_fdr = padj <= fdr_cutoff

        # Cutoff 2: a change of at least fold_change_cutoff (1.0) in
        # either direction. abs() drops the minus sign, so +1.5 and
        # -1.5 both pass.
        passes_fold_change = abs(fold_change) >= fold_change_cutoff

        # Significant only if it passes BOTH; then the sign of the fold
        # change says which way it went. (A fold change of 0 can never
        # pass cutoff 2, so "not above 0" here always means below 0.)
        if passes_fdr and passes_fold_change and fold_change > 0:
            directions.append(higher_label)
        elif passes_fdr and passes_fold_change:
            directions.append(lower_label)
        else:
            directions.append(not_significant_label)

    # Add the labels as a new column. The list is in the same order as
    # the table's rows, so each label lands on its own gene.
    results["direction"] = directions

    # Smallest adjusted p-value first (ties broken by the raw p-value);
    # genes with no padj go to the bottom. Then save the table.
    # "output / file name" joins the folder and the file name.
    results = results.sort_values(["padj", "pvalue"])
    results.to_csv(output / differential_expression_file, sep="\t")

    # ===============================================================
    # STEP 3: Make two expression tables for the plots
    # ===============================================================

    # Table A: normalized counts = raw count / the sample's size factor.
    # A deeply sequenced sample (like clone 20) has a bigger size
    # factor, so its counts are divided by more. After this, a gene's
    # numbers can be compared fairly across samples.
    normalized = layer_to_table(model, "normed_counts", metadata.index, counts.index)
    normalized.to_csv(output / normalized_counts_file, sep="\t")

    # Table B: VST ("variance-stabilizing transformation") expression.
    # In raw counts, highly expressed genes also vary the most, so they
    # would dominate heatmaps and correlations. VST is a log-like
    # transform that evens this out, so every gene counts about equally.
    # use_design=False: don't use the WT/KO groups here, so the
    # transform can't make the groups look more different than they are.
    model.vst(use_design=False)
    transformed = layer_to_table(model, "vst_counts", metadata.index, counts.index)
    transformed.to_csv(output / vst_expression_file, sep="\t")

    # ===============================================================
    # STEP 4: Save the sample QC table and the significant genes
    # ===============================================================

    # Add each sample's size factor to the sample table. Near 1 means
    # "average depth"; the deepest sample gets the largest (for our
    # data, clone 20's, about 1.19). index_label puts the word "sample"
    # at the top of the first column.
    metadata["size_factor"] = model.obs["size_factors"]
    metadata.to_csv(output / sample_qc_file, sep="\t", index_label=sample_column)

    # Keep only the genes that passed both cutoffs, and save them.
    significant = results[results["direction"] != not_significant_label]
    significant.to_csv(output / significant_genes_file, sep="\t")

    # e.g. "Genes passing both thresholds: 848"
    print("Genes passing both thresholds:", len(significant))

    # Send three tables back to run_analysis.py:
    #   results     -> test results for every gene in the model
    #   normalized  -> normalized counts (genes x samples)
    #   transformed -> VST expression (genes x samples)
    return results, normalized, transformed
