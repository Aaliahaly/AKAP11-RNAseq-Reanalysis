"""
===================================================================
run_gsea.py
===================================================================

WHAT THIS FILE DOES
    It runs GSEA (gene set enrichment analysis) with GSEApy.

    The idea in three steps:
      1. Put ALL genes in one ordered list, from "most up in KO" at
         the top to "most up in WT" (most down in KO) at the bottom.
      2. Take a GO term: a named group of genes that share a job,
         like "Translation" (~200 genes that make proteins).
      3. Walk down the list. If that term's genes bunch up near the
         top or the bottom, instead of being scattered randomly, the
         term is "enriched": its genes shift together.

    Unlike ORA (run_ora.py), GSEA uses every gene, not only the 848
    significant ones. So it can find a pathway whose genes each change
    only a little, but all in the same direction.

WHERE IT FITS
    Step 7, after the plots of single genes.

INPUTS
    results     the test results from fit_model.py
    gene_sets   which GO collection to test (from settings.py)
    output      the results folder

OUTPUT FILES
    09_gsea_gene_ranking.rnk   the ordered gene list
    09_gsea_pathways.tsv       results for every GO term tested

RETURNS (to run_analysis.py)
    pathways     the results table (also saved as 09_gsea_pathways.tsv)
    enrichment   GSEApy's full result, which also holds each term's
                 running enrichment score curve; plot_gsea_enrichment.py
                 draws the classic GSEA plots from it

HOW TO READ THE RESULTS
    NES (normalized enrichment score)
        positive -> the term's genes lean toward the top (up in KO)
        negative -> they lean toward the bottom (down in KO)
        The further from 0, the stronger.
    FDR q-val
        how significant the term is, corrected for testing thousands
        of terms. Below 0.05 counts as significant.

PYTHON IDEAS USED IN THIS FILE
    try: ... except Exception as error: ...
        "try" runs the code under it. If anything inside goes wrong,
        Python jumps to "except" instead of crashing, and "error"
        holds the problem's description.
===================================================================
"""

# gseapy ("gp"): a Python version of the GSEA method.
import gseapy as gp
import numpy as np
import pandas as pd

from settings import gsea_min_size    # 15
from settings import gsea_max_size    # 500
from settings import permutations     # 1000
from settings import seed             # 42
from settings import gsea_ranking_file    # "09_gsea_gene_ranking.rnk"
from settings import gsea_results_file    # "09_gsea_pathways.tsv"


def run_gsea(results, gene_sets, output):

    # ===============================================================
    # STEP 1: Choose which genes can go into the ranking
    # ===============================================================

    # We order genes by their Wald statistic (the "stat" column from
    # fit_model.py). It combines the size of the change with how sure
    # we are about it:
    #   big positive number -> strongly higher in KO
    #   big negative number -> strongly higher in WT
    #   near 0              -> little evidence of any change
    # This is better for ranking than the fold change alone, because a
    # big fold change from a noisy, low-count gene gets a small stat.

    # Drop genes that are missing a stat or a p-value.
    ranking = results.dropna(subset=["stat", "pvalue"])

    # GSEA matches genes to GO terms by gene NAME, so a gene without a
    # name can't be matched to anything. Drop those genes.
    has_gene_name = ranking["gene_name"] != ""
    ranking = ranking[has_gene_name]

    # Drop genes whose stat is infinite (not a normal number), which
    # can't be placed in an ordered list.
    stat_is_finite = np.isfinite(ranking["stat"])
    ranking = ranking[stat_is_finite]

    # ---------------------------------------------------------------
    # Remove repeated gene names
    # ---------------------------------------------------------------

    # A few gene names belong to more than one gene ID. GSEA needs
    # each name only once in the list, so we keep the copy with the
    # highest expression, which is usually the main, real gene.

    # Sort from highest to lowest average count (baseMean).
    # kind="stable" keeps the previous order for genes with exactly
    # the same value, so the result is the same every run.
    ranking = ranking.sort_values("baseMean", ascending=False, kind="stable")

    # drop_duplicates("gene_name") keeps only the FIRST row for each
    # name. After the sort above, that's the most highly expressed copy.
    ranking = ranking.drop_duplicates("gene_name")

    # ===============================================================
    # STEP 2: Put the genes in order and save the ranked list
    # ===============================================================

    # Sort by stat, biggest first: genes most up in KO at the top,
    # genes most up in WT at the bottom.
    ranking = ranking.sort_values("stat", ascending=False, kind="stable")

    # GSEA only needs two columns: the gene name and its score.
    ranking = ranking[["gene_name", "stat"]]

    # Save the list in the standard ".rnk" format used by GSEA tools:
    # two tab-separated columns, no header row, no row labels.
    #   index=False  -> don't write the gene IDs
    #   header=False -> don't write the column names
    # Example lines:
    #   OTX1     16.42
    #   CXCL12   11.15
    # The file name is in settings.py, because run_paper_gsea.py reads
    # this same file later.
    ranking.to_csv(
        output / gsea_ranking_file,
        sep="\t",
        index=False,
        header=False,
    )

    # ===============================================================
    # STEP 3: Run GSEA
    # ===============================================================

    # "try" means: attempt the code below. If anything goes wrong (for
    # example, the GO collection can't be downloaded), jump to the
    # "except" part and show a clear message, instead of crashing with
    # a long, confusing error.
    try:
        # gp.prerank() = GSEA on a list that is already ranked.
        enrichment = gp.prerank(
            # The ranked gene list we just made.
            rnk=ranking,

            # The GO collection to test: an Enrichr library name
            # (downloaded) or a .gmt file name (offline).
            gene_sets=gene_sets,

            # Only test GO terms that share 15 to 500 genes with our
            # list (gsea_min_size and gsea_max_size in settings.py).
            min_size=gsea_min_size,
            max_size=gsea_max_size,

            # To judge whether a term is really enriched, GSEA shuffles
            # the gene labels 1,000 times (permutations in settings.py)
            # and reruns the test each time.
            # That shows what results look like by pure chance.
            permutation_num=permutations,

            # Fix the random shuffles (seed in settings.py), so every
            # run gives the same results.
            seed=seed,

            # Don't let GSEApy write its own files or plots: we save
            # our own table below, and plot_pathways.py draws the plot.
            outdir=None,
            no_plot=True,
        )

    except Exception as error:
        # Build the message. f"..." fills in {error} and {output};
        # "\n" starts a new line.
        message = (
            f"GSEA failed: {error}\n"
            f"The other results are in {output}. "
            "Check the pathway collection and gene identifiers."
        )
        # SystemExit stops the whole program and prints the message.
        # Everything saved before this point is kept.
        raise SystemExit(message) from error

    # ===============================================================
    # STEP 4: Clean up and save the GSEA results
    # ===============================================================

    # GSEApy keeps its results table in enrichment.res2d.
    # .copy() gives us our own table to change.
    pathways = enrichment.res2d.copy()

    # GSEApy can store these two columns as text. pd.to_numeric() makes
    # sure they are numbers, so sorting and plotting work correctly.
    #   FDR q-val -> the corrected significance (smaller = stronger)
    #   NES       -> the normalized enrichment score
    #                (+ leans toward KO, - leans toward WT)
    pathways["FDR q-val"] = pd.to_numeric(pathways["FDR q-val"])
    pathways["NES"] = pd.to_numeric(pathways["NES"])

    # Put the most significant GO terms first.
    pathways = pathways.sort_values("FDR q-val")

    # Save every tested term. index=False: GSEApy's row numbers aren't
    # useful, so we leave them out.
    pathways.to_csv(output / gsea_results_file, sep="\t", index=False)

    # Send two things back to run_analysis.py:
    #   pathways   -> the table, for plot_pathways.py (bar chart)
    #   enrichment -> GSEApy's full result, for plot_gsea_enrichment.py
    #                 (the classic enrichment plots need its curves)
    return pathways, enrichment
