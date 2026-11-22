"""
===================================================================
run_paper_gsea.py
===================================================================

WHAT THIS FILE DOES
    It runs a SECOND GSEA, set up like the paper's Fig. 2F, so that
    figure can be redrawn and compared term by term:

      - GO Molecular Function terms (the paper's collection), instead
        of the GO Biological Process terms of our main GSEA
      - GO terms with up to 3,000 genes, because some of the paper's
        terms (like "RNA binding") have more than 1,000 genes, and our
        main GSEA stops at 500

    It reuses the ranked gene list from the main GSEA
    (09_gsea_gene_ranking.rnk), so the ONLY things that change are the
    GO collection and the size limits. Any difference between the two
    GSEA runs therefore comes from the gene sets, not the genes.

WHERE IT FITS
    Step 10. It runs only if the paper_reference folder exists, just
    before compare_to_paper.py and plot_paper_figures.py, which both
    use its results.

INPUT
    output   the results folder (it reads 09_gsea_gene_ranking.rnk)

OUTPUT FILE
    12_fig2F_gsea_go_mf.tsv   results for every GO term tested

    It needs internet, to download the GO collection. Without it, this
    step is skipped with a message and the rest of the analysis still
    runs; only the Fig. 2F plot is missing.

PYTHON IDEAS USED IN THIS FILE
    None    a special value meaning "nothing". This function returns
            None when the GSEA couldn't run, so the next steps can tell.
===================================================================
"""

import gseapy as gp
import pandas as pd

from settings import paper_gsea_gene_sets    # "GO_Molecular_Function_2025"
from settings import paper_gsea_min_size     # 15
from settings import paper_gsea_max_size     # 3000
from settings import permutations            # 1000
from settings import seed                    # 42


def run_paper_gsea(output):

    # ===============================================================
    # STEP 1: Read the ranked gene list that run_gsea.py saved
    # ===============================================================

    # The .rnk file has two tab-separated columns and no header row:
    #   gene name <tab> Wald statistic
    #
    #   header=None            -> the first line is data, not names
    #   names=[...]            -> so we give the two columns names
    ranking = pd.read_csv(
        output / "09_gsea_gene_ranking.rnk",
        sep="\t",
        header=None,
        names=["gene_name", "stat"],
    )

    # ===============================================================
    # STEP 2: Run GSEA with the paper's settings
    # ===============================================================

    # Tell the user what is about to run. Two f-strings next to each
    # other inside print( ) are joined into one line of text.
    print(f"Paper-style GSEA: testing {paper_gsea_gene_sets} "
          f"(gene sets of {paper_gsea_min_size} to {paper_gsea_max_size} genes).")

    # If the download fails (for example, no internet), we skip this
    # step instead of stopping the whole analysis. That's the
    # difference from run_gsea.py, which stops: the main GSEA is
    # essential, this one is extra.
    try:
        # Exactly the same call as in run_gsea.py; only the GO
        # collection and the size limits are different.
        enrichment = gp.prerank(
            rnk=ranking,
            gene_sets=paper_gsea_gene_sets,
            min_size=paper_gsea_min_size,
            max_size=paper_gsea_max_size,
            permutation_num=permutations,
            seed=seed,
            outdir=None,
            no_plot=True,
        )
    except Exception as error:
        # Print why it didn't run, then leave the function. Returning
        # None ("nothing") means no results were made.
        print(f"Paper-style GSEA skipped: {error}")
        return None

    # ===============================================================
    # STEP 3: Clean up and save the results
    # ===============================================================

    # The same clean-up as in run_gsea.py.
    pathways = enrichment.res2d.copy()

    # Make sure these two columns hold numbers, not text.
    pathways["FDR q-val"] = pd.to_numeric(pathways["FDR q-val"])
    pathways["NES"] = pd.to_numeric(pathways["NES"])

    # Most significant terms first.
    pathways = pathways.sort_values("FDR q-val")

    # Save every term. plot_paper_figures.py and compare_to_paper.py
    # read this file.
    pathways.to_csv(output / "12_fig2F_gsea_go_mf.tsv", sep="\t", index=False)

    return pathways
