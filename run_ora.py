"""
===================================================================
run_ora.py
===================================================================

WHAT THIS FILE DOES
    It runs ORA (over-representation analysis) with Reactome pathways,
    repeating the analysis in the paper's Fig. 2E.

    The idea, with a bag of marbles:
      - The BACKGROUND is every gene we tested: all the marbles in
        the bag (about 26,000).
      - Our SIGNIFICANT genes are the marbles we pulled out (848).
      - A PATHWAY is one set of genes: say, all the red marbles.

    If we pulled out many more red marbles than chance would predict,
    that pathway is "over-represented" among our significant genes.

    The chance of that is calculated with the hypergeometric test,
    which is exactly the "marbles from a bag" calculation.

WHERE IT FITS
    Step 8, after GSEA.

INPUTS
    results    the test results from fit_model.py
    output     the results folder
    (and the pathway file named in settings.py: reactome_human.gmt)

OUTPUT FILE
    10_reactome_ora.tsv   one row per pathway tested

THE NUMBERS FOR EACH PATHWAY, WITH THE MARBLE NAMES
    N  = genes in the background             (marbles in the bag)
    n  = significant genes                   (marbles pulled out)
    K  = the pathway's genes in the background  (red marbles in the bag)
    k  = the pathway's genes among the significant ones
                                              (red marbles pulled out)

    By chance we'd expect about n x K / N red marbles. The
    hypergeometric test gives the probability of getting k or more.

PYTHON IDEAS USED IN THIS FILE
    set         like a list, but each item appears only once, and
                checking "is x in it?" is very fast
    continue    skip the rest of this loop round and go to the next
    with open(file) as f:
                open a file, and close it automatically at the end
    a tuple     a fixed pair of values written with round brackets,
                like (0.03, 7). Python sorts pairs by the first value,
                and uses the second value to break ties.
===================================================================
"""

import pandas as pd

# The hypergeometric distribution from SciPy, Python's main
# scientific library. It does the "marbles from a bag" calculation.
from scipy.stats import hypergeom

from settings import ora_gene_sets    # "reactome_human.gmt"
from settings import ora_min_size     # 2
from settings import ora_max_size     # 2000


def read_gmt(gmt_file):
    """Read a .gmt pathway file into a dictionary.

    A .gmt file is a plain text file with one pathway per line, and
    tabs between the parts:

        pathway name <tab> description <tab> gene 1 <tab> gene 2 ...

    For example:

        GPCR Ligand Binding R-HSA-500792 <tab> <tab> ADRA1A <tab> ADRB2 ...

    We return a dictionary of pathway name -> list of its genes:

        {"GPCR Ligand Binding R-HSA-500792": ["ADRA1A", "ADRB2", ...], ...}
    """

    # Start with an empty dictionary and fill it one line at a time.
    pathways = {}

    # "with open(...) as file" opens the file for reading. The file is
    # closed again automatically when this block ends.
    with open(gmt_file) as file:

        # Looping over a file gives one line of text at a time.
        for line in file:

            # Each line ends with a line break ("\n"). Remove it, so
            # the last gene name doesn't end with an invisible "\n".
            line = line.rstrip("\n")

            # Split the line at every tab, giving a list of pieces:
            #   ["GPCR Ligand Binding R-HSA-500792", "", "ADRA1A", ...]
            parts = line.split("\t")

            # The first piece (position 0) is the pathway name. The
            # second (position 1) is a description, usually empty.
            # parts[2:] means "every piece from position 2 to the end":
            # the genes.
            name = parts[0]
            genes = parts[2:]

            # Some lines end with extra tabs, which create empty
            # "genes" (""). Keep only real gene names.
            real_genes = []
            for gene in genes:
                if gene != "":
                    real_genes.append(gene)

            # Store this pathway's genes under its name.
            pathways[name] = real_genes

    return pathways


def run_ora(results, output):

    # ===============================================================
    # STEP 1: Make the background list and the significant-gene list
    # ===============================================================

    # The BACKGROUND is every gene that went into the model and has a
    # gene name (pathways list genes by name, so a gene without one
    # can't be matched). The paper also used its tested genes as the
    # background, which is the fair choice: a gene we couldn't measure
    # could never have been "pulled out of the bag".
    #
    # A set keeps each name only once, and checking whether a name is
    # in it is very fast, which matters below.
    background = set()

    for gene_name in results["gene_name"]:
        if gene_name != "":
            background.add(gene_name)

    # Our SIGNIFICANT genes: the ones that passed both cutoffs in
    # fit_model.py.
    is_significant = results["direction"] != "Below thresholds"
    significant_genes = results[is_significant]

    significant = set()

    for gene_name in significant_genes["gene_name"]:
        if gene_name != "":
            significant.add(gene_name)

    # N = genes in the background  (marbles in the bag)
    # n = significant genes        (marbles pulled out)
    N = len(background)
    n = len(significant)

    print(f"ORA: testing {n:,} significant genes against {N:,} background genes.")

    # ===============================================================
    # STEP 2: Test every pathway
    # ===============================================================

    # Read the pathway file into a dictionary (see read_gmt above).
    pathways = read_gmt(ora_gene_sets)

    # We collect one row of results per pathway in this list. Each row
    # is a small dictionary; at the end the list becomes a table.
    rows = []

    # Looping over a dictionary gives its keys: the pathway names.
    for pathway_name in pathways:

        pathway_genes = pathways[pathway_name]

        # Only count this pathway's genes that are in our background.
        genes_in_background = []
        for gene in pathway_genes:
            if gene in background:
                genes_in_background.append(gene)

        # A pathway file can list the same gene twice in one pathway,
        # so turn the list into a set to count each gene once.
        unique_genes_in_background = set(genes_in_background)

        # K = how many of the pathway's genes are in the bag
        # (how many red marbles there are).
        K = len(unique_genes_in_background)

        # Skip pathways that are too small or too big. "continue"
        # jumps straight to the next pathway.
        if K < ora_min_size or K > ora_max_size:
            continue

        # Find the pathway's genes that are also significant.
        overlap_genes = []
        for gene in unique_genes_in_background:
            if gene in significant:
                overlap_genes.append(gene)

        # Sort them alphabetically, so the output is tidy and the same
        # every run.
        overlap_genes.sort()

        # k = how many significant genes are in this pathway
        # (how many red marbles we pulled out).
        k = len(overlap_genes)

        # The hypergeometric test: the probability of pulling out k OR
        # MORE red marbles just by luck.
        #
        # hypergeom.sf(x, ...) gives the probability of MORE than x.
        # "More than k - 1" is the same as "k or more", so we pass k - 1.
        # The other three numbers describe the bag:
        #   N = marbles in the bag, K = red ones, n = marbles pulled out.
        p_value = hypergeom.sf(k - 1, N, K, n)

        # Fold enrichment: how many times more pathway genes we found
        # than chance would predict.
        #   share of OUR significant genes that are in the pathway
        #   divided by
        #   share of ALL background genes that are in the pathway
        # Example: 52 of 848 significant genes (6.1%) are GPCR
        # signalling genes, but only about 2% of the background is,
        # so the fold enrichment is about 3.
        share_in_significant = k / n
        share_in_background = K / N
        fold_enrichment = share_in_significant / share_in_background

        # Store this pathway's results as one dictionary (one row).
        # ";".join(list) glues the gene names into one piece of text
        # with ";" between them: "ADRA1A;ADRB2;...".
        row = {
            "Term": pathway_name,
            "significant_genes_in_term": k,
            "background_genes_in_term": K,
            "fold_enrichment": fold_enrichment,
            "p_value": p_value,
            "genes": ";".join(overlap_genes),
        }
        rows.append(row)

    # If no pathway had the right size, there is nothing to correct or
    # save. Return an empty table; run_analysis.py then skips the plot.
    if len(rows) == 0:
        print("ORA skipped because no pathway had the right number of genes.")
        return pd.DataFrame()

    # Turn the list of row dictionaries into a table: each dictionary
    # becomes one row, and its keys become the column names.
    ora = pd.DataFrame(rows)

    # ===============================================================
    # STEP 3: Correct the p-values for testing many pathways
    # ===============================================================

    # We tested about 1,600 pathways. With that many tests, some will
    # get small p-values by luck. The Benjamini-Hochberg method adjusts
    # the p-values for this, giving an FDR ("false discovery rate")
    # for each pathway. See benjamini_hochberg() below for the steps.
    ora["fdr"] = benjamini_hochberg(list(ora["p_value"]))

    # Put the most significant pathways first: smallest FDR, then
    # smallest raw p-value to break ties. kind="stable" keeps the
    # original order of any rows that are still tied.
    ora = ora.sort_values(["fdr", "p_value"], kind="stable")

    # Put the columns in a sensible order.
    ora = ora[
        [
            "Term",
            "significant_genes_in_term",
            "background_genes_in_term",
            "fold_enrichment",
            "p_value",
            "fdr",
            "genes",
        ]
    ]

    # Save the table. index=False leaves out pandas' row numbers.
    ora.to_csv(output / "10_reactome_ora.tsv", sep="\t", index=False)

    # Count the pathways with FDR below 0.05:
    #   ora["fdr"] < 0.05  -> True/False per pathway
    #   .sum()             -> counts the Trues (True counts as 1)
    passes_fdr = ora["fdr"] < 0.05
    number_passing = passes_fdr.sum()

    print(f"ORA: {number_passing} Reactome pathways have FDR < 0.05.")

    return ora


def benjamini_hochberg(p_values):
    """Turn a list of p-values into Benjamini-Hochberg FDR values.

    The steps are:
      1. Rank the p-values from smallest (rank 1) to largest (rank m).
      2. For each one, FDR = p-value x m / rank.
      3. Going from the largest p-value down to the smallest, make sure
         an FDR is never bigger than the FDR of the next-larger p-value.
      4. An FDR can never be more than 1.

    A small example with m = 4 p-values:

        p-value   rank   p x 4 / rank   after step 3
        0.01      1      0.04           0.04
        0.02      2      0.04           0.04
        0.03      3      0.04           0.04
        0.50      4      0.50           0.50
    """

    # m = how many p-values (how many tests).
    m = len(p_values)

    # ---------------------------------------------------------------
    # Step 1: find the order of the p-values, smallest first.
    # ---------------------------------------------------------------

    # We need to know WHERE each p-value sits in the original list, so
    # we can put its FDR back in the same place at the end.
    #
    # Make a list of (p-value, position) pairs, e.g.
    #   [(0.03, 0), (0.50, 1), (0.01, 2), (0.02, 3)]
    pairs = []
    for position in range(m):
        pairs.append((p_values[position], position))

    # Sort the pairs. Python sorts pairs by their first value (the
    # p-value); if two p-values are equal, it uses the second value
    # (the position), so equal p-values keep their original order:
    #   [(0.01, 2), (0.02, 3), (0.03, 0), (0.50, 1)]
    pairs.sort()

    # Keep just the positions, in sorted order: [2, 3, 0, 1].
    # order[0] is the position of the smallest p-value, order[1] the
    # second smallest, and so on.
    order = []
    for pair in pairs:
        order.append(pair[1])

    # ---------------------------------------------------------------
    # Steps 2 to 4
    # ---------------------------------------------------------------

    # A list of m zeros to fill in. [0.0] * 4 makes [0.0, 0.0, 0.0, 0.0].
    fdr = [0.0] * m

    # Start at 1.0: this also does step 4 (an FDR never above 1).
    smallest_so_far = 1.0

    # Walk from the LARGEST p-value (rank m) down to the smallest
    # (rank 1). range(m, 0, -1) counts m, m-1, ..., 2, 1.
    for rank in range(m, 0, -1):

        # Python lists start at 0, so the p-value with rank r is at
        # order[r - 1].
        position = order[rank - 1]

        # Step 2: FDR = p-value x m / rank.
        value = p_values[position] * m / rank

        # Step 3: never bigger than the FDR we found just before
        # (for the next-larger p-value).
        if value < smallest_so_far:
            smallest_so_far = value

        # Put the result back in the p-value's original position.
        fdr[position] = smallest_so_far

    return fdr
