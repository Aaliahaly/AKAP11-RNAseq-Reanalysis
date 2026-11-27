"""
===================================================================
plot_target_gene.py
===================================================================

WHAT THIS FILE DOES
    It makes a bar plot for ONE gene you choose (for this project,
    AKAP11, the gene that was knocked out), with one bar per sample
    (08_target_gene.png). It also saves that gene's test results
    (08_target_gene_results.tsv).

    For AKAP11 this checks the knockout itself: is AKAP11 lower in
    every knockout clone?

WHERE IT FITS
    Step 6. It only runs if you typed a gene name as the last word of
    the command, e.g.  ... results_folder AKAP11

INPUTS
    results       the test results from fit_model.py
    normalized    the normalized counts from fit_model.py
    metadata      each sample's group
    symbol        the gene name you typed, e.g. "AKAP11"
    colors        a color per group (WT blue, KO red)
    output        the results folder

OUTPUT FILES
    08_target_gene.png            the bar plot
    08_target_gene_results.tsv    the gene's row(s) from the results

PYTHON IDEAS USED IN THIS FILE
    .str.upper()   make every text value in a column UPPERCASE
    .idxmax()      the row label of the largest value in a column
    .map(dict)     replace each value using a dictionary,
                   e.g. "WT" -> "steelblue"
===================================================================
"""

import matplotlib.pyplot as plt

from save_plot import save_plot
from settings import condition_column            # "condition"
from settings import target_gene_plot_file       # "08_target_gene.png"
from settings import target_gene_results_file    # "08_target_gene_results.tsv"


def plot_target_gene(results, normalized, metadata, symbol, colors, output):

    # ===============================================================
    # STEP 1: Find the gene in the results
    # ===============================================================

    # Make both sides UPPERCASE before comparing, so "akap11",
    # "Akap11" and "AKAP11" all match.
    gene_names_upper = results["gene_name"].str.upper()
    symbol_upper = symbol.upper()

    # True for rows whose gene name matches, False for every other row.
    is_match = gene_names_upper == symbol_upper

    # Keep only the matching rows (usually exactly one).
    matches = results.loc[is_match]

    # If nothing matched, the gene was probably removed by the
    # low-count filter in read_counts.py (or the name is misspelled).
    # Print a message and skip this plot.
    if len(matches) == 0:
        print(f"{symbol} was absent from the filtered count table.")
        return

    # Save the test results of every matching row: fold change,
    # p-value, adjusted p-value, and so on.
    matches.to_csv(output / target_gene_results_file, sep="\t")

    # ===============================================================
    # STEP 2: Get the values and colors for the bars
    # ===============================================================

    # Occasionally two gene IDs share one gene name. If that happens,
    # plot the one with the highest average count (baseMean).
    # .idxmax() gives the ROW LABEL (here, the gene ID) of the largest
    # value, e.g. "ENSG00000023516.9" for AKAP11.
    gene_id = matches["baseMean"].idxmax()

    # The sample names, in the order we want the bars (WT first).
    samples = metadata.index

    # From the normalized counts, take this gene's row and those
    # sample columns: six numbers, one per sample.
    expression = normalized.loc[gene_id, samples]

    # Look up each sample's bar color from its group.
    # .map(colors) replaces each value using the colors dictionary:
    #   "WT" -> "steelblue",  "KO" -> "firebrick"
    bar_colors = metadata[condition_column].map(colors)

    # ===============================================================
    # STEP 3: Draw the bars and save the figure
    # ===============================================================

    fig, ax = plt.subplots(figsize=(8, 4))

    # ax.bar(x positions, heights): one bar per sample, as tall as its
    # normalized count, in its group's color.
    ax.bar(samples, expression, color=bar_colors)

    # Tilt the sample names by 30 degrees so the long names don't
    # overlap each other.
    ax.tick_params(axis="x", labelrotation=30)

    ax.set_ylabel("Normalized count")

    # The title shows the gene name and its ID, e.g.
    # "AKAP11 expression (ENSG00000023516.9)".
    ax.set_title(f"{symbol} expression ({gene_id})")

    save_plot(fig, output / target_gene_plot_file)
