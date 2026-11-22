"""
===================================================================
plot_volcano.py
===================================================================

WHAT THIS FILE DOES
    It makes the volcano plot (05_volcano.png). Every dot is one gene:

        x-axis (left/right) = how much the gene changed
                              (log2 fold change, KO vs WT)
        y-axis (up/down)    = how confident we are that it changed
                              (-log10 of the adjusted p-value)

    Genes in the top-left and top-right corners changed a lot AND are
    statistically strong. The cloud of dots ends up shaped a bit like
    an erupting volcano, hence the name.

    At the end it also calls plot_ma() (in plot_ma.py) to make the
    MA plot, because both plots use the same prepared gene table.

WHERE IT FITS
    Step 4, after the PCA.

INPUTS
    results            the test results from fit_model.py
    reference_group    "WT"
    comparison_group   "KO"
    colors             a color per group (WT blue, KO red)
    output             the results folder

OUTPUT FILES
    05_volcano.png          the plot
    05_volcano_labels.tsv   the numbers for every gene named on the plot

WHY -log10 OF THE P-VALUE?
    Small p-values mean strong evidence, but on a normal axis 0.001
    and 0.000001 would both sit squashed near zero. -log10 turns them
    into numbers that grow as the evidence gets stronger:
        padj = 0.05     -> -log10 = 1.3   (our cutoff line)
        padj = 0.001    -> -log10 = 3
        padj = 1e-20    -> -log10 = 20

PYTHON IDEAS USED IN THIS FILE
    dictionary { key: value }
                a lookup table, used here to match each direction
                label to a color and a legend text
    .head(n)    the first n rows of a table
===================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from plot_ma import plot_ma
from save_plot import save_plot
from settings import fdr_cutoff            # 0.05
from settings import fold_change_cutoff    # 1.0
from settings import volcano_labels        # 8 names per side


def plot_volcano(results, reference_group, comparison_group, colors, output):

    # ===============================================================
    # STEP 1: Get the numbers for the x-axis and y-axis
    # ===============================================================

    # A few genes have no adjusted p-value or no fold change (NaN),
    # so they can't be placed on the plot. dropna(subset=[...]) drops
    # the rows that are empty in either of these two columns.
    # .copy() gives us our own table, so we can add a column to it.
    genes_with_values = results.dropna(subset=["padj", "log2FoldChange"])
    genes = genes_with_values.copy()

    # If nothing is left, there's nothing to draw: print a message and
    # leave the function early with "return".
    if len(genes) == 0:
        print("Volcano and MA plots skipped because no adjusted p-values were available.")
        return

    # We are about to take log10 of the adjusted p-values, and log10(0)
    # is impossible. A few genes have p-values so tiny the computer
    # stores them as 0, so we raise anything smaller than 1e-300
    # (that's 0.000...0001 with 300 zeros) up to 1e-300.
    # This only affects the plot; the "padj" column itself isn't changed.
    adjusted_p = genes["padj"].clip(lower=1e-300)

    # -log10: smaller p-values become bigger numbers (see the top of
    # this file), so stronger genes sit higher on the plot.
    log10_p = np.log10(adjusted_p)
    genes["minus_log10_fdr"] = -log10_p

    # ===============================================================
    # STEP 2: Draw the dots and the cutoff lines
    # ===============================================================

    # fit_model.py gave every gene one of three "direction" labels.
    # This dictionary gives each label a color:
    direction_colors = {
        "Below thresholds": "#b9bfc6",                      # light gray
        "Lower in comparison": colors[reference_group],     # WT's color (blue)
        "Higher in comparison": colors[comparison_group],   # KO's color (red)
    }

    # And this one gives each label the text shown in the legend.
    # A gene lower in KO is the same thing as higher in WT, so
    # "Lower in comparison" appears as "Higher in WT".
    direction_labels = {
        "Below thresholds": "Below thresholds",
        "Lower in comparison": "Higher in " + reference_group,
        "Higher in comparison": "Higher in " + comparison_group,
    }

    # A wide figure (12 x 7 inches), with room for gene names and a
    # legend on the right.
    fig, ax = plt.subplots(figsize=(12, 7))

    # Draw the three groups of genes one after another, so each group
    # gets its own color and its own line in the legend.
    # Looping over a dictionary goes through its keys:
    # "Below thresholds", then "Lower in comparison", then
    # "Higher in comparison". The gray group is drawn first, so the
    # colored dots end up on top of it.
    for direction in direction_colors:

        # Keep only the genes with this direction.
        is_this_direction = genes["direction"] == direction
        group = genes[is_this_direction]

        x_values = group["log2FoldChange"]
        y_values = group["minus_log10_fdr"]

        # ax.scatter() draws one dot per (x, y) pair.
        ax.scatter(
            x_values,
            y_values,
            color=direction_colors[direction],
            label=direction_labels[direction],
            s=12,        # dot size
            alpha=0.65,  # 65% opaque, so overlapping dots stay visible
        )

    # A dashed horizontal line at the p-value cutoff. The y-axis is
    # -log10(padj), so the line goes at -log10(0.05), which is about
    # 1.3. Dots above it have padj < 0.05.
    fdr_line_height = -np.log10(fdr_cutoff)
    ax.axhline(fdr_line_height, linestyle="--", color="gray", linewidth=0.8)

    # Two dashed vertical lines at the fold-change cutoff: one at -1
    # (halved) and one at +1 (doubled). Colored dots are always outside
    # these two lines AND above the horizontal one.
    left_cutoff = -fold_change_cutoff
    right_cutoff = fold_change_cutoff

    ax.axvline(left_cutoff, linestyle="--", color="gray", linewidth=0.8)
    ax.axvline(right_cutoff, linestyle="--", color="gray", linewidth=0.8)

    # ===============================================================
    # STEP 3: Choose which genes get a name label
    # ===============================================================

    # A True/False column: did the gene pass both cutoffs?
    genes["passes_thresholds"] = genes["direction"] != "Below thresholds"

    # Sort the genes so the best candidates for labels come first:
    #   1st: genes that passed (True) before genes that didn't (False)
    #   2nd: smaller adjusted p-value first
    #   3rd: smaller raw p-value first (breaks ties)
    # ascending=[False, True, True] gives the direction of each of the
    # three sorts. False = biggest first, and True counts as bigger
    # than False, so passing genes come first.
    ranked = genes.sort_values(
        ["passes_thresholds", "padj", "pvalue"],
        ascending=[False, True, True],
    )

    # Left side: genes with a negative fold change (higher in WT).
    # .head(8) takes the first 8 rows, which are the best 8.
    is_left_side = ranked["log2FoldChange"] < 0
    left = ranked[is_left_side].head(volcano_labels)

    # Right side: genes with a positive fold change (higher in KO).
    is_right_side = ranked["log2FoldChange"] > 0
    right = ranked[is_right_side].head(volcano_labels)

    # Stack the two small tables into one table of genes to label.
    candidates = pd.concat([left, right])

    # Write each chosen gene's label at its dot, and keep the text
    # objects in a list for adjust_text() below.
    labels = []

    for gene_id, gene in candidates.iterrows():
        x_position = gene["log2FoldChange"]
        y_position = gene["minus_log10_fdr"]
        gene_label = gene["label"]

        label = ax.text(x_position, y_position, gene_label, fontsize=8)
        labels.append(label)

    # ---------------------------------------------------------------
    # Axis titles and legend
    # ---------------------------------------------------------------

    # matplotlib draws text between $ signs as math, so "$\log_2$"
    # appears as "log" with a small 2 below it.
    # The r before the quotes makes it a "raw" string: the backslash
    # (\) is kept exactly as typed instead of being treated as a
    # special character.
    # The pieces are joined with + into e.g. "log2 fold change (KO / WT)".
    x_axis_title = (
        r"$\log_2$ fold change ("
        + comparison_group
        + " / "
        + reference_group
        + ")"
    )
    ax.set_xlabel(x_axis_title)
    ax.set_ylabel(r"$-\log_{10}$ adjusted p-value")
    ax.set_title("Differential expression volcano plot")

    # Put the legend just outside the plot's right edge.
    #   loc="upper left"          -> which corner of the legend to anchor
    #   bbox_to_anchor=(1.01, 1)  -> where to anchor it: just right of
    #                                the plot (x = 1.01) at the top (y = 1)
    #   frameon=False             -> no box around the legend
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)

    # Leave 12% empty space near the plot's edges.
    ax.margins(0.12)

    # Move overlapping labels apart; the dots stay where they are, and
    # short gray lines connect moved labels to their dots.
    # Only do this if there is at least one label.
    if len(labels) > 0:
        adjust_text(
            labels,
            x=candidates["log2FoldChange"],
            y=candidates["minus_log10_fdr"],
            ax=ax,
            arrowprops={"arrowstyle": "-", "color": "#65717a"},
            iter_lim=250,
        )

    # ===============================================================
    # STEP 4: Save the plot and the list of labeled genes
    # ===============================================================

    save_plot(fig, output / "05_volcano.png")

    # Save the numbers behind every labeled gene, so you can look them
    # up without reading them off the plot.
    columns_to_save = ["gene_name", "log2FoldChange", "padj", "direction"]
    labeled_genes = candidates[columns_to_save]
    labeled_genes.to_csv(output / "05_volcano_labels.tsv", sep="\t")

    # ===============================================================
    # STEP 5: Make the MA plot, with the same genes and colors
    # ===============================================================

    # Pass the prepared table and both dictionaries to plot_ma.py, so
    # the two plots use exactly the same colors and legend text.
    plot_ma(genes, direction_colors, direction_labels, output)
