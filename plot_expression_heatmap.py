"""
===================================================================
plot_expression_heatmap.py
===================================================================

WHAT THIS FILE DOES
    It makes a clustered heatmap of the top 25 significant genes
    (07_de_heatmap.png).

        each row    = one gene
        each column = one sample
        each color  = whether the gene is higher (red) or lower (blue)
                      in that sample than the gene's OWN average

    "Clustered" means the rows and columns are reordered so that
    similar genes sit next to each other, and similar samples sit next
    to each other. The trees (dendrograms) on the left and at the top
    show how they group: the shorter the branch joining two items, the
    more alike they are. A colored bar above the columns shows each
    sample's group (WT blue, KO red).

    It answers two questions:
      - Do the samples group by genotype? If the column tree splits
        WT from KO, the top genes separate the two groups cleanly.
      - Do the top genes change consistently in all three knockout
        clones, or is one clone driving the result?

WHERE IT FITS
    Step 5, after the volcano and MA plots.

INPUTS
    results            the test results from fit_model.py
                       (already sorted, smallest padj first)
    transformed        the VST expression table from fit_model.py
    metadata           each sample's group (for the color bar)
    colors             a color per group (WT blue, KO red)
    number_of_genes    how many genes to show (25, from settings.py)
    output             the results folder

OUTPUT FILE
    07_de_heatmap.png

WHY Z-SCORES?
    Genes have very different expression levels: one gene might sit
    around 15 on the VST scale and another around 5. Colored by raw
    value, the heatmap would only show which genes are high overall.
    A z-score compares each gene only with itself, so every row shows
    the same thing: in which samples is THIS gene higher or lower than
    usual. Clustering also needs z-scores, so that genes are grouped by
    their PATTERN across samples, not by how highly they're expressed.

PYTHON IDEAS USED IN THIS FILE
    max(a, b)   the larger of two numbers
    .loc[list_of_labels]
                pick rows by their labels, e.g. a list of gene IDs
    .map(dict)  replace each value using a dictionary,
                e.g. "WT" -> "steelblue"
===================================================================
"""

import matplotlib.pyplot as plt
import seaborn as sns

# Patch draws a small colored square, used for the legend that
# explains the group color bar.
from matplotlib.patches import Patch

from save_plot import save_plot
from settings import heatmap_z_limit         # 2
from settings import cluster_method          # "average"
from settings import cluster_metric          # "euclidean"
from settings import condition_column        # "condition"
from settings import not_significant_label   # "Below thresholds"
from settings import heatmap_plot_file       # "07_de_heatmap.png"


def plot_expression_heatmap(results, transformed, metadata, colors, number_of_genes, output):

    # ===============================================================
    # STEP 1: Choose which genes to show
    # ===============================================================

    # fit_model.py already sorted the results by adjusted p-value, so
    # the strongest genes are at the top of the table.

    # Keep the genes that passed both cutoffs...
    is_significant = results["direction"] != not_significant_label
    significant_genes = results[is_significant]

    # ...and take the first number_of_genes of them (25, from
    # heatmap_genes in settings.py).
    genes = significant_genes.head(number_of_genes)
    title = f"Top {len(genes)} differentially expressed genes"

    # If NO gene passed both cutoffs (possible with another dataset),
    # show the top-ranked genes instead, and change the title so the
    # reader knows these genes are not significant.
    if len(genes) == 0:
        genes_with_padj = results.dropna(subset=["padj"])
        genes = genes_with_padj.head(number_of_genes)
        title = f"Top {len(genes)} ranked genes (none passed both thresholds)"

    # Clustering needs at least 2 genes to compare. With fewer, there
    # is nothing to draw.
    if len(genes) < 2:
        print("Expression heatmap skipped because fewer than 2 genes were available.")
        return

    # ===============================================================
    # STEP 2: Get the expression values for those genes
    # ===============================================================

    # genes.index holds the gene IDs we chose. .loc[ ] picks those rows
    # out of the expression table, giving a small 25 x 6 table.
    expression = transformed.loc[genes.index]

    # The standard deviation ("sd", the typical spread) of each gene
    # across the 6 samples. axis=1 means "across each row".
    sd = expression.std(axis=1)

    # In STEP 3 we divide by this spread. A gene with exactly the same
    # value in every sample has a spread of 0, and dividing by 0 is
    # impossible, so keep only genes with a spread above 0.
    # (For real significant genes this never removes anything.)
    sd = sd[sd > 0]
    expression = expression.loc[sd.index]

    # ===============================================================
    # STEP 3: Turn each value into a z-score
    # ===============================================================

    # A z-score says how far a value is from the gene's own average,
    # measured in standard deviations:
    #
    #     z = (value - gene's average) / gene's standard deviation
    #
    #   z =  0  -> exactly this gene's average
    #   z =  2  -> 2 standard deviations above its average
    #   z = -1  -> 1 standard deviation below its average
    #
    # Example: a gene with VST values 8, 8, 8, 4, 4, 4 (WT then KO) has
    # average 6 and sd about 2.2, so its z-scores are about
    # 0.9, 0.9, 0.9, -0.9, -0.9, -0.9.

    # Each gene's average across the samples (axis=1 = across the row).
    gene_average = expression.mean(axis=1)

    # Subtract each gene's average from every value in its row, then
    # divide by that gene's standard deviation. axis=0 lines each
    # average (or sd) up with the row of the same gene ID.
    centered = expression.sub(gene_average, axis=0)
    z_scores = centered.div(sd, axis=0)

    # Replace the gene IDs (row labels) with the readable labels made
    # in read_counts.py, e.g. "ENSG00000105996.7" -> "HOXA2".
    z_scores.index = genes.loc[sd.index, "label"]

    # ===============================================================
    # STEP 4: Make the group color bar for the top of the heatmap
    # ===============================================================

    # One color per sample, from its group: .map(colors) swaps each
    # group name for its color, e.g. "WT" -> "steelblue".
    # The bar's name (shown next to it) is the group column's name.
    sample_colors = metadata[condition_column].map(colors)
    sample_colors.name = condition_column

    # ===============================================================
    # STEP 5: Draw the clustered heatmap
    # ===============================================================

    # Make the figure taller when there are more genes: about 0.3
    # inches per row, plus 3 inches for the trees, title and sample
    # names, but never shorter than 5 inches (visual only).
    height = max(5, 0.3 * len(z_scores) + 3)

    # sns.clustermap() clusters the rows and columns, draws the trees,
    # and draws the heatmap, all in one call.
    #   method, metric   -> how similarity is measured (settings.py:
    #                       average linkage, euclidean distance)
    #   cmap="vlag"      -> blue for low, white for middle, red for high
    #   center=0         -> white means exactly the gene's average
    #   vmin, vmax       -> colors stop getting darker past -2 and +2
    #                       (heatmap_z_limit), so one extreme value
    #                       can't wash out all the others
    #   col_colors       -> the group color bar above the columns
    #   yticklabels=True -> write every gene's name
    #   linewidths, linecolor -> thin white lines between cells, so
    #                       each cell is easy to see
    #   dendrogram_ratio -> how much of the figure the trees take
    #                       (left tree, top tree) (visual only)
    #   cbar_pos         -> where the color scale goes:
    #                       (left, bottom, width, height), as fractions
    #                       of the figure (visual only)
    grid = sns.clustermap(
        z_scores,
        method=cluster_method,
        metric=cluster_metric,
        cmap="vlag",
        center=0,
        vmin=-heatmap_z_limit,
        vmax=heatmap_z_limit,
        col_colors=sample_colors,
        yticklabels=True,
        linewidths=0.5,
        linecolor="white",
        figsize=(9, height),
        dendrogram_ratio=(0.15, 0.12),
        cbar_pos=(0.02, 0.83, 0.03, 0.12),
        cbar_kws={"label": f"Row z-score (±{heatmap_z_limit})"},
    )

    # ---------------------------------------------------------------
    # Finishing touches
    # ---------------------------------------------------------------

    # clustermap makes several drawing areas; ax_heatmap is the heatmap
    # itself. Remove the automatic axis titles ("label", "sample").
    grid.ax_heatmap.set_xlabel("")
    grid.ax_heatmap.set_ylabel("")

    # Gene names in italics (the usual way to write gene symbols) and
    # horizontal; sample names tilted so they don't overlap.
    for gene_label in grid.ax_heatmap.get_yticklabels():
        gene_label.set_style("italic")
    grid.ax_heatmap.tick_params(axis="y", labelrotation=0)
    grid.ax_heatmap.tick_params(axis="x", labelrotation=45)

    # Line up the END of each tilted sample name with its column
    # (ha = "horizontal alignment"), so it's clear which name is which.
    for sample_label in grid.ax_heatmap.get_xticklabels():
        sample_label.set_horizontalalignment("right")

    # A legend for the group color bar: one colored square per group.
    # Looping over the colors dictionary gives each group name
    # (e.g. "WT", "KO") and we look up its color.
    legend_squares = []
    for group in colors:
        legend_squares.append(Patch(color=colors[group], label=group))

    # Put the legend to the right of the top tree, out of the way.
    grid.ax_col_dendrogram.legend(
        handles=legend_squares,
        title=condition_column,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
    )

    # One title over the whole figure. y=1.02 puts it just above the
    # top tree.
    grid.figure.suptitle(title, y=1.02)

    save_plot(grid.figure, output / heatmap_plot_file)
