"""
===================================================================
plot_expression_heatmap.py
===================================================================

WHAT THIS FILE DOES
    It makes a heatmap of the top 25 significant genes
    (07_de_heatmap.png).

        each row    = one gene
        each column = one sample
        each color  = whether the gene is higher (red) or lower (blue)
                      in that sample than the gene's OWN average

    It answers: "Do the top genes change consistently in all three
    knockout clones, or is one clone driving the result?"

WHERE IT FITS
    Step 5, after the volcano and MA plots.

INPUTS
    results            the test results from fit_model.py
                       (already sorted, smallest padj first)
    transformed        the VST expression table from fit_model.py
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
    usual.

PYTHON IDEAS USED IN THIS FILE
    max(a, b)   the larger of two numbers
    .loc[list_of_labels]
                pick rows by their labels, e.g. a list of gene IDs
===================================================================
"""

import matplotlib.pyplot as plt
import seaborn as sns

from save_plot import save_plot


def plot_expression_heatmap(results, transformed, number_of_genes, output):

    # ===============================================================
    # STEP 1: Choose which genes to show
    # ===============================================================

    # fit_model.py already sorted the results by adjusted p-value, so
    # the strongest genes are at the top of the table.

    # Keep the genes that passed both cutoffs...
    is_significant = results["direction"] != "Below thresholds"
    significant_genes = results[is_significant]

    # ...and take the first 25 of them.
    genes = significant_genes.head(number_of_genes)
    title = "Top genes passing the differential expression thresholds"

    # If NO gene passed both cutoffs (possible with another dataset),
    # show the top-ranked genes instead, and change the title so the
    # reader knows these genes are not significant.
    if len(genes) == 0:
        genes_with_padj = results.dropna(subset=["padj"])
        genes = genes_with_padj.head(number_of_genes)
        title = "Top-ranked genes (none passed both differential expression thresholds)"

    # If there are still no genes, there is nothing to draw.
    if len(genes) == 0:
        print("Expression heatmap skipped because no adjusted p-values were available.")
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
    has_spread = sd > 0
    sd = sd[has_spread]

    if len(sd) == 0:
        return

    # Keep the same genes in the expression table.
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

    # Subtract each gene's average from every value in its row.
    # axis=0 lines each average up with the row of the same gene ID.
    centered = expression.sub(gene_average, axis=0)

    # Divide every value in a row by that gene's standard deviation.
    scaled = centered.div(sd, axis=0)

    # Replace the gene IDs (row labels) with the readable labels made
    # in read_counts.py, e.g. "ENSG00000105996.7" -> "HOXA2".
    scaled.index = genes.loc[sd.index, "label"]

    # ===============================================================
    # STEP 4: Draw and save the heatmap
    # ===============================================================

    # Make the figure taller when there are more genes: about 0.28
    # inches per row, plus 2 inches for the title and sample names.
    # max(4, ...) makes sure the figure is at least 4 inches tall,
    # even with only a few genes.
    number_of_rows = len(scaled)
    height = max(4, 0.28 * number_of_rows + 2)

    fig, ax = plt.subplots(figsize=(9, height))

    # Draw the heatmap.
    #   cmap="vlag"       -> blue for low, white for middle, red for high
    #   center=0          -> white means exactly the gene's average
    #   vmin=-2, vmax=2   -> the colors stop getting darker past -2 and
    #                        +2, so one extreme value can't wash out
    #                        all the others
    #   yticklabels=True  -> write every gene's name (by default
    #                        seaborn may skip some to save space)
    sns.heatmap(
        scaled,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        yticklabels=True,
        cbar_kws={"label": "Row z-score (clipped at ±2)"},
    )

    ax.set_title(title)

    # Keep the gene names horizontal (not tilted), so they're easy to
    # read.
    ax.tick_params(axis="y", labelrotation=0)

    save_plot(fig, output / "07_de_heatmap.png")
