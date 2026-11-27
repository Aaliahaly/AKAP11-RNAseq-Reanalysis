"""
===================================================================
plot_pca.py
===================================================================

WHAT THIS FILE DOES
    It makes two plots that compare the SAMPLES with each other
    (not the genes):

    1. A PCA plot (03_pca.png)
       Each sample is described by ~79,000 numbers (one per gene),
       which is impossible to draw. PCA ("principal component
       analysis") finds the two directions in which the samples differ
       the most, and gives each sample a position along them: PC1 and
       PC2. Then each sample can be drawn as one dot.
         - dots close together = samples with similar expression
         - dots far apart      = samples with different expression

    2. A sample correlation heatmap (04_sample_correlation.png)
       A grid comparing every sample with every other sample.
       A number near 1 means the two samples are very similar.

    Both plots answer: "Do the samples group by genotype (WT vs KO),
    and is any sample an outlier?"

WHERE IT FITS
    Step 3, after the model (fit_model.py).

INPUTS
    count_file    the original count table (ALL genes, unfiltered)
    transformed   the VST table from fit_model.py (filtered genes)
    metadata      each sample's group
    colors        a color for each group (WT blue, KO red)
    output        the results folder

OUTPUT FILES
    03_pca.png, 03_pca_scores.tsv    the PCA plot and each sample's
                                     PC1/PC2 position
    04_sample_correlation.png        the correlation heatmap

PYTHON IDEAS USED IN THIS FILE
    fig, ax = plt.subplots(...)
                  make an empty figure (fig) with one drawing area (ax)
    for name, row in table.iterrows():
                  go through a table one row at a time; each time,
                  "name" is the row's label and "row" its values
===================================================================
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# seaborn: a library built on matplotlib that makes nicer statistical
# plots (scatter plots by group, heatmaps) with less code.
import seaborn as sns

# adjust_text: moves text labels apart so they don't overlap.
from adjustText import adjust_text

# PCA from scikit-learn, Python's main machine-learning library.
from sklearn.decomposition import PCA

from save_plot import save_plot
from settings import condition_column               # "condition"
from settings import sample_column                  # "sample"
from settings import label_line_color               # gray for label lines
from settings import pca_scores_file                # "03_pca_scores.tsv"
from settings import pca_plot_file                  # "03_pca.png"
from settings import sample_correlation_plot_file   # "04_sample_correlation.png"


def plot_pca(count_file, transformed, metadata, colors, output):

    # ===============================================================
    # STEP 1: Read the original count table (every gene, no filtering)
    # ===============================================================

    # We read the count table again here, because this PCA uses ALL
    # 78,932 genes, not only the 26,288 that passed the filter.
    table = pd.read_csv(count_file, sep="\t")

    # The sample names are the row labels of the metadata table:
    # WT_rep1, WT_rep2, WT_rep3, AKAP11_KO_clone16, ...
    samples = list(metadata.index)

    # Keep only the sample columns, in the same order as the metadata.
    # (This drops the gene_id and gene_name columns.)
    counts = table[samples]

    # ===============================================================
    # STEP 2: Turn counts into log2 counts per million (log2 CPM)
    # ===============================================================

    # Some samples were sequenced more deeply than others (clone 20 has
    # about twice as many reads), so EVERY gene has more reads in them.
    # To compare samples fairly, we rescale each sample as if it had
    # exactly one million reads: "counts per million" (CPM).

    # Add up each column. axis=0 means "add DOWN each column", which
    # gives one total per sample.
    sample_totals = counts.sum(axis=0)

    # We are about to divide by these totals, and dividing by zero is
    # impossible, so check every total is above zero first.
    for sample_name in samples:
        if sample_totals[sample_name] <= 0:
            raise ValueError("Every sample needs a positive count total for PCA.")

    # Divide every count by its own sample's total.
    # axis="columns" lines up each total with the column that has the
    # same sample name. The result is each gene's FRACTION of the
    # sample's reads, e.g. 0.0003.
    fraction_of_sample = counts.div(sample_totals, axis="columns")

    # Multiply by one million to get counts per million (0.0003 -> 300).
    # (The underscores in 1_000_000 are just for readability.)
    counts_per_million = fraction_of_sample * 1_000_000

    # Gene expression ranges from 0 to hundreds of thousands of CPM.
    # Without a log, the few hugely expressed genes would decide the
    # whole PCA. log2 squeezes big numbers together:
    #   log2(1) = 0,  log2(1,024) = 10,  log2(1,048,576) = 20
    #
    # log2(0) is impossible (minus infinity), so we add 1 first.
    # A gene with 0 CPM then becomes log2(0 + 1) = 0.
    counts_per_million_plus_one = counts_per_million + 1
    log_expression = np.log2(counts_per_million_plus_one)

    # ===============================================================
    # STEP 3: Run PCA
    # ===============================================================

    # PCA treats each ROW as one item to place, so we need one row per
    # sample. .T swaps rows and columns: 6 rows (samples) x 78,932
    # columns (genes).
    pca_input = log_expression.T

    # Ask for 2 principal components:
    #   PC1 = the single direction in which the samples differ the most
    #   PC2 = the next biggest difference, independent of PC1
    # Why 2? A flat plot has exactly two axes, so 2 is all we can draw.
    # (This is not a tunable analysis choice; the plot needs it.)
    pca = PCA(n_components=2)

    # fit_transform() does two things at once:
    #   "fit"       -> find PC1 and PC2 from the data
    #   "transform" -> give each sample its position (PC1, PC2)
    # The result is a grid with 6 rows and 2 columns.
    coordinates = pca.fit_transform(pca_input)

    # Wrap the grid in a table, with sample names and column names.
    scores = pd.DataFrame(
        coordinates,
        index=samples,
        columns=["PC1", "PC2"],
    )

    # Add each sample's group; the plot uses it to color the dots.
    scores[condition_column] = metadata[condition_column]

    # Save the positions. index_label names the first column "sample".
    scores.to_csv(output / pca_scores_file, sep="\t", index_label=sample_column)

    # ===============================================================
    # STEP 4: Draw the PCA plot
    # ===============================================================

    # An empty figure, 8 inches wide and 5 inches tall.
    #   fig -> the whole figure (what we save)
    #   ax  -> the drawing area inside it (where the dots go)
    fig, ax = plt.subplots(figsize=(8, 5))

    # Draw one dot per sample.
    #   data=scores     -> take the values from this table
    #   x="PC1"/y="PC2" -> which columns go on which axis
    #   hue=...         -> color the dots by group (the "condition" column)
    #   palette=colors  -> use our colors: {"WT": "steelblue", "KO": "firebrick"}
    #   s=95            -> dot size
    sns.scatterplot(
        data=scores,
        x="PC1",
        y="PC2",
        hue=condition_column,
        palette=colors,
        s=95,
        ax=ax,
    )

    # Write each sample's name next to its dot.
    # We keep every text label in a list, so adjust_text() can move
    # them apart afterwards if they overlap.
    labels = []

    # One loop per sample. sample_name is the row label ("WT_rep1"),
    # row holds that row's values (PC1, PC2, condition).
    for sample_name, row in scores.iterrows():
        x_position = row["PC1"]
        y_position = row["PC2"]

        # ax.text(x, y, text) writes text at that point of the plot.
        label = ax.text(x_position, y_position, sample_name, fontsize=9)
        labels.append(label)

    # explained_variance_ratio_ tells us how much of all the variation
    # between samples each PC captures, as a fraction. For our data:
    # PC1 = 0.469 (46.9%) and PC2 = 0.308 (30.8%).
    variance = pca.explained_variance_ratio_
    pc1_fraction = variance[0]   # [0] = the first item (Python counts from 0)
    pc2_fraction = variance[1]   # [1] = the second item

    # ":.1%" formats a fraction as a percentage with 1 decimal place,
    # so 0.469 is written as "46.9%".
    ax.set_xlabel(f"PC1 ({pc1_fraction:.1%})")
    ax.set_ylabel(f"PC2 ({pc2_fraction:.1%})")

    # Title, e.g. "PCA of all 78,932 genes". ":," adds the comma.
    number_of_genes = len(counts)
    ax.set_title(f"PCA of all {number_of_genes:,} genes")

    # Leave 20% empty space around the dots so labels near the edge
    # still fit inside the plot. (Visual only.)
    ax.margins(0.2)

    # Move overlapping labels apart. The dots themselves don't move.
    #   x=, y=      -> where the dots are, so labels avoid covering them
    #   arrowprops  -> draw a thin gray line from a moved label back
    #                  to its dot ("-" = a plain line, no arrowhead)
    #   iter_lim    -> try at most 200 rounds of moving
    adjust_text(
        labels,
        x=scores["PC1"],
        y=scores["PC2"],
        ax=ax,
        arrowprops={"arrowstyle": "-", "color": label_line_color},
        iter_lim=200,
    )

    save_plot(fig, output / pca_plot_file)

    # ===============================================================
    # STEP 5: Draw the sample correlation heatmap
    # ===============================================================

    # This plot uses the filtered, variance-stabilized (VST) table from
    # fit_model.py, not the raw counts. It is a separate function
    # (below) to keep each part short.
    plot_sample_correlation(transformed, output)


def plot_sample_correlation(transformed, output):

    # Each column of "transformed" is one sample.
    # .corr() compares every pair of columns and gives their Pearson
    # correlation, a number between -1 and 1:
    #   1    -> the two samples rise and fall together perfectly
    #   0    -> no relationship
    # A sample compared with itself always gives exactly 1.
    # The result is a 6 x 6 table: every sample against every sample.
    correlations = transformed.corr()

    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw the 6 x 6 table as colored squares.
    #   annot=True     -> write the number inside each square
    #   fmt=".3f"      -> with 3 decimal places (0.984)
    #   cmap="Blues"   -> light blue for lower, dark blue for higher
    #   square=True    -> make every cell a square
    #   cbar_kws=...   -> the title of the color bar on the right
    sns.heatmap(
        correlations,
        annot=True,
        fmt=".3f",
        cmap="Blues",
        square=True,
        ax=ax,
        cbar_kws={"label": "Pearson correlation"},
    )

    ax.set_title("Sample correlations across filtered genes")

    save_plot(fig, output / sample_correlation_plot_file)
