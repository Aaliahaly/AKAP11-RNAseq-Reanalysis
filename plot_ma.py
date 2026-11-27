"""
===================================================================
plot_ma.py
===================================================================

WHAT THIS FILE DOES
    It makes the MA plot (06_ma.png). Every dot is one gene:

        x-axis (left/right) = how highly the gene is expressed on
                              average (mean normalized count, log scale)
        y-axis (up/down)    = how much it changed (log2 fold change)

    Dots above the middle line are higher in KO; dots below it are
    higher in WT.

    The name comes from the two axes in the original version of this
    plot: M (the log ratio, our fold change) against A (the average).

WHY MAKE IT, IF WE HAVE A VOLCANO PLOT?
    It shows something the volcano plot hides: whether changes depend
    on how highly a gene is expressed. Genes with very low counts can
    show huge fold changes just by chance (2 reads vs 0 reads is an
    "infinite" change), so the low-expression end of an MA plot always
    fans out. A healthy analysis has significant genes spread across
    the whole range, not only at the low end.

WHERE IT FITS
    Called at the end of plot_volcano.py, which passes it the gene
    table it has already prepared.

INPUTS
    plot_data          the gene table from plot_volcano.py
    direction_colors   the color for each direction label
    direction_labels   the legend text for each direction label
    output             the results folder

OUTPUT FILE
    06_ma.png

PYTHON IDEAS USED IN THIS FILE
    for name, group in table.groupby("column"):
        split a table into smaller tables, one per value of a column,
        and loop over them. "name" is the value, "group" is the
        smaller table with just those rows.
===================================================================
"""

import matplotlib.pyplot as plt

from save_plot import save_plot
from settings import cutoff_line_color   # gray
from settings import ma_plot_file        # "06_ma.png"


def plot_ma(plot_data, direction_colors, direction_labels, output):

    # An empty figure, 8 inches wide and 5 inches tall.
    fig, ax = plt.subplots(figsize=(8, 5))

    # .groupby("direction") splits the genes into groups by their
    # direction label, and the loop runs once per group:
    #   direction -> the group's label, e.g. "Higher in comparison"
    #   genes     -> a table with only the genes in that group
    # The groups come out in alphabetical order of their labels.
    for direction, genes in plot_data.groupby("direction"):

        # x-axis: the gene's average normalized count over all
        # samples (DESeq2 calls this "baseMean").
        mean_count = genes["baseMean"]

        # y-axis: the gene's log2 fold change, KO vs WT.
        fold_change = genes["log2FoldChange"]

        # Look up this group's color and legend text, so they match
        # the volcano plot exactly.
        color = direction_colors[direction]
        label = direction_labels[direction]

        # Draw the dots.
        #   s=10      -> dot size
        #   alpha=0.6 -> 60% opaque, so overlapping dots stay visible
        ax.scatter(
            mean_count,
            fold_change,
            color=color,
            label=label,
            s=10,
            alpha=0.6,
        )

    # Average counts range from below 10 to over 1,000,000. On a normal
    # axis almost every gene would be squashed at the left edge. A log
    # scale gives each factor of 10 the same width (10, 100, 1,000...),
    # so both ends are visible.
    ax.set_xscale("log")

    ax.set_xlabel("Mean normalized count")
    ax.set_ylabel(r"$\log_2$ fold change")
    ax.set_title("MA plot")

    # A flat line at 0 means "no change": KO and WT the same.
    # (0 is not a setting: it is the definition of "no change" on a
    # log scale, because log2(1) = 0.)
    ax.axhline(0, color=cutoff_line_color, linewidth=0.8)

    # Show the legend, without a box around it.
    ax.legend(frameon=False)

    save_plot(fig, output / ma_plot_file)
