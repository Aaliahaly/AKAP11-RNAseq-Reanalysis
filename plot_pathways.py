"""
===================================================================
plot_pathways.py
===================================================================

WHAT THIS FILE DOES
    It makes the GSEA bar plot (09_gsea_pathways.png). Each bar is
    one GO term:

        bars to the RIGHT of zero -> the term's genes lean toward KO
                                     (positive NES, drawn in KO red)
        bars to the LEFT of zero  -> they lean toward WT
                                     (negative NES, drawn in WT blue)

    The longer the bar, the stronger the shift. Each label also shows
    the term's FDR q-value (its corrected significance).

WHERE IT FITS
    Step 7, right after run_gsea.py.

INPUTS
    pathways           the GSEA table from run_gsea.py
                       (already sorted, most significant first)
    reference_group    "WT"
    comparison_group   "KO"
    colors             a color per group
    output             the results folder

OUTPUT FILE
    09_gsea_pathways.png

PYTHON IDEAS USED IN THIS FILE
    text.replace("_", " ")   swap every "_" in a piece of text for " "
    f"{value:.3g}"           write a number with up to 3 significant
                             digits, e.g. 0.0123
===================================================================
"""

import matplotlib.pyplot as plt
import pandas as pd

from save_plot import save_plot
from settings import gsea_top                 # 8 terms per direction
from settings import gsea_collection_label    # e.g. "GO Biological Process"


def plot_pathways(pathways, reference_group, comparison_group, colors, output):

    # ===============================================================
    # STEP 1: Choose the GO terms to show
    # ===============================================================

    # run_gsea.py already sorted the table so the most significant
    # terms (smallest FDR q-value) come first. So taking the "first 8"
    # of a group means taking its 8 most significant terms.

    # Terms that lean toward KO have a positive NES.
    # Keep those rows, then take the top 8.
    leans_toward_comparison = pathways["NES"] > 0
    positive = pathways[leans_toward_comparison].head(gsea_top)

    # Terms that lean toward WT have a negative NES.
    leans_toward_reference = pathways["NES"] < 0
    negative = pathways[leans_toward_reference].head(gsea_top)

    # Stack the two small tables into one (up to 16 rows).
    top = pd.concat([positive, negative])

    # Sort by NES, from most negative to most positive. matplotlib
    # draws the first bar at the BOTTOM, so the most negative (WT)
    # terms end up at the bottom and the most positive (KO) at the top.
    top = top.sort_values("NES")

    # ===============================================================
    # STEP 2: Make a label and pick a color for each bar
    # ===============================================================

    # Two lists, built in the same order as the rows of "top":
    #   labels     -> the text next to each bar
    #   bar_colors -> each bar's color
    labels = []
    bar_colors = []

    # One loop per GO term. row_number is GSEApy's row number (not
    # needed); row holds that term's values.
    for row_number, row in top.iterrows():

        # Some collections write term names with underscores
        # ("CYTOPLASMIC_TRANSLATION"). Replace them with spaces so the
        # names are easier to read.
        term_name = row["Term"]
        term_name = term_name.replace("_", " ")

        q_value = row["FDR q-val"]

        # Add the q-value after the term name.
        #
        # With 1,000 permutations, GSEApy can't measure a q-value
        # smaller than about 0.001. When it's smaller than that, GSEApy
        # writes exactly 0. A q-value of 0 isn't really possible, so we
        # label it "estimated" instead of claiming it's exactly zero.
        if q_value == 0:
            label = f"{term_name} (q=0 (estimated))"
        else:
            label = f"{term_name} (q={q_value:.3g})"

        labels.append(label)

        # KO's color for a positive NES, WT's color for a negative one.
        if row["NES"] > 0:
            bar_colors.append(colors[comparison_group])
        else:
            bar_colors.append(colors[reference_group])

    # ===============================================================
    # STEP 3: Draw and save the plot
    # ===============================================================

    # Make the figure taller when there are more bars: about 0.36
    # inches per bar, plus 1.5 inches for the title and axis.
    # max(4, ...) keeps it at least 4 inches tall.
    number_of_bars = len(top)
    height = max(4, 0.36 * number_of_bars + 1.5)

    fig, ax = plt.subplots(figsize=(9, height))

    # barh() draws HORIZONTAL bars ("h" for horizontal). That leaves
    # room on the left for the long GO term names.
    #   labels      -> one label per bar (the y positions)
    #   top["NES"]  -> each bar's length
    ax.barh(labels, top["NES"], color=bar_colors)

    # A thin line at zero separates the WT side from the KO side.
    ax.axvline(0, color="gray", linewidth=0.8)

    # The x-axis title, e.g.
    # "Normalized enrichment score (− WT, + KO)"
    x_axis_title = (
        "Normalized enrichment score (− "
        + reference_group
        + ", + "
        + comparison_group
        + ")"
    )
    ax.set_xlabel(x_axis_title)

    # The plot title. gsea_collection_label comes from settings.py,
    # e.g. "Top GO Biological Process terms by GSEA FDR".
    ax.set_title("Top " + gsea_collection_label + " terms by GSEA FDR")

    save_plot(fig, output / "09_gsea_pathways.png")
