"""
===================================================================
plot_ora.py
===================================================================

WHAT THIS FILE DOES
    It makes the Reactome ORA bar plot (10_reactome_ora.png), drawn
    like the paper's Fig. 2E.

    Each bar is one Reactome pathway. A LONGER bar means a SMALLER
    FDR, which means stronger evidence that the pathway has more of
    our significant genes than we'd expect by chance.

WHERE IT FITS
    Step 10. run_analysis.py calls it right after run_ora.py, but only
    if run_ora.py tested at least one pathway.

INPUTS
    ora     the table made by run_ora.py, one row per pathway, already
            sorted with the smallest FDR first. Columns used here:
                Term                       the pathway's name
                significant_genes_in_term  how many of our significant
                                           genes are in the pathway
                fdr                        the adjusted p-value
    output  the results folder

OUTPUT FILE
    10_reactome_ora.png

PYTHON IDEAS USED IN THIS FILE
    table.head(10)     the first 10 rows of a table
    text.split(x)[0]   cut a text at x and keep the piece BEFORE it
    text.strip()       remove spaces from the start and end of a text
    -np.log10(x)       turns a small number into a big one:
                       0.1 -> 1, 0.01 -> 2, 0.001 -> 3
===================================================================
"""

import matplotlib.pyplot as plt
import numpy as np

from save_plot import save_plot
from settings import ora_top            # how many pathways to show (10)
from settings import fdr_cutoff         # 0.05
from settings import term_id_markers    # [" R-HSA-", "(GO:"]
from settings import neutral_color      # "slategray"
from settings import cutoff_line_color  # gray
from settings import ora_plot_file      # "10_reactome_ora.png"
from settings import ora_collection_label   # "Reactome"


def remove_term_id(name):
    """Cut a pathway name at its database ID, if it has one.

    For example
        "Extracellular Matrix Organization R-HSA-1474244"
    becomes
        "Extracellular Matrix Organization"

    The ID markers to look for are listed in settings.py
    (term_id_markers), so this works for Reactome and GO names, and
    for any other database you add to that list.
    """

    # Try each marker in turn. text.split(marker) cuts the text at the
    # marker and gives a list of pieces; [0] is the piece BEFORE it.
    # If the marker isn't in the text, split() gives back the whole
    # text as the only piece, so nothing changes.
    for marker in term_id_markers:
        name = name.split(marker)[0]

    # .strip() removes any spaces left at the start or end.
    return name.strip()


def plot_ora(ora, output):

    # ===============================================================
    # STEP 1: Choose the pathways to show
    # ===============================================================

    # run_ora.py already sorted the table, smallest FDR first, so the
    # first ora_top rows (10, from settings.py) are the strongest
    # pathways.
    top = ora.head(ora_top)

    # If there are no rows, there's nothing to draw. "return" leaves
    # the function straight away.
    if len(top) == 0:
        print("ORA plot skipped because no pathways were tested.")
        return

    # Matplotlib draws the FIRST bar at the BOTTOM of the plot. We want
    # the strongest pathway at the TOP, so we reverse the row order:
    # make a list of the row labels, reverse it, and pick the rows in
    # that new order with .loc[].
    reversed_index = list(top.index)
    reversed_index.reverse()
    top = top.loc[reversed_index]

    # ===============================================================
    # STEP 2: Make the labels and the bar lengths
    # ===============================================================

    labels = []        # the text beside each bar
    bar_lengths = []   # how long each bar is

    # Go through the pathways one row at a time.
    for row_number, row in top.iterrows():

        # Pathway names end with a database ID, for example
        #   "Extracellular Matrix Organization R-HSA-1474244"
        # Remove the ID so the label is shorter and easier to read
        # (see remove_term_id() at the top of this file).
        name = remove_term_id(row["Term"])

        # Add the number of our significant genes in the pathway, e.g.
        #   "Extracellular Matrix Organization (34 genes)"
        k = row["significant_genes_in_term"]
        labels.append(f"{name} ({k} genes)")

        # The bar length is -log10(FDR), so a smaller FDR gives a
        # longer bar. For example, FDR 0.001 gives a bar of length 3.
        bar_lengths.append(-np.log10(row["fdr"]))

    # ===============================================================
    # STEP 3: Draw and save the plot
    # ===============================================================

    # Give each bar 0.4 inches of height, plus 1.5 inches for the
    # title and axis, but never make the plot shorter than 4 inches.
    height = max(4, 0.4 * len(top) + 1.5)

    # A figure 9 inches wide; "ax" is the area we draw in.
    fig, ax = plt.subplots(figsize=(9, height))

    # barh = horizontal bars: one label and one length per bar.
    ax.barh(labels, bar_lengths, color=neutral_color)

    # A dashed line at FDR = fdr_cutoff (0.05 in settings.py;
    # -log10(0.05) is about 1.3). Bars that reach past it pass the
    # significance cutoff.
    ax.axvline(-np.log10(fdr_cutoff), linestyle="--", color=cutoff_line_color, linewidth=0.8)

    # The r"..." text with $...$ is written in math style, so it shows
    # as "-log10 FDR" with a proper subscript.
    ax.set_xlabel(r"$-\log_{10}$ FDR")
    # The collection's name comes from settings.py, so the title stays
    # right if you switch to another pathway file.
    ax.set_title("Top " + ora_collection_label + " pathways among significant genes (ORA)")

    # save_plot() saves the figure as a PNG and closes it.
    save_plot(fig, output / ora_plot_file)
