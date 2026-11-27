"""
===================================================================
plot_gsea_enrichment.py
===================================================================

WHAT THIS FILE DOES
    It draws classic GSEA enrichment plots (09_gsea_enrichment_plots.png)
    for the top GO terms in each direction: the figure style of the
    original GSEA software, used in most papers that report GSEA.

    Each panel shows one GO term, in three stacked parts:

      1. The green curve: the "running enrichment score". GSEA walks
         down our ranked gene list (most up in KO on the left, most
         down in KO on the right). Each time it meets a gene of this
         term, the curve steps UP; for every other gene it drifts a
         little DOWN. The highest (or lowest) point of the curve is the
         term's enrichment score (ES).
           - peak near the LEFT  -> the term's genes are mostly up in KO
           - dip near the RIGHT  -> the term's genes are mostly down in KO
      2. The black tick marks: where each gene of the term sits in the
         ranked list. Ticks bunched at one end = a strong result.
      3. The colored bar: the ranked list itself, red where genes are
         higher in KO, blue where they are lower in KO.

    The left column shows the top terms LOWER in KO, the right column
    the top terms HIGHER in KO. Each title gives the term's NES and
    FDR q-value.

WHERE IT FITS
    Step 7, right after run_gsea.py (and next to plot_pathways.py,
    which draws the summary bar chart of the same GSEA).

INPUTS
    pathways           the GSEA table from run_gsea.py
                       (already sorted, most significant first)
    enrichment         GSEApy's full result, from run_gsea.py. It holds
                       each term's curve and tick positions.
    reference_group    "WT"
    comparison_group   "KO"
    output             the results folder

OUTPUT FILE
    09_gsea_enrichment_plots.png

PYTHON IDEAS USED IN THIS FILE
    plt.subplots(rows, columns, gridspec_kw={"height_ratios": [...]})
                    make a grid of drawing areas, some taller than
                    others
    axes[row, column]
                    the drawing area in that row and column of the grid
    list * 3        repeat a list: [3, 1] * 2 -> [3, 1, 3, 1]
    textwrap.fill(text, 40)
                    break a long text into lines of at most 40 characters
===================================================================
"""

import textwrap

import matplotlib.pyplot as plt
import numpy as np

from save_plot import save_plot
from settings import gsea_enrichment_plot_terms   # 3 terms per direction
from settings import gsea_collection_label        # e.g. "GO Biological Process"
from settings import enrichment_curve_color       # green
from settings import cutoff_line_color            # gray
from settings import term_id_markers              # [" R-HSA-", "(GO:"]
from settings import gsea_enrichment_plot_file    # "09_gsea_enrichment_plots.png"


def short_term_name(term):
    """Remove the database ID from a term name, e.g.
    "Cytoplasmic Translation (GO:0002181)" -> "Cytoplasmic Translation".
    The ID markers are listed in settings.py (term_id_markers)."""

    for marker in term_id_markers:
        term = term.split(marker)[0]
    return term.strip()


def draw_one_term(curve_ax, ticks_ax, bar_ax, term_row, enrichment,
                  reference_group, comparison_group):
    """Draw the three parts of one enrichment plot."""

    # ---------------------------------------------------------------
    # The data for this term, from GSEApy
    # ---------------------------------------------------------------

    term = term_row["Term"]

    # enrichment.results is a dictionary: term name -> its details.
    #   "RES"  -> the running enrichment score: one number per gene in
    #             the ranked list (the green curve)
    #   "hits" -> the positions in the ranked list of this term's genes
    #             (the black ticks)
    details = enrichment.results[term]
    running_score = details["RES"]
    hit_positions = details["hits"]

    # enrichment.ranking holds the score of every gene in the ranked
    # list (the Wald statistic), from most up in KO to most down.
    ranking = enrichment.ranking

    # x positions 0, 1, 2, ... one per gene in the ranked list.
    positions = np.arange(len(running_score))

    # ---------------------------------------------------------------
    # Part 1: the running enrichment score curve
    # ---------------------------------------------------------------

    curve_ax.plot(positions, running_score, color=enrichment_curve_color, linewidth=2)

    # A thin line at 0, so you can see whether the curve goes up or down.
    curve_ax.axhline(0, color=cutoff_line_color, linewidth=0.8)

    # The title: the term name (broken onto two lines if long), then
    # its NES and FDR. GSEApy reports an FDR of 0 when it is smaller
    # than its permutations can measure, so we label it "estimated"
    # instead of claiming it is exactly zero.
    name = textwrap.fill(short_term_name(term), 40)
    if term_row["FDR q-val"] == 0:
        fdr_text = "FDR ≈ 0 (estimated)"
    else:
        fdr_text = f"FDR = {term_row['FDR q-val']:.2g}"
    curve_ax.set_title(f"{name}\nNES = {term_row['NES']:.2f}, {fdr_text}", fontsize=10)

    curve_ax.set_ylabel("Enrichment score")
    curve_ax.set_xlim(0, len(positions))

    # The x-axis numbers are shown on the colored bar below instead.
    curve_ax.set_xticks([])

    # ---------------------------------------------------------------
    # Part 2: a black tick for every gene of the term
    # ---------------------------------------------------------------

    # vlines() draws vertical lines: one at each hit position, from
    # the bottom (0) to the top (1) of this small drawing area.
    ticks_ax.vlines(hit_positions, 0, 1, color="black", linewidth=0.5)
    ticks_ax.set_xlim(0, len(positions))
    ticks_ax.set_ylim(0, 1)

    # This area only shows ticks, so hide its axes and grid.
    ticks_ax.axis("off")

    # ---------------------------------------------------------------
    # Part 3: the ranked list as a red-to-blue bar
    # ---------------------------------------------------------------

    # We color the bar by each gene's PLACE in the list, not by its raw
    # score: a few genes have huge scores, and coloring by score would
    # make almost the whole bar pale. Instead:
    #   - genes with a positive score (up in KO) fade from 1 (deep red)
    #     at the top of the list down to 0 (white)
    #   - genes with a negative score (down in KO) go from 0 (white)
    #     to -1 (deep blue) at the bottom of the list
    # So white marks exactly where the list crosses from "up" to "down".
    # np.linspace(a, b, n) makes n evenly spaced numbers from a to b.
    number_up = int((ranking > 0).sum())
    number_down = len(ranking) - number_up
    up_part = np.linspace(1, 0, number_up)
    down_part = np.linspace(0, -1, number_down)

    # np.concatenate() joins the two parts into one list, and
    # reshape(1, -1) turns it into a grid of 1 row (-1 means "as many
    # columns as needed"), which is what imshow() draws as colors.
    bar_values = np.concatenate([up_part, down_part]).reshape(1, -1)

    # cmap="RdBu_r" goes red (1) -> white (0) -> blue (-1).
    # aspect="auto" stretches the bar to fill its area.
    bar_ax.imshow(bar_values, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)

    # Label the two ends of the bar with the group names (from the
    # command line), like "Up in KO" ... "Down in KO".
    bar_ax.set_yticks([])
    bar_ax.set_xticks([0, len(positions) - 1])
    bar_ax.set_xticklabels([f"Up in {comparison_group}", f"Down in {comparison_group}"], fontsize=9)

    # Line up the end labels with the ends of the bar.
    tick_labels = bar_ax.get_xticklabels()
    tick_labels[0].set_horizontalalignment("left")
    tick_labels[1].set_horizontalalignment("right")
    bar_ax.grid(False)


def plot_gsea_enrichment(pathways, enrichment, reference_group, comparison_group, output):

    # ===============================================================
    # STEP 1: Choose the GO terms to draw
    # ===============================================================

    # run_gsea.py already sorted the table so the most significant
    # terms come first. A negative NES means the term's genes are
    # lower in KO; a positive NES means higher in KO.
    # .head(n) keeps the first n rows (3, from settings.py).
    lower_terms = pathways[pathways["NES"] < 0].head(gsea_enrichment_plot_terms)
    higher_terms = pathways[pathways["NES"] > 0].head(gsea_enrichment_plot_terms)

    # The grid needs as many rows as the longer of the two lists.
    number_of_rows = max(len(lower_terms), len(higher_terms))

    if number_of_rows == 0:
        print("GSEA enrichment plots skipped because no GO terms were tested.")
        return

    # ===============================================================
    # STEP 2: Make the grid of drawing areas
    # ===============================================================

    # Each term needs three stacked areas: the curve (tall), the ticks
    # (short) and the colored bar (thin), plus an empty gap before the
    # next term, with room for the bar's labels and the next title.
    # Their heights relative to each other (visual only):
    #   [3, 0.6, 0.3, 1.6]  ->  curve, ticks, bar, gap
    # "* number_of_rows" repeats that pattern once per row of terms.
    heights_for_one_term = [3, 0.6, 0.3, 1.6]
    height_ratios = heights_for_one_term * number_of_rows

    # 2 columns: lower in KO on the left, higher in KO on the right.
    # The figure is about 3.5 inches tall per row of terms (visual only).
    fig, axes = plt.subplots(
        len(height_ratios),
        2,
        figsize=(12, 4 * number_of_rows),
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.05, "wspace": 0.25},
    )

    # ===============================================================
    # STEP 3: Draw each term in its place
    # ===============================================================

    # The two columns, each with its list of terms and a heading.
    columns = [
        (0, lower_terms, f"Lower in {comparison_group} (higher in {reference_group})"),
        (1, higher_terms, f"Higher in {comparison_group}"),
    ]

    for column, terms, heading in columns:

        # Put this column's terms in a plain list, so we can take them
        # by number (term_rows[0], term_rows[1], ...). .iterrows()
        # gives each row's label (not needed) and the row itself.
        term_rows = []
        for row_label, term_row in terms.iterrows():
            term_rows.append(term_row)

        for row in range(number_of_rows):

            # The 4 grid rows that belong to term number "row":
            # row 0 uses grid rows 0-3, row 1 uses 4-7, and so on.
            first = row * len(heights_for_one_term)
            curve_ax = axes[first, column]
            ticks_ax = axes[first + 1, column]
            bar_ax = axes[first + 2, column]
            gap_ax = axes[first + 3, column]

            # The gap is only empty space.
            gap_ax.axis("off")

            # If this column has fewer terms than the other, leave the
            # spare areas empty.
            if row >= len(term_rows):
                curve_ax.axis("off")
                ticks_ax.axis("off")
                bar_ax.axis("off")
                continue

            draw_one_term(curve_ax, ticks_ax, bar_ax, term_rows[row], enrichment,
                          reference_group, comparison_group)

        # The column heading, just above the first panel's title.
        # xy=(0.5, 1.3) means: half-way across the panel, and 1.3
        # panel-heights up from its bottom edge.
        axes[0, column].annotate(heading, xy=(0.5, 1.3), xycoords="axes fraction",
                                 ha="center", fontsize=12, weight="bold")

    # One title for the whole figure, naming the GO collection.
    fig.suptitle(f"GSEA enrichment plots: top {gsea_collection_label} terms in each direction",
                 fontsize=14, y=0.97)

    save_plot(fig, output / gsea_enrichment_plot_file)
