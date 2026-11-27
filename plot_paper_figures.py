"""
===================================================================
plot_paper_figures.py
===================================================================

WHAT THIS FILE DOES
    It redraws the paper's figures from OUR results, panel by panel,
    so they can be put side by side with the published ones.
    (Which paper, and every value taken from it, is read from the
    paper_reference folder; for this project that is
    Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557.)

      Fig. 1B  the knocked-out gene in each clone, relative to WT
      Fig. 2A  volcano plot, with the genes the paper labels
      Fig. 2B  number of significant genes, ours next to the paper's
      Fig. 2C  clustered heatmap of all significant genes
      Fig. 2D  our top down and top up genes, next to the paper's
      Fig. 2E  ORA (already drawn by plot_ora.py, so here we only
               describe it)
      Fig. 2F  GSEA with GO Molecular Function terms, as dot plots

    Some panels (for this paper Fig. 1C, Fig. 1D and Fig. 3) come from
    Western blots, microscope images and DNA methylation sequencing,
    not RNA-seq, so they can't be redrawn from this data. They are
    listed, with the reason, in paper_reference/paper_not_reproducible.tsv
    and the report says so.

    It also writes a report, 12_paper_figures_report.md, that says for
    every panel what matches, what differs, and why.

WHERE IT FITS
    Step 13, the last step. run_analysis.py calls it if the
    paper_reference folder exists. You can also run it on its own:

        python3 plot_paper_figures.py akap11_go_bp_results WT KO

    (the results folder, the reference group, the comparison group).

INPUTS (all read from files, so it can run on its own)
    from the results folder: the 01_, 02_, 10_, 11_ and 12_fig2F_
        files (their names are in settings.py, PART 10)
    from the paper_reference folder: values copied from the paper
        (see paper_reference/README.md)

OUTPUT FILES (names set in settings.py, PART 10)
    12_fig1B_target_gene_relative.png / .tsv
    12_fig2A_volcano.png
    12_fig2B_deg_counts.png
    12_fig2C_heatmap.png
    12_fig2D_top_genes.png / .tsv
    12_fig2F_gsea_go_mf.png
    12_paper_figures_report.md

WHY NO RESULT IS TYPED INTO THE REPORT TEXT
    Every number, gene name and sample name in the report is filled
    in from our results or from the paper_reference files when the
    code runs. Sentences that make a claim ("matches", "differs") are
    only written when the data actually show it. So if the data, a
    setting, or even the paper changes, the report stays correct.

HOW THE FILE IS ORGANIZED
    PART 1   small helper functions
    PART 2   one function per figure panel. Each one draws its figure
             (if it can) and returns two things:
                 notes   -> a list of sentences for the report
                 status  -> e.g. "Reproduced with differences"
    PART 3   plot_paper_figures(): calls every panel function, then
             writes the report
    PART 4   the code that runs when you start this file on its own

PYTHON IDEAS USED IN THIS FILE
    return notes, status
             a function can send back two values at once; the caller
             receives them with  notes, status = some_function(...)
    textwrap.fill(text, 38)
             break a long text into lines of at most 38 characters
    "{:.2e}" a number in scientific notation, e.g. 3.78e-60
    all(list_of_True_False)
             True only if EVERY item in the list is True
===================================================================
"""

import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text

# LinearSegmentedColormap lets us make our own color scale from a
# list of colors (used for the Fig. 2C heatmap).
from matplotlib.colors import LinearSegmentedColormap

# Helpers that live in other files of this project, borrowed here so
# each one is written only once:
#   simplify_name()          makes GO term names comparable
#   read_paper_info()        reads paper_info.tsv
#   read_paper_numbers()     reads paper_summary.tsv
#   read_methods()           reads paper_methods.tsv
#   paper_citation()         builds the citation line
#   fold_change_correlation() the paper-vs-ours correlation
#   remove_term_id()         cuts the database ID off a term name
from compare_to_paper import simplify_name
from compare_to_paper import read_paper_info
from compare_to_paper import read_paper_numbers
from compare_to_paper import read_methods
from compare_to_paper import paper_citation
from compare_to_paper import fold_change_correlation
from plot_ora import remove_term_id
from save_plot import save_plot

# Cutoffs, labels and options (see settings.py for the reasons).
from settings import fdr_cutoff                  # 0.05
from settings import fold_change_cutoff          # 1.0
from settings import permutations                # 1000
from settings import paper_gsea_gene_sets        # "GO_Molecular_Function_2025"
from settings import paper_gsea_collection_label # "GO Molecular Function"
from settings import condition_column            # "condition"
from settings import higher_label                # "Higher in comparison"
from settings import lower_label                 # "Lower in comparison"
from settings import not_significant_label       # "Below thresholds"
from settings import smallest_plotted_p_value    # 1e-300
from settings import heatmap_z_limit             # 2
from settings import cluster_method              # "average"
from settings import cluster_metric              # "euclidean"
from settings import ora_top                     # 10
from settings import ora_collection_label        # "Reactome"
from settings import agree_direction_percent     # 80
from settings import report_examples             # 3

# The paper-style colors (PART 8 of settings.py).
from settings import paper_down_color
from settings import paper_up_color
from settings import paper_gray_color
from settings import paper_wt_color
from settings import paper_ko_color
from settings import paper_highlight_color
from settings import paper_heatmap_colors
from settings import label_line_color

# The paper_reference files this script reads (PART 9).
from settings import paper_genes_file
from settings import paper_volcano_labels_file
from settings import paper_ora_file
from settings import paper_gsea_file
from settings import paper_fig1b_file
from settings import paper_fig2c_columns_file
from settings import paper_not_reproducible_file

# Our results files, read and written here (PART 10).
from settings import sample_qc_file
from settings import differential_expression_file
from settings import normalized_counts_file
from settings import vst_expression_file
from settings import ora_results_file
from settings import ora_plot_file
from settings import gene_comparison_file
from settings import term_comparison_file
from settings import fold_change_comparison_plot_file
from settings import paper_gsea_results_file
from settings import fig1b_table_file
from settings import fig1b_plot_file
from settings import fig2a_plot_file
from settings import fig2b_plot_file
from settings import fig2c_plot_file
from settings import fig2d_table_file
from settings import fig2d_plot_file
from settings import fig2f_plot_file
from settings import figures_report_file


# ###################################################################
# PART 1: Helper functions
# ###################################################################

def read_table(path):
    """Read one of our .tsv results files, keeping empty gene names as ""."""

    # index_col=0 uses the first column (gene IDs) as row labels.
    # keep_default_na=False + na_values=[""] treat ONLY empty cells as
    # missing, so a gene called "NA" is kept as a name.
    table = pd.read_csv(path, sep="\t", index_col=0, keep_default_na=False, na_values=[""])

    # Empty gene names were read as missing (NaN); turn them into "".
    if "gene_name" in table.columns:
        table["gene_name"] = table["gene_name"].fillna("")

    return table


def find_gene(results, symbol):
    """Return the gene ID for a gene name, or None if it isn't there.

    If two gene IDs share the name, return the more highly expressed one.
    """

    # Compare in UPPERCASE, so letter case can't cause a missed match.
    names_upper = results["gene_name"].str.upper()
    is_match = names_upper == symbol.upper()
    matches = results[is_match]

    if len(matches) == 0:
        return None

    # .idxmax() = the row label (gene ID) of the largest baseMean.
    return matches["baseMean"].idxmax()


def samples_in_group(metadata, group):
    """Return the names of the samples in one group, in table order.

    For example samples_in_group(metadata, "WT") ->
        ["WT_rep1", "WT_rep2", "WT_rep3"]
    The group name comes from the command line, never from this file.
    """

    names = []
    for sample in metadata.index:
        if metadata.loc[sample, condition_column] == group:
            names.append(sample)
    return names


def top_gene_names(results, how_many):
    """The names of our most significant genes, e.g. "CYP26C1, HOXA2, ...".

    Used in the report wherever it says which genes come out on top,
    so the names are always the ones in the current results.
    """

    # Keep genes that passed both cutoffs, smallest adjusted p first.
    is_significant = results["direction"] != not_significant_label
    ranked = results[is_significant].sort_values(["padj", "pvalue"])

    # .head(n) = the first n rows; use the plot label (gene name, or
    # the gene ID if there's no name).
    names = list(ranked.head(how_many)["label"])
    return ", ".join(names)


def read_gene_comparison(output):
    """Read compare_to_paper.py's gene-by-gene table, if it exists.

    Returns three numbers about the paper's top genes:
        number found in our results,
        number of those changing in the same direction,
        the correlation of the two sets of fold changes.
    If the file is missing, returns 0, 0 and np.nan ("missing").
    """

    path = output / gene_comparison_file
    if not path.exists():
        return 0, 0, np.nan

    table = pd.read_csv(path, sep="\t")
    found_genes = table[table["found_in_our_results"]]

    number_found = len(found_genes)
    number_same = found_genes["same_direction"].sum()
    correlation = fold_change_correlation(found_genes)
    return number_found, number_same, correlation


def read_sample_labels(reference_folder, samples):
    """Short labels for the samples, e.g. "AKAP11_KO_clone20" -> "KO clone 20".

    The labels come from paper_fig2c_columns.tsv. Any sample not in
    that file keeps its own name as its label, so nothing breaks with
    another dataset.
    """

    labels = {}

    # Start with every sample labeled by its own name...
    for sample in samples:
        labels[sample] = sample

    # ...then use the short label from the file where there is one.
    path = reference_folder / paper_fig2c_columns_file
    if path.exists():
        table = pd.read_csv(path, sep="\t", dtype=str)
        for row_number, row in table.iterrows():
            labels[row["sample"]] = row["label"]

    return labels


# ###################################################################
# PART 2: One function per figure panel
# ###################################################################

# ===================================================================
# Fig. 1B: the knocked-out gene in each clone, relative to WT
# ===================================================================
#
# The paper measured the knocked-out gene (AKAP11, named in
# paper_info.tsv) by RT-qPCR and showed each clone relative to WT
# (WT = 1). We do the same with RNA-seq: each clone's normalized count
# of that gene divided by the WT average.

def plot_fig1b(results, normalized, metadata, info, reference_group, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: Our numbers
    # ---------------------------------------------------------------

    # The gene and the method the paper used, from paper_info.tsv.
    gene = info["fig1b_gene"]
    method = info["fig1b_method"]

    gene_id = find_gene(results, gene)

    # If the gene isn't in the results, return a note and a status, and
    # skip the rest of this function.
    if gene_id is None:
        return [f"{gene} is not in the filtered results, so Fig. 1B was skipped."], "skipped"

    # The gene's normalized count in each sample (one row of the table).
    counts = normalized.loc[gene_id]

    # The reference (WT) samples. The group name comes from the
    # command line, so this works for any group name.
    wt_samples = samples_in_group(metadata, reference_group)

    # The WT average becomes our "1.0".
    wt_mean = counts[wt_samples].mean()

    # ---------------------------------------------------------------
    # STEP 2: The paper's numbers, and its clone order and labels
    # ---------------------------------------------------------------

    # paper_fig1b_qpcr.tsv lists the clones in the paper's bar order,
    # each with its value read off the figure and its bar label:
    #   sample              relative_expression   label
    #   AKAP11_KO_clone21   0.43                  KO clone 21
    paper_table = pd.read_csv(reference_folder / paper_fig1b_file, sep="\t")

    clones = list(paper_table["sample"])
    clone_labels = list(paper_table["label"])

    # Every clone in the file must also be in our data.
    for clone in clones:
        if clone not in counts.index:
            return [f"{clone} (from {paper_fig1b_file}) is not in our samples, "
                    "so Fig. 1B was skipped."], "skipped"

    # The paper's value for each clone, as a dictionary:
    #   {"AKAP11_KO_clone21": 0.43, ...}
    paper = {}
    for row_number, row in paper_table.iterrows():
        paper[row["sample"]] = row["relative_expression"]

    # Ours: each clone relative to WT, e.g. 0.55 = 55% of the WT level.
    ours = {}
    for clone in clones:
        ours[clone] = counts[clone] / wt_mean

    # ---------------------------------------------------------------
    # STEP 3: Save the numbers side by side
    # ---------------------------------------------------------------

    # The WT row is 1.0 in both, by definition (everything is divided
    # by the WT average, and the paper set WT = 1).
    rows = []
    rows.append({
        "sample": f"{reference_group} (mean of {len(wt_samples)})",
        "this_project_rnaseq": 1.0,
        "paper": 1.0,
    })
    for clone in clones:
        rows.append({
            "sample": clone,
            "this_project_rnaseq": round(ours[clone], 3),   # 3 decimal places
            "paper": paper[clone],
        })
    numbers = pd.DataFrame(rows)
    numbers.to_csv(output / fig1b_table_file, sep="\t", index=False)

    # ---------------------------------------------------------------
    # STEP 4: Draw two bars per group: ours (filled) and the paper's
    #         (striped)
    # ---------------------------------------------------------------

    # The bar labels: the reference group's name, then the clones'.
    labels = [reference_group] + clone_labels

    # The bar heights, in the same order as the labels.
    our_values = [1.0]
    paper_values = [1.0]
    for clone in clones:
        our_values.append(ours[clone])
        paper_values.append(paper[clone])

    # Each group sits at x = 0, 1, 2, ... Each bar is 0.38 wide (visual
    # only); our bar goes half a bar to the LEFT of the group's center
    # and the paper's half a bar to the RIGHT, so they sit side by side.
    width = 0.38
    positions = []
    our_positions = []
    paper_positions = []
    for i in range(len(labels)):
        positions.append(i)
        our_positions.append(i - width / 2)
        paper_positions.append(i + width / 2)

    # One color per group: green for WT, blue for every clone.
    # [color] * 3 makes a list of the same color 3 times.
    bar_colors = [paper_wt_color] + [paper_ko_color] * len(clones)

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Our bars: filled with color.
    ax.bar(our_positions, our_values, width, color=bar_colors,
           label="This project: RNA-seq, whole gene")

    # The paper's bars: white with a colored outline and stripes
    # (hatch="//"), so the two kinds are easy to tell apart.
    ax.bar(paper_positions, paper_values, width, color="white",
           edgecolor=bar_colors, hatch="//", linewidth=1.5,
           label=f"Paper: {method} (read off Fig. 1B)")

    # Show each WT replicate as a black dot on the WT bar, like the
    # paper does. All the dots share the WT bar's x position.
    wt_relative = counts[wt_samples] / wt_mean
    wt_dot_x = []
    for sample in wt_samples:
        wt_dot_x.append(our_positions[0])
    ax.scatter(wt_dot_x, wt_relative, color="black", s=18, zorder=3)

    # The top of the y-axis: 30% above the tallest bar or dot, to leave
    # room for the numbers above the bars and for the legend. It is
    # worked out from the data, so no bar can ever be cut off.
    tallest = max(our_values + paper_values + list(wt_relative))
    y_top = tallest * 1.3

    # Write each bar's value just above it (2% of the axis height
    # higher), e.g. "0.55".
    gap = y_top * 0.02
    for i in range(len(labels)):
        ax.text(our_positions[i], our_values[i] + gap, f"{our_values[i]:.2f}", ha="center", fontsize=9)
    for i in range(len(labels)):
        ax.text(paper_positions[i], paper_values[i] + gap, f"{paper_values[i]:.2f}", ha="center", fontsize=9)

    # Put the group names under the group centers.
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)

    ax.set_ylim(0, y_top)
    ax.set_ylabel(f"{gene} relative to {reference_group} ({reference_group} = 1)")
    ax.set_title(f"Fig. 1B: {gene} in each knockout clone")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    save_plot(fig, output / fig1b_plot_file)

    # ---------------------------------------------------------------
    # STEP 5: What does it show? Notes for the report
    # ---------------------------------------------------------------

    # A label for each clone, for the sentences below:
    #   {"AKAP11_KO_clone21": "KO clone 21", ...}
    label_of = {}
    for position in range(len(clones)):
        label_of[clones[position]] = clone_labels[position]

    # Which clone keeps the most of the gene, in our data and in the
    # paper's? Start with the first clone and replace it whenever a
    # later clone is higher.
    our_highest = clones[0]
    for clone in clones:
        if ours[clone] > ours[our_highest]:
            our_highest = clone

    paper_highest = clones[0]
    for clone in clones:
        if paper[clone] > paper[paper_highest]:
            paper_highest = clone

    # Is the gene lower than WT (below 1) in every clone, in both?
    # all() is True only if every item is True.
    our_all_lower = all([ours[clone] < 1 for clone in clones])
    paper_all_lower = all([paper[clone] < 1 for clone in clones])

    # The two lists of values as text: "KO clone 21 = 0.55, ..."
    our_parts = []
    paper_parts = []
    for clone in clones:
        our_parts.append(f"{label_of[clone]} = {ours[clone]:.2f}")
        paper_parts.append(f"{label_of[clone]} = {paper[clone]:.2f}")

    notes = []
    notes.append(f"Our RNA-seq levels relative to {reference_group}: " + ", ".join(our_parts)
                 + f". The paper ({method}): " + ", ".join(paper_parts) + ".")

    # Only claim what the numbers show.
    if our_all_lower and paper_all_lower:
        matches_text = f"Matches: {gene} is lower than in {reference_group} in every knockout clone in both"
        if our_highest == paper_highest:
            matches_text = matches_text + f", and {label_of[our_highest]} keeps the most {gene} in both"
        notes.append(matches_text + ".")
    else:
        notes.append(f"Differs: {gene} is not lower than in {reference_group} in every clone in both "
                     "datasets; check the table "
                     f"{fig1b_table_file}.")

    # Compare the average drop in the two methods.
    our_average = np.mean([ours[clone] for clone in clones])
    paper_average = np.mean([paper[clone] for clone in clones])
    if our_average > paper_average:
        notes.append(
            f"Differs: the drop is smaller in RNA-seq (average {our_average:.2f} vs {paper_average:.2f}). "
            f"The paper measured {method}, the region targeted by the knockout edit, so it sees only intact copies. RNA-seq counts reads along the "
            "whole gene, including the edited copy, which is still transcribed. A heterozygous knockout "
            "also keeps one working copy, so about 0.5 is the expected level for the intact region."
        )
        notes.append("To check the targeted region directly, you could count reads on it in the "
                     "STAR BAM files from Part 1.")
    else:
        notes.append(f"Our average level ({our_average:.2f}) is as low as or lower than the paper's "
                     f"({paper_average:.2f}).")

    status = "Reproduced with differences"
    return notes, status


# ===================================================================
# Fig. 2A: volcano plot with the paper's labeled genes
# ===================================================================
#
# Same idea as our 05_volcano.png, but in the paper's style and with
# the genes the paper labels (paper_volcano_labels.tsv), so you can
# see where they fall in our data.

# A note on "methods" (used by several panels below): it is the table
# paper_methods.tsv as a dictionary, read by read_methods() in
# compare_to_paper.py. Its keys are the step names in the file's first
# column, so methods["Statistics"]["paper"] is the paper's statistics
# method ("limma-voom, robust fit") and methods["Statistics"]
# ["this_project"] is ours ("PyDESeq2").

def plot_fig2a(results, methods, reference_group, comparison_group, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: The x and y values for every gene
    # ---------------------------------------------------------------

    # Drop genes with no adjusted p-value or no fold change: they
    # can't be placed on the plot.
    genes = results.dropna(subset=["padj", "log2FoldChange"]).copy()

    # y = -log10(adjusted p-value). clip() raises values below
    # smallest_plotted_p_value (settings.py) up to it, because
    # log10(0) is impossible.
    clipped_padj = genes["padj"].clip(lower=smallest_plotted_p_value)
    genes["minus_log10_padj"] = -np.log10(clipped_padj)

    # ---------------------------------------------------------------
    # STEP 2: Find the paper's labeled genes in our results
    # ---------------------------------------------------------------

    label_table = pd.read_csv(reference_folder / paper_volcano_labels_file, sep="\t")
    paper_labels = list(label_table["gene_name"])

    labeled_ids = []   # gene IDs we can plot
    missing = []       # not in our results at all
    no_padj = []       # in our results, but without an adjusted p-value

    for symbol in paper_labels:
        gene_id = find_gene(genes, symbol)
        if gene_id is not None:
            labeled_ids.append(gene_id)
        elif find_gene(results, symbol) is not None:
            # In our results, but DESeq2 gave it no adjusted p-value.
            no_padj.append(symbol)
        else:
            missing.append(symbol)

    # Nothing to draw if none of the paper's genes are in our results.
    if len(labeled_ids) == 0:
        return ["None of the paper's labeled genes are in our results, so Fig. 2A was skipped."], "skipped"

    # The rows of the genes we can plot.
    labeled = genes.loc[labeled_ids]
    found = len(labeled)

    # How many of them are significant here?
    is_significant = labeled["direction"] != not_significant_label
    significant = is_significant.sum()

    # Each labeled gene's color: blue if down in KO, red if up.
    labeled_colors = []
    for fold_change in labeled["log2FoldChange"]:
        if fold_change < 0:
            labeled_colors.append(paper_down_color)
        else:
            labeled_colors.append(paper_up_color)

    # ---------------------------------------------------------------
    # STEP 3: Draw two panels
    # ---------------------------------------------------------------

    # Left: the whole plot. Right: a zoom on the area where the paper's
    # genes sit, so their names have room.
    # plt.subplots(1, 2, ...) makes 1 row of 2 drawing areas; "axes" is
    # a list of the two. width_ratios makes the right one a bit wider.
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.25]})

    # One color per direction label (labels from settings.py).
    color_of = {
        not_significant_label: paper_gray_color,
        lower_label: paper_down_color,
        higher_label: paper_up_color,
    }

    # Draw the same dots in both panels.
    for ax in axes:

        # All genes, one color per direction. Looping over the
        # dictionary goes through its keys in the order written above,
        # so the gray dots are drawn first and the colored ones on top.
        for direction in color_of:
            is_this_direction = genes["direction"] == direction
            group = genes[is_this_direction]
            ax.scatter(group["log2FoldChange"], group["minus_log10_padj"], color=color_of[direction],
                       s=8, alpha=0.6)

        # The paper's labeled genes on top: bigger, with a black
        # outline. zorder=3 draws them above everything else.
        ax.scatter(labeled["log2FoldChange"], labeled["minus_log10_padj"],
                   color=labeled_colors,
                   s=40, edgecolors="black", linewidths=0.8, zorder=3)

        # The cutoff lines, from settings.py: the adjusted p-value
        # cutoff, and minus/plus the fold-change cutoff.
        ax.axhline(-np.log10(fdr_cutoff), linestyle="--", color="black", linewidth=0.8)
        ax.axvline(-fold_change_cutoff, linestyle="--", color="black", linewidth=0.8)
        ax.axvline(fold_change_cutoff, linestyle="--", color="black", linewidth=0.8)

        ax.set_xlabel("log2(FC)")
        ax.set_ylabel("-log10(adj p-value)")

    # axes[0] is the left panel, axes[1] the right one.
    axes[0].set_title("All genes (outlined = labeled in the paper)")

    # The zoom area is worked out from the labeled genes themselves:
    #   x: from the lowest to the highest labeled fold change, plus
    #      15% of that range on each side for the names
    #   y: from 0 to 15% above the highest labeled gene
    # So the zoom always fits the genes, whatever the data.
    lowest_x = labeled["log2FoldChange"].min()
    highest_x = labeled["log2FoldChange"].max()
    x_padding = (highest_x - lowest_x) * 0.15

    zoom = axes[1]
    zoom.set_xlim(lowest_x - x_padding, highest_x + x_padding)
    zoom.set_ylim(0, labeled["minus_log10_padj"].max() * 1.15)
    zoom.set_title("Zoom on the paper's labeled genes")

    # Write each labeled gene's name, in its color, in bold italics.
    texts = []
    for gene_id, gene in labeled.iterrows():
        if gene["log2FoldChange"] < 0:
            color = paper_down_color
        else:
            color = paper_up_color
        text = zoom.text(gene["log2FoldChange"], gene["minus_log10_padj"], gene["gene_name"],
                         fontsize=9, style="italic", color=color, weight="bold")
        texts.append(text)

    # Move the names apart so they don't overlap.
    #   expand  -> how much empty space to keep around each name
    #   iter_lim=1000 -> try harder than usual, as there are many names
    adjust_text(texts, x=list(labeled["log2FoldChange"]), y=list(labeled["minus_log10_padj"]), ax=zoom,
                arrowprops={"arrowstyle": "-", "color": label_line_color, "lw": 0.6},
                expand=(1.4, 1.6), iter_lim=1000)

    # suptitle = one title over both panels. The group names come from
    # the command line.
    fig.suptitle(f"Fig. 2A: {comparison_group} vs {reference_group}", y=1.0)
    save_plot(fig, output / fig2a_plot_file)

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # Are the labeled genes on the same side (sign of the fold change)
    # as in the paper? The paper's fold changes are in paper_degs.tsv.
    paper_genes = pd.read_csv(reference_folder / paper_genes_file, sep="\t")
    paper_fold_change = {}
    for row_number, row in paper_genes.iterrows():
        paper_fold_change[row["gene_name"].upper()] = row["log2FoldChange"]

    checked = 0
    same_side = 0
    for gene_id, gene in labeled.iterrows():
        name = gene["gene_name"].upper()
        if name in paper_fold_change:
            checked = checked + 1
            # np.sign() gives -1, 0 or +1; equal signs = same side.
            if np.sign(gene["log2FoldChange"]) == np.sign(paper_fold_change[name]):
                same_side = same_side + 1

    # Our most significant gene: the first row after sorting by padj.
    # .iloc[0] means "the row at position 0", the first one.
    sorted_genes = genes.sort_values("padj")
    our_top = sorted_genes.iloc[0]

    # The paper's most significant gene, from its top-gene table.
    # .idxmin() = the row label of the smallest padj.
    paper_top = paper_genes.loc[paper_genes["padj"].idxmin()]
    paper_top_height = -np.log10(paper_top["padj"])

    # The paper-vs-ours correlation from compare_to_paper.py.
    number_found, number_same, correlation = read_gene_comparison(output)

    notes = []
    notes.append(
        f"Matches: the same layout, cutoffs and colors. Of the {len(paper_labels)} genes the paper labels, "
        f"{found} can be plotted here and {significant} of those are significant; "
        f"{same_side} of the {checked} that are in the paper's gene table are on the same side as in the paper."
    )
    notes.append(
        f"Differs: our y-axis reaches {our_top['minus_log10_padj']:.0f} (top gene {our_top['gene_name']}), "
        f"the paper's about {paper_top_height:.0f} (top gene {paper_top['gene_name']}). Our model gives much "
        "smaller p-values for the strongest genes, so the genes the paper labels sit lower down in our plot, "
        f"under genes such as {top_gene_names(results, report_examples)}."
    )
    notes.append(
        f"Why: the paper's model ({methods['Statistics']['paper']}) shrinks extreme statistics; ours "
        f"({methods['Statistics']['this_project']}, a Wald test) does not. The two methods agree on which "
        f"genes change and in which direction ({number_same} of {number_found} of the paper's top genes go the "
        f"same way, r = {correlation:.2f}), but not on how small the smallest p-values are."
    )
    if len(missing) > 0:
        notes.append(f"Not in our results (renamed between {methods['Annotation']['paper']} (paper) and "
                     f"{methods['Annotation']['this_project']} (ours), or removed by our low-count filter): "
                     + ", ".join(missing) + ".")
    if len(no_padj) > 0:
        notes.append("In our results but not plotted, because DESeq2 gave them no adjusted p-value (it leaves out "
                     "some low-count or outlier genes when correcting p-values): " + ", ".join(no_padj) + ".")

    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2B: number of significant genes
# ===================================================================
#
# The paper shows one stacked bar: up genes (red) with down genes
# (blue) on top. We draw the paper's bar and ours side by side.

def plot_fig2b(results, paper_numbers, methods, reference_group, comparison_group, output):

    # ---------------------------------------------------------------
    # STEP 1: The numbers
    # ---------------------------------------------------------------

    is_down = results["direction"] == lower_label
    our_down = is_down.sum()

    is_up = results["direction"] == higher_label
    our_up = is_up.sum()

    # The paper's numbers, from paper_summary.tsv. int() turns 421.0
    # into 421, because these are counts of genes.
    paper_down = int(paper_numbers["lower_in_comparison"])
    paper_up = int(paper_numbers["higher_in_comparison"])

    # ---------------------------------------------------------------
    # STEP 2: Draw two stacked bars
    # ---------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(4.5, 5))

    names = ["Paper", "This project"]
    ups = [paper_up, our_up]
    downs = [paper_down, our_down]

    # The red "up" part first, from 0 up to the number of up genes.
    ax.bar(names, ups, color=paper_up_color, label="Upregulated", width=0.55)

    # The blue "down" part on top of it: bottom=ups starts each blue
    # bar where its red bar ends.
    ax.bar(names, downs, bottom=ups, color=paper_down_color, label="Downregulated", width=0.55)

    # The y-axis goes 15% above the taller bar, so the "total" text
    # fits above it. Worked out from the data, so it fits any numbers.
    tallest = max(paper_up + paper_down, our_up + our_down)
    y_top = tallest * 1.15

    # Write the numbers: in the middle of each colored part (white
    # text), and the total just above each bar (2% of the axis height
    # above it). len(names) is the number of bars (2).
    for i in range(len(names)):
        ax.text(i, ups[i] / 2, str(ups[i]), ha="center", va="center", color="white", fontsize=11)
        ax.text(i, ups[i] + downs[i] / 2, str(downs[i]), ha="center", va="center", color="white", fontsize=11)
        ax.text(i, ups[i] + downs[i] + y_top * 0.02, f"total {ups[i] + downs[i]}", ha="center", fontsize=10)

    ax.set_ylabel("Number of DEGs")
    ax.set_ylim(0, y_top)
    ax.set_title(f"Fig. 2B: {comparison_group} vs {reference_group}")

    # Put the legend beside the bars, so it can't cover the totals.
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
    save_plot(fig, output / fig2b_plot_file)

    # ---------------------------------------------------------------
    # STEP 3: Notes for the report
    # ---------------------------------------------------------------

    # Did the paper use the same two cutoffs as settings.py? The
    # paper's cutoffs are in paper_summary.tsv.
    paper_fdr = paper_numbers["fdr_cutoff"]
    paper_fold = paper_numbers["fold_change_cutoff"]
    if paper_fdr == fdr_cutoff and paper_fold == fold_change_cutoff:
        cutoff_text = (f"with the same cutoffs (adjusted p < {fdr_cutoff}, "
                       f"|log2FC| ≥ {fold_change_cutoff:g})")
    else:
        cutoff_text = (f"but with different cutoffs (paper: adjusted p < {paper_fdr}, |log2FC| ≥ {paper_fold:g}; "
                       f"ours: adjusted p < {fdr_cutoff}, |log2FC| ≥ {fold_change_cutoff:g})")

    notes = []
    notes.append(f"Matches: {our_down + our_up} significant genes here, {paper_down + paper_up} in the paper, "
                 + cutoff_text + ".")
    notes.append(f"Differs: the split is {our_down} down / {our_up} up here and {paper_down} / {paper_up} in the paper.")

    # Every pipeline step where our method differs from the paper's
    # (from paper_methods.tsv). Looping over the dictionary gives the
    # step names. We leave out steps that are the same in both, and
    # steps that come AFTER the gene list is made (the file's
    # affects_gene_list column says "no" for those, e.g. GSEA), since
    # they can't change the number of significant genes.
    differences = []
    for step in methods:
        ours = methods[step]["this_project"]
        theirs = methods[step]["paper"]
        affects_genes = methods[step]["affects_gene_list"] == "yes"
        if ours != theirs and affects_genes:
            differences.append(f"{step.lower()} ({ours} vs {theirs})")
    notes.append("Why: genes close to the cutoffs can tip either way. These steps differ from the paper's: "
                 + "; ".join(differences)
                 + ". Each one moves a few genes across a cutoff.")
    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2C: clustered heatmap of all significant genes
# ===================================================================
#
# Like our 07_de_heatmap.png, but in the paper's style: ALL
# significant genes instead of the top heatmap_genes (50), only the
# rows (genes) clustered, and the columns kept in the paper's order.

def plot_fig2c(results, vst, metadata, info, methods, reference_group, comparison_group,
               reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: The column order, and z-scores of all significant genes
    # ---------------------------------------------------------------

    is_significant = results["direction"] != not_significant_label
    significant = results[is_significant]

    # The samples in the paper's column order, from
    # paper_fig2c_columns.tsv (for this paper, clone 21 comes before
    # clone 20). If that file is missing, or lists a sample we don't
    # have, use our own sample order instead.
    our_samples = list(metadata.index)
    order = our_samples

    columns_file = reference_folder / paper_fig2c_columns_file
    if columns_file.exists():
        paper_order = list(pd.read_csv(columns_file, sep="\t", dtype=str)["sample"])
        all_present = all([sample in our_samples for sample in paper_order])
        if all_present:
            order = paper_order

    # Short column labels, e.g. "KO clone 16" (see read_sample_labels).
    label_of = read_sample_labels(reference_folder, our_samples)
    column_labels = [label_of[sample] for sample in order]

    # The VST values of the significant genes, in that column order.
    expression = vst.loc[significant.index, order]

    # Row z-scores, as in plot_expression_heatmap.py:
    # (value - gene's average) / gene's standard deviation.
    # First drop any gene with no spread (it would mean dividing by 0).
    spread = expression.std(axis=1)
    expression = expression[spread > 0]
    gene_average = expression.mean(axis=1)
    centered = expression.sub(gene_average, axis=0)
    z = centered.div(expression.std(axis=1), axis=0)

    # ---------------------------------------------------------------
    # STEP 2: Draw the clustered heatmap
    # ---------------------------------------------------------------

    # A color scale like the paper's (colors in settings.py): purple
    # (low), white, red (high).
    colors = LinearSegmentedColormap.from_list("paper", paper_heatmap_colors)

    # sns.clustermap() clusters the rows and draws the heatmap.
    #   col_cluster=False   -> keep the columns in the paper's order
    #   row_cluster=True    -> reorder the genes so similar patterns
    #                          sit together
    #   method, metric      -> how "similar" is measured (cluster_method
    #                          and cluster_metric in settings.py: the
    #                          average straight-line distance)
    #   vmin, vmax          -> the colors stop at -2 and +2
    #                          (heatmap_z_limit in settings.py)
    #   yticklabels=False   -> no gene names (848 would be unreadable)
    #   dendrogram_ratio    -> how much room the trees get (visual)
    #   cbar_pos            -> where the color bar goes: (left, bottom,
    #                          width, height), here outside on the right
    grid = sns.clustermap(
        z, col_cluster=False, row_cluster=True, method=cluster_method, metric=cluster_metric,
        cmap=colors, vmin=-heatmap_z_limit, vmax=heatmap_z_limit, yticklabels=False, figsize=(6, 9),
        dendrogram_ratio=(0.18, 0.04), cbar_kws={"label": "Expression (Z-score)"},
        cbar_pos=(1.02, 0.72, 0.03, 0.15),
    )

    # clustermap makes several drawing areas; ax_heatmap is the heatmap
    # itself. Give its columns the short labels.
    grid.ax_heatmap.set_xticklabels(column_labels, rotation=90)
    grid.ax_heatmap.set_xlabel("")
    grid.ax_heatmap.set_ylabel("")
    grid.figure.suptitle(f"Fig. 2C: all {len(z)} significant genes, clustered", y=1.01)
    save_plot(grid.figure, output / fig2c_plot_file)

    # ---------------------------------------------------------------
    # STEP 3: How closely does each clone follow the shared pattern?
    # ---------------------------------------------------------------

    # The reference (WT) and comparison (KO) samples. The KO samples
    # are kept in the heatmap's column order.
    wt = samples_in_group(metadata, reference_group)
    ko_set = samples_in_group(metadata, comparison_group)
    ko = [sample for sample in order if sample in ko_set]

    # For each clone: its expression minus the WT average, for every
    # significant gene. Positive = higher than WT.
    wt_mean = vst.loc[significant.index, wt].mean(axis=1)

    change = {}
    for clone in ko:
        change[clone] = vst.loc[significant.index, clone] - wt_mean

    # A table with one column per clone.
    change = pd.DataFrame(change)

    # Measure 1: the share of significant genes where this clone
    # changes in the same direction as the overall result.
    # np.sign() turns every number into -1 (negative), 0 or +1
    # (positive). Two numbers go "the same direction" if their signs
    # are equal. The .mean() of a True/False column is the share of
    # Trues (0 to 1), so x 100 makes a percentage.
    overall_sign = np.sign(significant["log2FoldChange"])
    same_direction = {}
    for clone in ko:
        clone_sign = np.sign(change[clone])
        agrees = clone_sign == overall_sign
        same_direction[clone] = agrees.mean() * 100

    # Measure 2: how similar the clones' changes are to each other.
    # .corr() gives the correlation of every pair of clones.
    between = change.corr()

    # The clone that agrees least with the overall result: start with
    # the first clone, and switch whenever a later clone is lower.
    least_similar = ko[0]
    for clone in ko:
        if same_direction[clone] < same_direction[least_similar]:
            least_similar = clone

    # The other clones.
    others = []
    for clone in ko:
        if clone != least_similar:
            others.append(clone)

    # The most deeply sequenced KO sample: DESeq2 gives it the largest
    # size factor (saved in 01_sample_qc.tsv by fit_model.py).
    # .idxmax() = the row label (sample name) of the largest value.
    deepest = metadata.loc[ko, "size_factor"].idxmax()

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # e.g. "KO clone 16 100.0%, KO clone 21 99.8%, KO clone 20 96.6%"
    share_parts = []
    for clone in ko:
        share_parts.append(f"{label_of[clone]} {same_direction[clone]:.1f}%")
    share_text = ", ".join(share_parts)

    # Correlations of the least similar clone with each other clone.
    least_parts = []
    for clone in others:
        least_parts.append(f"{between.loc[least_similar, clone]:.2f}")

    # Correlations among the other clones (every pair once). The inner
    # loop starts after the outer one, so each pair is listed once.
    other_parts = []
    for first in range(len(others)):
        for second in range(first + 1, len(others)):
            a = others[first]
            b = others[second]
            other_parts.append(f"{label_of[a]} vs {label_of[b]} r = {between.loc[a, b]:.2f}")

    notes = []
    notes.append(
        f"{len(z)} genes split into two blocks, high in {reference_group} and high in the knockouts, "
        "as in the paper's heatmap."
    )
    notes.append("Share of these genes that change in the same direction as the overall result: " + share_text + ".")

    similarity_text = "Similarity of the clones' changes: "
    if len(other_parts) > 0:
        similarity_text = similarity_text + "; ".join(other_parts) + "; "
    similarity_text = (similarity_text + f"{label_of[least_similar]} vs the others r = "
                       + " and ".join(least_parts) + ".")
    notes.append(similarity_text)

    # Compare with the sample that looks most different in the paper's
    # heatmap (from paper_info.tsv).
    paper_different = info["fig2c_most_different_sample"]
    least_text = f"{label_of[least_similar]} follows the shared pattern least closely"
    if least_similar == paper_different:
        least_text = least_text + ", and the paper's heatmap shows the same: its column looks different from the other clones"
    else:
        least_text = (least_text + f", while in the paper's heatmap {label_of.get(paper_different, paper_different)} "
                      "looks the most different")
    if least_similar == deepest:
        size_factor = metadata.loc[deepest, "size_factor"]
        least_text = (least_text + f". It is also the most deeply sequenced knockout sample "
                      f"(DESeq2 size factor {size_factor:.2f}), which points to a technical difference")
    notes.append(least_text + ".")

    notes.append("Differs: the row order and exact colors, because the paper clustered its own gene list, "
                 f"with its own normalization ({methods['Statistics']['paper']}) and clustering settings, "
                 "which the methods don't state.")

    # "Reproduced" only if every clone follows the overall direction
    # for at least agree_direction_percent (80%) of the genes.
    all_clones_agree = all([same_direction[clone] >= agree_direction_percent for clone in ko])
    if all_clones_agree:
        status = "Reproduced"
    else:
        status = "Reproduced with differences"
    return notes, status


# ===================================================================
# Fig. 2D: top down and top up genes
# ===================================================================
#
# The paper shows two tables: its top down and top up genes (30 of
# each, the rows of paper_degs.tsv), ranked by adjusted p-value. We
# make the same two tables from our results, with the same number of
# genes, and shade the genes that are also in the paper's tables.

def add_top_rows(rows, table, direction, paper_set):
    """Add one row per gene in "table" to the list "rows"."""

    rank = 1
    for gene_id, gene in table.iterrows():
        rows.append({
            "direction": direction,
            "rank": rank,
            "gene_id": gene_id,
            "gene": gene["label"],
            "log2FoldChange": round(gene["log2FoldChange"], 3),
            "padj": f"{gene['padj']:.2e}",                               # e.g. 3.78e-60
            "in_paper_list": gene["gene_name"].upper() in paper_set,     # True / False
        })
        rank = rank + 1


def draw_top_table(ax, part, title, color):
    """Draw one table of top genes in the drawing area "ax"."""

    # The table's cells: one list per row, [gene, log2FC, padj].
    cells = []
    for row_number, row in part.iterrows():
        cells.append([row["gene"], f"{row['log2FoldChange']:.3f}", row["padj"]])

    # The header row. len(header) is the number of columns (3); we use
    # it below instead of typing 3.
    header = ["Gene", "Log2FC", "Padj"]

    # A table has no axes, so hide them.
    ax.axis("off")

    # ax.table() draws a table.
    #   colLabels -> the header row
    #   colWidths -> each column's share of the width (visual)
    table = ax.table(cellText=cells, colLabels=header, loc="upper center",
                     cellLoc="center", colWidths=[0.45, 0.25, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # Make the rows 25% taller, so the text has room.
    table.scale(1, 1.25)

    # Color the header row's text. table[row, column] is one cell;
    # row 0 is the header.
    for col in range(len(header)):
        table[0, col].set_text_props(color=color, weight="bold")

    # Shade the genes that are also in the paper's list. The first
    # gene is in table row 1 (row 0 is the header).
    table_row = 1
    for row_number, row in part.iterrows():
        if row["in_paper_list"]:
            for col in range(len(header)):
                table[table_row, col].set_facecolor(paper_highlight_color)
        table_row = table_row + 1

    ax.set_title(title, fontsize=12, weight="bold")


def plot_fig2d(results, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: The paper's top down and top up genes, as sets
    # ---------------------------------------------------------------

    paper = pd.read_csv(reference_folder / paper_genes_file, sep="\t")
    paper_down = set()
    paper_up = set()
    for row_number, row in paper.iterrows():
        if row["log2FoldChange"] < 0:
            paper_down.add(row["gene_name"].upper())
        elif row["log2FoldChange"] > 0:
            paper_up.add(row["gene_name"].upper())

    # How many genes to show in each of our tables: as many as the
    # paper's table has (30 down and 30 up for this paper).
    number_down = len(paper_down)
    number_up = len(paper_up)

    # ---------------------------------------------------------------
    # STEP 2: Our top down and top up genes
    # ---------------------------------------------------------------

    # Sort by adjusted p-value (then raw p-value), like the paper.
    ranked = results.sort_values(["padj", "pvalue"])

    is_down = ranked["direction"] == lower_label
    down = ranked[is_down].head(number_down)

    is_up = ranked["direction"] == higher_label
    up = ranked[is_up].head(number_up)

    # One table: the down genes, then the up genes.
    rows = []
    add_top_rows(rows, down, "down", paper_down)
    add_top_rows(rows, up, "up", paper_up)
    top = pd.DataFrame(rows)
    top.to_csv(output / fig2d_table_file, sep="\t", index=False)

    # ---------------------------------------------------------------
    # STEP 3: Draw the two tables side by side
    # ---------------------------------------------------------------

    # The figure gets taller with more rows: about 0.3 inches per row
    # plus 1 inch for the titles, but at least 4 inches (visual only).
    longest = max(number_down, number_up)
    height = max(4, 0.3 * longest + 1)
    fig, axes = plt.subplots(1, 2, figsize=(10, height))

    down_part = top[top["direction"] == "down"]
    draw_top_table(axes[0], down_part, "Top downregulated genes", paper_down_color)

    up_part = top[top["direction"] == "up"]
    draw_top_table(axes[1], up_part, "Top upregulated genes", paper_up_color)

    fig.suptitle("Fig. 2D (shaded = also in the paper's list)", y=0.995)
    save_plot(fig, output / fig2d_plot_file)

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # How many of our top genes are also in the paper's lists?
    # .sum() of a True/False column counts the Trues.
    down_overlap = down_part["in_paper_list"].sum()
    up_overlap = up_part["in_paper_list"].sum()

    # The paper's genes in our results, from compare_to_paper.py.
    number_found, number_same, correlation = read_gene_comparison(output)

    notes = []
    notes.append(
        f"Overlap with the paper's lists: {down_overlap} of our top {len(down_part)} down genes and "
        f"{up_overlap} of our top {len(up_part)} up genes are also in the paper's lists (shaded in the figure)."
    )

    # "Mostly different" only if fewer than half of our genes overlap.
    total_overlap = down_overlap + up_overlap
    total_shown = len(down_part) + len(up_part)
    if total_overlap < total_shown / 2:
        notes.append(
            f"Differs: the two lists are mostly different genes, even though {number_same} of the "
            f"{number_found} paper genes we found change the same way here (r = {correlation:.2f}). Both tables "
            "rank genes by adjusted p-value, and the two statistical methods order the strongest genes "
            f"differently: ours puts very strong, consistent changes such as {top_gene_names(results, report_examples)} "
            "at the top."
        )
    else:
        notes.append(f"Matches: most of our top genes are also in the paper's lists "
                     f"({total_overlap} of {total_shown}).")
    notes.append(
        f"The gene-by-gene comparison of the paper's {len(paper)} genes is in {gene_comparison_file} and "
        f"{fold_change_comparison_plot_file}."
    )
    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2E: ORA (the plot is made by plot_ora.py)
# ===================================================================
#
# Our 10_reactome_ora.png is already drawn like the paper's Fig. 2E,
# so here we only write the notes for the report.

def describe_fig2e(paper_numbers, methods, reference_folder, output):

    ora_file = output / ora_results_file
    if not ora_file.exists():
        return [f"{ora_results_file} not found, so Fig. 2E was skipped."], "skipped"

    # ---------------------------------------------------------------
    # STEP 1: Which of the paper's pathways are significant here?
    #         (worked out by compare_to_paper.py)
    # ---------------------------------------------------------------

    term_file = output / term_comparison_file
    found_text = ""
    if term_file.exists():
        terms = pd.read_csv(term_file, sep="\t")

        number_ora_terms = 0
        number_agree = 0
        not_significant = []

        for row_number, row in terms.iterrows():
            # Only the ORA rows (the file also holds GO terms).
            # compare_to_paper.py starts their name with the
            # collection label, e.g. "Reactome ORA (Fig. 2E)".
            if row["analysis"].startswith(ora_collection_label):
                number_ora_terms = number_ora_terms + 1
                if row["agrees"] == True:
                    number_agree = number_agree + 1
                else:
                    not_significant.append(row["paper_term"])

        found_text = f"{number_agree} of the paper's {number_ora_terms} pathways are significant here."
        if len(not_significant) > 0:
            found_text = found_text + " Not significant: " + ", ".join(not_significant) + "."

    # ---------------------------------------------------------------
    # STEP 2: Our top pathways, and how many are in the paper's list
    # ---------------------------------------------------------------

    # The table is sorted by FDR, so the first rows are the strongest.
    ora = pd.read_csv(ora_file, sep="\t")
    paper_ora = pd.read_csv(reference_folder / paper_ora_file, sep="\t")

    # The paper's pathways, simplified so they can be compared.
    paper_simple = set()
    for term in paper_ora["term"]:
        paper_simple.add(simplify_name(term))

    # Our top pathways: as many as the paper lists.
    our_top = ora.head(len(paper_ora))
    in_paper_list = 0
    for term in our_top["Term"]:
        if simplify_name(term) in paper_simple:
            in_paper_list = in_paper_list + 1

    # The names of our first few pathways, without their IDs.
    our_first = []
    for term in list(ora["Term"])[:report_examples]:
        our_first.append(remove_term_id(term))

    # Our strongest pathway's -log10(FDR), and the paper's (read off
    # its figure, in paper_summary.tsv).
    first_fdr = list(ora["fdr"])[0]
    best = -np.log10(first_fdr)
    paper_best = paper_numbers["ora_top_minus_log10_fdr"]

    # ---------------------------------------------------------------
    # STEP 3: Notes for the report
    # ---------------------------------------------------------------

    notes = []
    notes.append(f"Our version is {ora_plot_file}, drawn the same way: the top {ora_top} pathways by FDR "
                 "among the significant genes.")
    notes.append(
        f"{in_paper_list} of our top {len(our_top)} pathways are also in the paper's top {len(paper_ora)}. "
        f"Our top {report_examples}: " + ", ".join(our_first) + ". " + found_text
    )
    notes.append(
        f"Differs: our strongest pathway reaches -log10(FDR) = {best:.1f}, the paper's about {paper_best:g} "
        f"(read off its figure). Our gene list is a little different, our ORA "
        f"({methods['ORA']['this_project']}) tests against all genes we measured, and the paper's "
        f"({methods['ORA']['paper']}) used its own background. These change the p-values more than the order "
        "of the pathways."
    )
    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2F: GSEA with GO Molecular Function, as dot plots
# ===================================================================
#
# The paper shows two dot plots (its top terms down, and its top terms
# up):
#   x position  -> gene ratio
#   dot size    -> number of genes
#   dot color   -> -log10(adjusted p-value)

def draw_gsea_dots(fig, ax, gsea, sign, how_many, title, floor):
    """Draw one dot plot: the most significant terms with NES of this sign."""

    # Keep the terms whose NES has this sign (-1 = down, +1 = up).
    has_this_sign = np.sign(gsea["NES"]) == sign
    part = gsea[has_this_sign]

    # The "how_many" most significant of them (as many as the paper
    # shows in that direction).
    part = part.sort_values("FDR q-val", kind="stable")
    part = part.head(how_many)

    # matplotlib draws the first item at the BOTTOM, so reverse the
    # order to put the most significant term at the top.
    reversed_index = list(part.index)
    reversed_index.reverse()
    part = part.loc[reversed_index]

    # The term names: without the database ID, in capitals like the
    # paper, and broken onto two lines if longer than 38 characters
    # (visual only: it keeps the labels from squeezing the plot).
    names = []
    for term in part["Term"]:
        names.append(textwrap.fill(remove_term_id(term).upper(), 38))

    # One dot per term.
    #   s=...      -> dot size grows with the number of genes
    #   c=...      -> dot color from -log10(FDR)
    #   cmap       -> blue (weak) to red (strong)
    #   vmin, vmax -> the color scale runs from 0 to the measurable limit
    points = ax.scatter(part["gene_ratio"], names, s=part["count"] * 2 + 10, c=part["minus_log10_fdr"],
                        cmap="coolwarm", vmin=0, vmax=-np.log10(floor), edgecolors="none")

    ax.set_xlabel("Gene ratio (leading-edge genes / genes in term)")
    ax.set_title(title, fontsize=11, weight="bold")

    # The color bar beside the plot.
    colorbar = fig.colorbar(points, ax=ax)
    colorbar.set_label("-log10(FDR)")

    ax.margins(x=0.15, y=0.08)
    ax.tick_params(axis="y", labelsize=9)


def plot_fig2f(paper_numbers, methods, reference_group, comparison_group, reference_folder, output):

    # This panel needs the paper-style GSEA from run_paper_gsea.py.
    gsea_file = output / paper_gsea_results_file
    if not gsea_file.exists():
        notes = [f"Skipped: the paper-style GSEA ({paper_gsea_gene_sets}) didn't run, usually because the gene "
                 "sets couldn't be downloaded. Run the analysis with internet access to make this panel."]
        return notes, "Not run"

    gsea = pd.read_csv(gsea_file, sep="\t")
    paper_terms = pd.read_csv(reference_folder / paper_gsea_file, sep="\t")

    # ---------------------------------------------------------------
    # STEP 1: Gene ratio and gene count for every term
    # ---------------------------------------------------------------

    # GSEApy's "Tag %" column looks like "16/38": 16 "leading-edge"
    # genes (the genes that drive the enrichment) out of 38 genes in
    # the term. Split it at "/" into the two numbers.
    counts = []
    sizes = []
    for tag in gsea["Tag %"]:
        pieces = tag.split("/")
        counts.append(int(pieces[0]))
        sizes.append(int(pieces[1]))

    gene_ratios = []
    for i in range(len(counts)):
        gene_ratios.append(counts[i] / sizes[i])

    gsea["count"] = counts
    gsea["gene_ratio"] = gene_ratios

    # With 1,000 permutations (settings.py), the smallest FDR GSEApy
    # can report is about 1/1,000. An FDR of exactly 0 means "smaller
    # than that", so we show it at that limit (clip raises anything
    # smaller up to it).
    floor = 1 / permutations
    clipped_fdr = gsea["FDR q-val"].clip(lower=floor)
    gsea["minus_log10_fdr"] = -np.log10(clipped_fdr)

    # ---------------------------------------------------------------
    # STEP 2: Draw the two dot plots, one above the other
    # ---------------------------------------------------------------

    # How many terms to show in each direction: as many as the paper
    # shows (10 down and 10 up for this paper). Comparing the
    # direction column with a label gives True/False; .sum() counts
    # the Trues.
    paper_down = (paper_terms["direction"] == lower_label).sum()
    paper_up = (paper_terms["direction"] == higher_label).sum()

    comparison_text = f"{comparison_group} vs {reference_group}"

    fig, axes = plt.subplots(2, 1, figsize=(10, 11))
    draw_gsea_dots(fig, axes[0], gsea, -1, paper_down, comparison_text + " (downregulated)", floor)
    draw_gsea_dots(fig, axes[1], gsea, 1, paper_up, comparison_text + " (upregulated)", floor)

    fig.suptitle(f"Fig. 2F: GSEA with {paper_gsea_collection_label} (dot size = leading-edge genes)", y=1.0)

    # tight_layout() spaces the two plots so their labels don't overlap.
    fig.tight_layout()
    save_plot(fig, output / fig2f_plot_file)

    # ---------------------------------------------------------------
    # STEP 3: Compare with the paper's terms
    # ---------------------------------------------------------------

    # Simplified names, to match the paper's spelling.
    simple_names = []
    for term in gsea["Term"]:
        simple_names.append(simplify_name(term))
    gsea["simple"] = simple_names

    lookup = gsea.drop_duplicates("simple")
    lookup = lookup.set_index("simple")

    found = 0
    same = 0
    same_and_significant = 0
    missing = []

    for row_number, row in paper_terms.iterrows():
        key = simplify_name(row["term"])

        # Not in our results: remember it (in Title Case) and move on.
        if key not in lookup.index:
            missing.append(row["term"].title())
            continue

        found = found + 1
        ours = lookup.loc[key]

        if ours["NES"] > 0:
            our_direction = higher_label
        else:
            our_direction = lower_label

        if our_direction == row["direction"]:
            same = same + 1
            if ours["FDR q-val"] < fdr_cutoff:
                same_and_significant = same_and_significant + 1

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # Count the significant terms, down and up.
    number_significant = 0
    number_down = 0
    number_up = 0
    for row_number, row in gsea.iterrows():
        if row["FDR q-val"] < fdr_cutoff:
            number_significant = number_significant + 1
            if row["NES"] < 0:
                number_down = number_down + 1
            if row["NES"] > 0:
                number_up = number_up + 1

    # The biggest term tested.
    largest = max(sizes)

    # The top of the paper's color scale, read off its figure
    # (paper_summary.tsv).
    paper_scale = paper_numbers["gsea_top_minus_log10_padj"]

    notes = []
    notes.append(
        f"{number_significant} {paper_gsea_collection_label} terms are significant here ({number_down} down, "
        f"{number_up} up). The largest term tested has {largest:,} genes."
    )
    notes.append(
        f"Of the paper's {len(paper_terms)} terms, {found} are in our results; {same} of those change in the "
        f"same direction, and {same_and_significant} are also significant here."
    )
    if len(missing) > 0:
        notes.append("Not in our results: " + ", ".join(missing) + ". GO term names and gene sets differ between "
                     f"the GO version the paper used and the {paper_gsea_gene_sets} library we used.")
    notes.append(f"Differs: the paper's colour scale reaches -log10(padj) of about {paper_scale:g}; ours stops at "
                 f"{-np.log10(floor):.0f}. The paper's GSEA ({methods['GSEA']['paper']}) estimates very small "
                 f"p-values precisely; GSEApy with {permutations:,} permutations can't report an FDR much below "
                 f"1/{permutations:,}. The gene ratio and dot size are our closest equivalents to the paper's.")
    return notes, "Reproduced with differences"


# ###################################################################
# PART 3: Make every panel and write the report
# ###################################################################

def plot_paper_figures(output, reference_folder, reference_group, comparison_group):

    output = Path(output)
    reference_folder = Path(reference_folder)

    # ---------------------------------------------------------------
    # STEP 1: Read our results and the paper's facts
    # ---------------------------------------------------------------

    results = read_table(output / differential_expression_file)
    normalized = pd.read_csv(output / normalized_counts_file, sep="\t", index_col=0)
    vst = pd.read_csv(output / vst_expression_file, sep="\t", index_col=0)
    metadata = pd.read_csv(output / sample_qc_file, sep="\t", index_col=0)

    # The two groups must be in the sample table, or nothing below
    # can work. (This catches a typo in the group names.)
    groups = list(metadata[condition_column])
    if reference_group not in groups or comparison_group not in groups:
        sys.exit(f"The groups {reference_group} and {comparison_group} must both be in {sample_qc_file}.")

    info = read_paper_info(reference_folder)
    paper_numbers = read_paper_numbers(reference_folder)
    methods = read_methods(output, reference_folder)

    # ---------------------------------------------------------------
    # STEP 2: Make every panel
    # ---------------------------------------------------------------

    # One dictionary per panel:
    #   panel   -> its name in the paper
    #   what    -> what it shows
    #   figure  -> our file for it ("none" if it can't be made)
    #   notes   -> the sentences for the report
    #   status  -> e.g. "Reproduced with differences"
    panels = []

    notes, status = plot_fig1b(results, normalized, metadata, info, reference_group, reference_folder, output)
    panels.append({"panel": "Fig. 1B", "what": f"{info['fig1b_gene']} per clone ({info['fig1b_method']})",
                   "figure": fig1b_plot_file, "notes": notes, "status": status})

    notes, status = plot_fig2a(results, methods, reference_group, comparison_group, reference_folder, output)
    panels.append({"panel": "Fig. 2A", "what": "Volcano plot",
                   "figure": fig2a_plot_file, "notes": notes, "status": status})

    notes, status = plot_fig2b(results, paper_numbers, methods, reference_group, comparison_group, output)
    panels.append({"panel": "Fig. 2B", "what": "Number of significant genes",
                   "figure": fig2b_plot_file, "notes": notes, "status": status})

    notes, status = plot_fig2c(results, vst, metadata, info, methods, reference_group, comparison_group,
                               reference_folder, output)
    panels.append({"panel": "Fig. 2C", "what": "Clustered heatmap of significant genes",
                   "figure": fig2c_plot_file, "notes": notes, "status": status})

    notes, status = plot_fig2d(results, reference_folder, output)
    panels.append({"panel": "Fig. 2D", "what": "Top down and up genes",
                   "figure": fig2d_plot_file, "notes": notes, "status": status})

    notes, status = describe_fig2e(paper_numbers, methods, reference_folder, output)
    panels.append({"panel": "Fig. 2E", "what": ora_collection_label + " ORA",
                   "figure": ora_plot_file, "notes": notes, "status": status})

    notes, status = plot_fig2f(paper_numbers, methods, reference_group, comparison_group, reference_folder, output)
    panels.append({"panel": "Fig. 2F", "what": "GSEA, " + paper_gsea_collection_label,
                   "figure": fig2f_plot_file, "notes": notes, "status": status})

    # The panels that RNA-seq can't redraw, with the reason for each,
    # from paper_not_reproducible.tsv.
    not_reproducible_file = reference_folder / paper_not_reproducible_file
    if not_reproducible_file.exists():
        table = pd.read_csv(not_reproducible_file, sep="\t", dtype=str)
        for row_number, row in table.iterrows():
            panels.append({"panel": row["panel"], "what": row["what"], "figure": "none",
                           "notes": [row["reason"]], "status": "Not reproducible from RNA-seq"})

    # Put the panels in the paper's order: "Fig. 1B" < "Fig. 1C" <
    # "Fig. 2A" < ... < "Fig. 3A–E" when sorted as text.
    # sorted(list, key=...) sorts by the value the key function gives
    # for each item; here, each panel's name.
    def panel_name(panel):
        return panel["panel"]

    panels = sorted(panels, key=panel_name)

    # ---------------------------------------------------------------
    # STEP 3: Write the report
    # ---------------------------------------------------------------

    # The title and a summary table with one row per panel.
    report = []
    report.append("# Reproducing the paper's figures")
    report.append("")
    report.append(paper_citation(info))
    report.append("")
    report.append("Every panel of the paper's figures, what this project made for it, and why any differences "
                  "are expected.")
    report.append("")
    report.append("| Panel | What it shows | Our figure | Status |")
    report.append("|---|---|---|---|")

    for panel in panels:
        # Show the file name in `code style`, or a dash if there's none.
        if panel["figure"] != "none":
            shown = f"`{panel['figure']}`"
        else:
            shown = "—"
        report.append(f"| {panel['panel']} | {panel['what']} | {shown} | {panel['status']} |")
    report.append("")

    # Then one section per panel, with its notes as bullet points.
    for panel in panels:
        report.append(f"## {panel['panel']}: {panel['what']}")
        report.append("")
        for note in panel["notes"]:
            report.append(f"- {note}")
        report.append("")

    # encoding="utf-8" saves symbols such as "≥" and "—" correctly on
    # every computer.
    with open(output / figures_report_file, "w", encoding="utf-8") as file:
        for line in report:
            file.write(line + "\n")

    # A short summary in Terminal.
    print("Paper figures:")
    for panel in panels:
        print(f"  {panel['panel']}: {panel['status']}")


# ###################################################################
# PART 4: Run this file on its own
# ###################################################################

# This part runs only when you start this file from Terminal:
#     python3 plot_paper_figures.py akap11_go_bp_results WT KO
# It does NOT run when run_analysis.py imports this file.
if __name__ == "__main__":

    # Save plots to files instead of opening windows, and use the same
    # plot style as run_analysis.py.
    import matplotlib
    matplotlib.use("Agg")
    sns.set_theme(style="whitegrid")

    from settings import paper_reference_folder

    # We need exactly 4 words: the script name, the results folder,
    # the reference group and the comparison group.
    if len(sys.argv) != 4:
        sys.exit("Usage: python3 plot_paper_figures.py results_folder reference_group comparison_group")

    plot_paper_figures(sys.argv[1], paper_reference_folder, sys.argv[2], sys.argv[3])
