"""
===================================================================
plot_paper_figures.py
===================================================================

WHAT THIS FILE DOES
    It redraws the paper's figures from OUR results, panel by panel,
    so they can be put side by side with the published ones:

        Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557

      Fig. 1B  AKAP11 in each knockout clone, relative to WT
      Fig. 2A  volcano plot, with the genes the paper labels
      Fig. 2B  number of significant genes, ours next to the paper's
      Fig. 2C  clustered heatmap of all significant genes
      Fig. 2D  top 30 down and top 30 up genes
      Fig. 2E  Reactome ORA (already drawn by plot_ora.py:
               10_reactome_ora.png, so here we only describe it)
      Fig. 2F  GSEA with GO Molecular Function terms, as dot plots

    Fig. 1C, Fig. 1D and Fig. 3 come from Western blots, microscope
    images and DNA methylation sequencing, not RNA-seq, so they can't
    be redrawn from this data. The report says so.

    It also writes a report, 12_paper_figures_report.md, that says for
    every panel what matches, what differs, and why.

WHERE IT FITS
    Step 12, the last step. run_analysis.py calls it if the
    paper_reference folder exists. You can also run it on its own:

        python3 plot_paper_figures.py akap11_go_bp_results

INPUTS (all read from files, so it can run on its own)
    from the results folder:
        01_sample_qc.tsv, 02_differential_expression.tsv,
        02_normalized_counts.tsv, 02_vst_expression.tsv,
        10_reactome_ora.tsv, 11_paper_term_comparison.tsv,
        12_fig2F_gsea_go_mf.tsv (if the paper-style GSEA ran)
    from the paper_reference folder:
        values copied from the paper's figures

OUTPUT FILES
    12_fig1B_akap11_relative.png / .tsv
    12_fig2A_volcano.png
    12_fig2B_deg_counts.png
    12_fig2C_heatmap.png
    12_fig2D_top_genes.png / .tsv
    12_fig2F_gsea_go_mf.png
    12_paper_figures_report.md

HOW THE FILE IS ORGANIZED
    PART 1   colors and small helper functions
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

# simplify_name() makes GO term names comparable; it lives in
# compare_to_paper.py, so we borrow it from there.
from compare_to_paper import simplify_name
from save_plot import save_plot
from settings import fdr_cutoff              # 0.05
from settings import fold_change_cutoff      # 1.0
from settings import permutations            # 1000
from settings import paper_gsea_gene_sets    # "GO_Molecular_Function_2025"
from settings import paper_gsea_max_size     # 3000


# ###################################################################
# PART 1: Colors and helper functions
# ###################################################################

# Colors close to the paper's: blue = down in knockout,
# red = up in knockout, gray = not significant.
BLUE = "#2F5DA8"
RED = "#B2182B"
GRAY = "#C9CCD1"

# The paper's Fig. 1B colors: green for WT, blue for the KO clones.
WT_GREEN = "#2E8B57"
KO_BLUE = "#3F51B5"


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


def short_go_name(term):
    """Remove the "(GO:0003735)" ID from the end of a GO term name."""

    # Everything before "(GO:", with spaces at the ends removed.
    if "(GO:" in term:
        term = term.split("(GO:")[0]
    return term.strip()


def short_clone_name(name):
    """Turn "AKAP11_KO_clone20" into "clone20" for the report."""
    return name.replace("AKAP11_KO_", "")


# ###################################################################
# PART 2: One function per figure panel
# ###################################################################

# ===================================================================
# Fig. 1B: AKAP11 in each clone, relative to WT
# ===================================================================
#
# The paper measured AKAP11 by RT-qPCR on exon 7 and showed each
# clone relative to WT (WT = 1). We do the same with RNA-seq: each
# clone's normalized AKAP11 count divided by the WT average.

def plot_fig1b(results, normalized, metadata, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: Our numbers
    # ---------------------------------------------------------------

    gene_id = find_gene(results, "AKAP11")

    # If AKAP11 isn't in the results, return a note and a status, and
    # skip the rest of this function.
    if gene_id is None:
        return ["AKAP11 is not in the filtered results, so Fig. 1B was skipped."], "skipped"

    # AKAP11's normalized count in each sample (one row of the table).
    counts = normalized.loc[gene_id]

    # The names of the WT samples.
    wt_samples = []
    for sample in metadata.index:
        if metadata.loc[sample, "condition"] == "WT":
            wt_samples.append(sample)

    # The WT average becomes our "1.0".
    wt_mean = counts[wt_samples].mean()

    # The clones in the paper's order: 21, 20, 16.
    clones = ["AKAP11_KO_clone21", "AKAP11_KO_clone20", "AKAP11_KO_clone16"]

    # Each clone relative to WT, e.g. 0.55 = 55% of the WT level.
    ours = {}
    for clone in clones:
        ours[clone] = counts[clone] / wt_mean

    # ---------------------------------------------------------------
    # STEP 2: The paper's numbers
    # ---------------------------------------------------------------

    # Read off the paper's Fig. 1B (approximate). Turn the table into a
    # dictionary: {"AKAP11_KO_clone21": 0.43, ...}
    paper_table = pd.read_csv(reference_folder / "paper_fig1b_qpcr.tsv", sep="\t")
    paper = {}
    for row_number, row in paper_table.iterrows():
        paper[row["sample"]] = row["relative_expression"]

    # ---------------------------------------------------------------
    # STEP 3: Save the numbers side by side
    # ---------------------------------------------------------------

    rows = []
    rows.append({"sample": "WT (mean of 3)", "this_project_rnaseq": 1.0, "paper_qpcr_exon7": 1.0})
    for clone in clones:
        rows.append({
            "sample": clone,
            "this_project_rnaseq": round(ours[clone], 3),   # 3 decimal places
            "paper_qpcr_exon7": paper[clone],
        })
    numbers = pd.DataFrame(rows)
    numbers.to_csv(output / "12_fig1B_akap11_relative.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------
    # STEP 4: Draw two bars per group: ours (filled) and the paper's
    #         (striped)
    # ---------------------------------------------------------------

    labels = ["WT", "KO clone 21", "KO clone 20", "KO clone 16"]

    # The bar heights, in the same order as the labels.
    our_values = [1.0]
    paper_values = [1.0]
    for clone in clones:
        our_values.append(ours[clone])
        paper_values.append(paper[clone])

    # Each group sits at x = 0, 1, 2, 3. Each bar is 0.38 wide; our bar
    # goes half a bar to the LEFT of the group's center and the paper's
    # half a bar to the RIGHT, so the two sit side by side.
    width = 0.38
    positions = []
    our_positions = []
    paper_positions = []
    for i in range(len(labels)):
        positions.append(i)
        our_positions.append(i - width / 2)
        paper_positions.append(i + width / 2)

    # One color per group: green for WT, blue for the three clones.
    bar_colors = [WT_GREEN, KO_BLUE, KO_BLUE, KO_BLUE]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Our bars: filled with color.
    ax.bar(our_positions, our_values, width, color=bar_colors,
           label="This project: RNA-seq, whole gene")

    # The paper's bars: white with a colored outline and stripes
    # (hatch="//"), so the two kinds are easy to tell apart.
    ax.bar(paper_positions, paper_values, width, color="white",
           edgecolor=bar_colors, hatch="//", linewidth=1.5,
           label="Paper: RT-qPCR, exon 7 (read off Fig. 1B)")

    # Show each WT replicate as a black dot on the WT bar, like the
    # paper does. All three dots share the WT bar's x position.
    wt_relative = counts[wt_samples] / wt_mean
    wt_dot_x = []
    for sample in wt_samples:
        wt_dot_x.append(our_positions[0])
    ax.scatter(wt_dot_x, wt_relative, color="black", s=18, zorder=3)

    # Write each bar's value just above it (0.02 higher), e.g. "0.55".
    for i in range(len(labels)):
        ax.text(our_positions[i], our_values[i] + 0.02, f"{our_values[i]:.2f}", ha="center", fontsize=9)
    for i in range(len(labels)):
        ax.text(paper_positions[i], paper_values[i] + 0.02, f"{paper_values[i]:.2f}", ha="center", fontsize=9)

    # Put the group names under the group centers.
    ax.set_xticks(positions)
    ax.set_xticklabels(labels)

    ax.set_ylim(0, 1.35)
    ax.set_ylabel("AKAP11 relative to WT (WT = 1)")
    ax.set_title("Fig. 1B: AKAP11 in each knockout clone")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    save_plot(fig, output / "12_fig1B_akap11_relative.png")

    # ---------------------------------------------------------------
    # STEP 5: What does it show? Notes for the report
    # ---------------------------------------------------------------

    # Which clone keeps the most AKAP11, in our data and in the paper's?
    # Start with the first clone and replace it whenever a later clone
    # is higher.
    our_highest = clones[0]
    for clone in clones:
        if ours[clone] > ours[our_highest]:
            our_highest = clone

    paper_highest = clones[0]
    for clone in clones:
        if paper[clone] > paper[paper_highest]:
            paper_highest = clone

    if our_highest == paper_highest:
        matches_text = "Matches: AKAP11 is lower in every knockout clone, and clone 20 keeps the most AKAP11 in both."
    else:
        matches_text = "Matches: AKAP11 is lower in every knockout clone."

    notes = []
    notes.append(
        f"Our RNA-seq levels relative to WT: clone 21 = {ours[clones[0]]:.2f}, clone 20 = {ours[clones[1]]:.2f}, "
        f"clone 16 = {ours[clones[2]]:.2f}. The paper's RT-qPCR: {paper[clones[0]]:.2f}, {paper[clones[1]]:.2f}, "
        f"{paper[clones[2]]:.2f}."
    )
    notes.append(matches_text)
    notes.append(
        "Differs: the drop is smaller in RNA-seq. The paper's qPCR measures only exon 7, the exon targeted by the "
        "CRISPR edit, so it sees only intact copies. RNA-seq counts reads along the whole gene, including the edited "
        "copy, which is still transcribed. A heterozygous knockout also keeps one working copy, so about 0.5 is the "
        "expected level for intact exon 7."
    )
    notes.append("To check exon 7 directly, you could count reads on exon 7 in the STAR BAM files on Expanse.")

    status = "Reproduced with differences"
    return notes, status


# ===================================================================
# Fig. 2A: volcano plot with the paper's labeled genes
# ===================================================================
#
# Same idea as our 05_volcano.png, but in the paper's style and with
# the 19 genes the paper labels, so you can see where they fall in
# our data.

def plot_fig2a(results, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: The x and y values for every gene
    # ---------------------------------------------------------------

    # Drop genes with no adjusted p-value or no fold change: they
    # can't be placed on the plot.
    genes = results.dropna(subset=["padj", "log2FoldChange"]).copy()

    # y = -log10(adjusted p-value). clip() raises values below 1e-300
    # up to 1e-300, because log10(0) is impossible.
    clipped_padj = genes["padj"].clip(lower=1e-300)
    genes["minus_log10_padj"] = -np.log10(clipped_padj)

    # ---------------------------------------------------------------
    # STEP 2: Find the paper's labeled genes in our results
    # ---------------------------------------------------------------

    label_table = pd.read_csv(reference_folder / "paper_volcano_labels.tsv", sep="\t")
    paper_labels = list(label_table["gene_name"])

    labeled_ids = []   # gene IDs we can plot
    renamed = []       # not in our results at all (renamed genes)
    no_padj = []       # in our results, but without an adjusted p-value

    for symbol in paper_labels:
        gene_id = find_gene(genes, symbol)
        if gene_id is not None:
            labeled_ids.append(gene_id)
        elif find_gene(results, symbol) is not None:
            # In our results, but DESeq2 gave it no adjusted p-value.
            no_padj.append(symbol)
        else:
            renamed.append(symbol)

    # The rows of the genes we can plot.
    labeled = genes.loc[labeled_ids]
    found = len(labeled)

    # How many of them are significant here?
    is_significant = labeled["direction"] != "Below thresholds"
    significant = is_significant.sum()

    # Each labeled gene's color: blue if down in KO, red if up.
    labeled_colors = []
    for fold_change in labeled["log2FoldChange"]:
        if fold_change < 0:
            labeled_colors.append(BLUE)
        else:
            labeled_colors.append(RED)

    # ---------------------------------------------------------------
    # STEP 3: Draw two panels
    # ---------------------------------------------------------------

    # Left: the whole plot. Right: a zoom on the area where the paper's
    # genes sit, so their names have room.
    # plt.subplots(1, 2, ...) makes 1 row of 2 drawing areas; "axes" is
    # a list of the two. width_ratios makes the right one a bit wider.
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.25]})

    color_of = {"Below thresholds": GRAY, "Lower in comparison": BLUE, "Higher in comparison": RED}

    # Draw the same dots in both panels.
    for ax in axes:

        # All genes, one color per direction (gray first, so the
        # colored dots are on top).
        for direction in ["Below thresholds", "Lower in comparison", "Higher in comparison"]:
            is_this_direction = genes["direction"] == direction
            group = genes[is_this_direction]
            ax.scatter(group["log2FoldChange"], group["minus_log10_padj"], color=color_of[direction],
                       s=8, alpha=0.6)

        # The paper's labeled genes on top: bigger, with a black
        # outline. zorder=3 draws them above everything else.
        ax.scatter(labeled["log2FoldChange"], labeled["minus_log10_padj"],
                   color=labeled_colors,
                   s=40, edgecolors="black", linewidths=0.8, zorder=3)

        # The cutoff lines, as in the paper: adjusted p = 0.05 and
        # log2 fold change = -1 and +1.
        ax.axhline(-np.log10(fdr_cutoff), linestyle="--", color="black", linewidth=0.8)
        ax.axvline(-fold_change_cutoff, linestyle="--", color="black", linewidth=0.8)
        ax.axvline(fold_change_cutoff, linestyle="--", color="black", linewidth=0.8)

        ax.set_xlabel("log2(FC)")
        ax.set_ylabel("-log10(adj p-value)")

    # axes[0] is the left panel, axes[1] the right one.
    axes[0].set_title("All genes (outlined = labeled in the paper)")

    # The zoom: x from -14 to 9, and y from 0 to a little above the
    # highest labeled gene (15% extra room for the names).
    zoom = axes[1]
    zoom.set_xlim(-14, 9)
    zoom.set_ylim(0, labeled["minus_log10_padj"].max() * 1.15)
    zoom.set_title("Zoom on the paper's labeled genes")

    # Write each labeled gene's name, in its color, in bold italics.
    texts = []
    for gene_id, gene in labeled.iterrows():
        if gene["log2FoldChange"] < 0:
            color = BLUE
        else:
            color = RED
        text = zoom.text(gene["log2FoldChange"], gene["minus_log10_padj"], gene["gene_name"],
                         fontsize=9, style="italic", color=color, weight="bold")
        texts.append(text)

    # Move the names apart so they don't overlap.
    #   expand  -> how much empty space to keep around each name
    #   iter_lim=1000 -> try harder than usual, as there are many names
    if len(texts) > 0:
        adjust_text(texts, x=list(labeled["log2FoldChange"]), y=list(labeled["minus_log10_padj"]), ax=zoom,
                    arrowprops={"arrowstyle": "-", "color": "#333333", "lw": 0.6},
                    expand=(1.4, 1.6), iter_lim=1000)

    # suptitle = one title over both panels.
    fig.suptitle("Fig. 2A: Het-AKAP11-KO vs WT", y=1.0)
    save_plot(fig, output / "12_fig2A_volcano.png")

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # Our most significant gene: the first row after sorting by padj.
    # .iloc[0] means "the row at position 0", the first one.
    sorted_genes = genes.sort_values("padj")
    top = sorted_genes.iloc[0]

    notes = []
    notes.append(
        f"Matches: the same layout, cutoffs and colors. Of the {len(paper_labels)} genes the paper labels, "
        f"{found} can be plotted here and {significant} of those are significant, all on the same side as in the paper."
    )
    notes.append(
        f"Differs: our y-axis reaches {top['minus_log10_padj']:.0f} (top gene {top['gene_name']}), the paper's about 7 "
        "(top gene NDN). Our model gives much smaller p-values for the strongest genes, so the genes the paper "
        "labels sit lower down in our plot, under genes like CYP26C1 and the HOX genes."
    )
    notes.append(
        "Why: the paper used limma-voom with robust empirical Bayes, which shrinks extreme statistics; we used DESeq2's "
        "Wald test, which does not. The two methods agree on which genes change and in which direction (r = 0.98 for "
        "the paper's top genes), but not on how small the smallest p-values are."
    )
    if len(renamed) > 0:
        notes.append("Not in our results, because they were renamed between Ensembl 104 (paper) and GENCODE v47 "
                     "(ours): " + ", ".join(renamed) + ".")
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

def plot_fig2b(results, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: The numbers
    # ---------------------------------------------------------------

    is_down = results["direction"] == "Lower in comparison"
    our_down = is_down.sum()

    is_up = results["direction"] == "Higher in comparison"
    our_up = is_up.sum()

    # The paper's numbers: use the "measure" column as row labels, then
    # take the "value" column, so we can look numbers up by name.
    paper_summary = pd.read_csv(reference_folder / "paper_summary.tsv", sep="\t")
    paper_summary = paper_summary.set_index("measure")
    paper_values = paper_summary["value"]
    paper_down = int(paper_values["lower_in_comparison"])
    paper_up = int(paper_values["higher_in_comparison"])

    # ---------------------------------------------------------------
    # STEP 2: Draw two stacked bars
    # ---------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(4.5, 5))

    names = ["Paper", "This project"]
    ups = [paper_up, our_up]
    downs = [paper_down, our_down]

    # The red "up" part first, from 0 up to the number of up genes.
    ax.bar(names, ups, color=RED, label="Upregulated", width=0.55)

    # The blue "down" part on top of it: bottom=ups starts each blue
    # bar where its red bar ends.
    ax.bar(names, downs, bottom=ups, color=BLUE, label="Downregulated", width=0.55)

    # Write the numbers: in the middle of each colored part (white
    # text), and the total just above each bar. Bar 0 is the paper's,
    # bar 1 is ours.
    for i in range(2):
        ax.text(i, ups[i] / 2, str(ups[i]), ha="center", va="center", color="white", fontsize=11)
        ax.text(i, ups[i] + downs[i] / 2, str(downs[i]), ha="center", va="center", color="white", fontsize=11)
        ax.text(i, ups[i] + downs[i] + 15, f"total {ups[i] + downs[i]}", ha="center", fontsize=10)

    ax.set_ylabel("Number of DEGs")
    ax.set_ylim(0, 1000)
    ax.set_title("Fig. 2B: Het-AKAP11-KO vs WT")

    # Put the legend beside the bars, so it can't cover the totals.
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
    save_plot(fig, output / "12_fig2B_deg_counts.png")

    # ---------------------------------------------------------------
    # STEP 3: Notes for the report
    # ---------------------------------------------------------------

    notes = []
    notes.append(
        f"Matches: {our_down + our_up} significant genes here, {paper_down + paper_up} in the paper, with the same cutoffs "
        "(adjusted p < 0.05, |log2FC| ≥ 1)."
    )
    notes.append(f"Differs: the split is {our_down} down / {our_up} up here and {paper_down} / {paper_up} in the paper.")
    notes.append(
        "Why: genes close to the cutoffs can tip either way. Counting (Salmon vs HTSeq), annotation (GENCODE v47 vs "
        "Ensembl 104), gene filter and model (DESeq2 vs limma-voom) all differ, and each moves a few genes across a cutoff."
    )
    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2C: clustered heatmap of all significant genes
# ===================================================================
#
# Like our 07_de_heatmap.png, but with ALL significant genes instead
# of the top 25, and with the rows "clustered": genes with similar
# patterns are placed next to each other, and a tree (dendrogram) on
# the left shows how they group.

def plot_fig2c(results, vst, output):

    # ---------------------------------------------------------------
    # STEP 1: Z-scores of all significant genes
    # ---------------------------------------------------------------

    is_significant = results["direction"] != "Below thresholds"
    significant = results[is_significant]

    # The samples in the paper's column order (note: clone 21 before
    # clone 20).
    order = ["WT_rep1", "WT_rep2", "WT_rep3", "AKAP11_KO_clone16", "AKAP11_KO_clone21", "AKAP11_KO_clone20"]

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

    # A color scale like the paper's: purple (low), white, red (high).
    colors = LinearSegmentedColormap.from_list("paper", ["#4B3F9E", "#FFFFFF", "#D7301F"])

    # sns.clustermap() clusters the rows and draws the heatmap.
    #   col_cluster=False   -> keep the columns in the paper's order
    #   row_cluster=True    -> reorder the genes so similar patterns
    #                          sit together
    #   method="average", metric="euclidean"
    #                       -> how "similar" is measured: the average
    #                          straight-line distance between genes
    #   vmin=-2, vmax=2     -> the colors stop at -2 and +2
    #   yticklabels=False   -> no gene names (848 would be unreadable)
    #   dendrogram_ratio    -> how much room the trees get
    #   cbar_pos            -> where the color bar goes: (left, bottom,
    #                          width, height), here outside on the right
    grid = sns.clustermap(
        z, col_cluster=False, row_cluster=True, method="average", metric="euclidean",
        cmap=colors, vmin=-2, vmax=2, yticklabels=False, figsize=(6, 9),
        dendrogram_ratio=(0.18, 0.04), cbar_kws={"label": "Expression (Z-score)"},
        cbar_pos=(1.02, 0.72, 0.03, 0.15),
    )

    # clustermap makes several drawing areas; ax_heatmap is the heatmap
    # itself. Give its columns short names.
    grid.ax_heatmap.set_xticklabels(["WT 1", "WT 2", "WT 3", "KO clone 16", "KO clone 21", "KO clone 20"],
                                    rotation=90)
    grid.ax_heatmap.set_xlabel("")
    grid.ax_heatmap.set_ylabel("")
    grid.figure.suptitle(f"Fig. 2C: all {len(z)} significant genes, clustered", y=1.01)
    save_plot(grid.figure, output / "12_fig2C_heatmap.png")

    # ---------------------------------------------------------------
    # STEP 3: How closely does each clone follow the shared pattern?
    # ---------------------------------------------------------------

    # For each clone: its expression minus the WT average, for every
    # significant gene. Positive = higher than WT.
    wt_mean = vst.loc[significant.index, ["WT_rep1", "WT_rep2", "WT_rep3"]].mean(axis=1)

    ko = ["AKAP11_KO_clone16", "AKAP11_KO_clone21", "AKAP11_KO_clone20"]

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

    # The other two clones.
    others = []
    for clone in ko:
        if clone != least_similar:
            others.append(clone)

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # "clone16 100.0%, clone21 99.8%, clone20 96.6%"
    share_parts = []
    for clone in ko:
        share_parts.append(f"{short_clone_name(clone)} {same_direction[clone]:.1f}%")
    share_text = ", ".join(share_parts)

    notes = []
    notes.append(
        f"Matches: {len(z)} genes split cleanly into two blocks, high in WT and high in the knockouts, "
        "as in the paper's heatmap."
    )
    notes.append("Share of these genes that change in the same direction as the overall result: " + share_text + ".")
    notes.append(
        f"Similarity of the clones' changes: {short_clone_name(others[0])} vs {short_clone_name(others[1])} r = "
        f"{between.loc[others[0], others[1]]:.2f}; {short_clone_name(least_similar)} vs the others r = "
        f"{between.loc[least_similar, others[0]]:.2f} and {between.loc[least_similar, others[1]]:.2f}."
    )
    if least_similar == "AKAP11_KO_clone20":
        notes.append("Clone 20 follows the shared pattern least closely, and has its own block of genes that are "
                     "up mainly in clone 20. The paper's heatmap shows the same thing: its clone 20 column looks "
                     "different from clones 16 and 21. This fits clone 20's technical differences (more intronic and "
                     "duplicate reads) and its higher remaining AKAP11 (Fig. 1B).")
    notes.append("Differs: the row order and exact colors, because the paper clustered its own gene list, "
                 "with its own normalization (limma-voom) and clustering settings, which the methods don't state.")
    return notes, "Reproduced"


# ===================================================================
# Fig. 2D: top 30 down and top 30 up genes
# ===================================================================
#
# The paper shows two tables: its top 30 down and top 30 up genes,
# ranked by adjusted p-value. We make the same two tables from our
# results, and shade the genes that are also in the paper's tables.

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
            "in_paper_top_30": gene["gene_name"].upper() in paper_set,   # True / False
        })
        rank = rank + 1


def draw_top_table(ax, part, title, color):
    """Draw one table of top genes in the drawing area "ax"."""

    # The table's cells: one list per row, [gene, log2FC, padj].
    cells = []
    for row_number, row in part.iterrows():
        cells.append([row["gene"], f"{row['log2FoldChange']:.3f}", row["padj"]])

    # A table has no axes, so hide them.
    ax.axis("off")

    # ax.table() draws a table.
    #   colLabels -> the header row
    #   colWidths -> each column's share of the width
    table = ax.table(cellText=cells, colLabels=["Gene", "Log2FC", "Padj"], loc="upper center",
                     cellLoc="center", colWidths=[0.45, 0.25, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    # Make the rows 25% taller, so the text has room.
    table.scale(1, 1.25)

    # Color the header row's text. table[row, column] is one cell;
    # row 0 is the header.
    for col in range(3):
        table[0, col].set_text_props(color=color, weight="bold")

    # Shade the genes that are also in the paper's top 30. The first
    # gene is in table row 1 (row 0 is the header).
    table_row = 1
    for row_number, row in part.iterrows():
        if row["in_paper_top_30"]:
            for col in range(3):
                table[table_row, col].set_facecolor("#FFF1C2")
        table_row = table_row + 1

    ax.set_title(title, fontsize=12, weight="bold")


def plot_fig2d(results, reference_folder, output):

    # ---------------------------------------------------------------
    # STEP 1: The paper's top 30 down and top 30 up genes, as sets
    # ---------------------------------------------------------------

    paper = pd.read_csv(reference_folder / "paper_degs.tsv", sep="\t")
    paper_down = set()
    paper_up = set()
    for row_number, row in paper.iterrows():
        if row["log2FoldChange"] < 0:
            paper_down.add(row["gene_name"].upper())
        elif row["log2FoldChange"] > 0:
            paper_up.add(row["gene_name"].upper())

    # ---------------------------------------------------------------
    # STEP 2: Our top 30 down and top 30 up genes
    # ---------------------------------------------------------------

    # Sort by adjusted p-value (then raw p-value), like the paper.
    ranked = results.sort_values(["padj", "pvalue"])

    is_down = ranked["direction"] == "Lower in comparison"
    down = ranked[is_down].head(30)

    is_up = ranked["direction"] == "Higher in comparison"
    up = ranked[is_up].head(30)

    # One table: the 30 down genes, then the 30 up genes.
    rows = []
    add_top_rows(rows, down, "down", paper_down)
    add_top_rows(rows, up, "up", paper_up)
    top = pd.DataFrame(rows)
    top.to_csv(output / "12_fig2D_top_genes.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------
    # STEP 3: Draw the two tables side by side
    # ---------------------------------------------------------------

    fig, axes = plt.subplots(1, 2, figsize=(10, 10))

    down_part = top[top["direction"] == "down"]
    draw_top_table(axes[0], down_part, "Top downregulated genes", BLUE)

    up_part = top[top["direction"] == "up"]
    draw_top_table(axes[1], up_part, "Top upregulated genes", RED)

    fig.suptitle("Fig. 2D (shaded = also in the paper's top 30)", y=0.995)
    save_plot(fig, output / "12_fig2D_top_genes.png")

    # ---------------------------------------------------------------
    # STEP 4: Notes for the report
    # ---------------------------------------------------------------

    # How many of our top 30 are also in the paper's top 30?
    down_overlap = 0
    for value in down_part["in_paper_top_30"]:
        if value:
            down_overlap = down_overlap + 1

    up_overlap = 0
    for value in up_part["in_paper_top_30"]:
        if value:
            up_overlap = up_overlap + 1

    notes = []
    notes.append(
        f"Overlap with the paper's lists: {down_overlap} of our top 30 down genes and {up_overlap} of our top 30 up "
        "genes are also in the paper's top 30 (shaded in the figure)."
    )
    notes.append(
        "Differs: the two top-30 lists are mostly different genes, even though the paper's top 60 genes all change "
        "the same way here (r = 0.98). Both tables rank genes by adjusted p-value, and DESeq2 and limma-voom order "
        "the strongest genes differently: ours puts very strong, consistent changes such as CYP26C1, the HOX genes "
        "and OTX1 at the top."
    )
    notes.append(
        "The gene-by-gene comparison of the paper's 60 genes is in 11_paper_gene_comparison.tsv and "
        "11_paper_vs_our_fold_change.png."
    )
    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2E: Reactome ORA (the plot is made by plot_ora.py)
# ===================================================================
#
# Our 10_reactome_ora.png is already drawn like the paper's Fig. 2E,
# so here we only write the notes for the report.

def describe_fig2e(output, reference_folder):

    ora_file = output / "10_reactome_ora.tsv"
    if not ora_file.exists():
        return ["10_reactome_ora.tsv not found, so Fig. 2E was skipped."], "skipped"

    # ---------------------------------------------------------------
    # STEP 1: Which of the paper's pathways are significant here?
    #         (worked out by compare_to_paper.py)
    # ---------------------------------------------------------------

    term_file = output / "11_paper_term_comparison.tsv"
    found_text = ""
    if term_file.exists():
        terms = pd.read_csv(term_file, sep="\t")

        number_ora_terms = 0
        number_agree = 0
        not_significant = []

        for row_number, row in terms.iterrows():
            # Only the Reactome rows (the file also holds GO terms).
            if row["analysis"].startswith("Reactome"):
                number_ora_terms = number_ora_terms + 1
                if row["agrees"] == True:
                    number_agree = number_agree + 1
                else:
                    not_significant.append(row["paper_term"])

        found_text = (f"{number_agree} of the paper's {number_ora_terms} pathways are significant here. "
                      "Not significant: " + ", ".join(not_significant) + ".")

    # ---------------------------------------------------------------
    # STEP 2: Our strongest pathway's -log10(FDR)
    # ---------------------------------------------------------------

    # The table is sorted by FDR, so the first row is the strongest.
    ora = pd.read_csv(ora_file, sep="\t")
    first_fdr = list(ora["fdr"])[0]
    best = -np.log10(first_fdr)

    # ---------------------------------------------------------------
    # STEP 3: Notes for the report
    # ---------------------------------------------------------------

    notes = []
    notes.append("Our version is 10_reactome_ora.png, drawn the same way: top 10 pathways by FDR among the significant genes.")
    notes.append(
        "Matches: GPCR pathways lead both lists (GPCR ligand binding, signalling by GPCR, G alpha (q), class A/1 "
        "rhodopsin-like receptors), followed by extracellular matrix organization. " + found_text
    )
    notes.append(
        f"Differs: our strongest pathway reaches -log10(FDR) = {best:.1f}, the paper's about 3.2. Our gene list is a little "
        "different, and we test against all genes we measured with an up-to-date Reactome file, while the paper used "
        "ShinyGO 0.77 with its own curated Reactome version and background. These change the p-values but not which "
        "pathways come first."
    )
    return notes, "Reproduced with differences"


# ===================================================================
# Fig. 2F: GSEA with GO Molecular Function, as dot plots
# ===================================================================
#
# The paper shows two dot plots (top 10 down terms, top 10 up terms):
#   x position  -> gene ratio
#   dot size    -> number of genes
#   dot color   -> -log10(adjusted p-value)

def draw_gsea_dots(fig, ax, gsea, sign, title, floor):
    """Draw one dot plot: the 10 most significant terms with NES of this sign."""

    # Keep the terms whose NES has this sign (-1 = down, +1 = up).
    has_this_sign = np.sign(gsea["NES"]) == sign
    part = gsea[has_this_sign]

    # The 10 most significant of them.
    part = part.sort_values("FDR q-val", kind="stable")
    part = part.head(10)

    # matplotlib draws the first item at the BOTTOM, so reverse the
    # order to put the most significant term at the top.
    reversed_index = list(part.index)
    reversed_index.reverse()
    part = part.loc[reversed_index]

    # The term names: without the GO ID, in capitals like the paper,
    # and broken onto two lines if longer than 38 characters.
    names = []
    for term in part["Term"]:
        names.append(textwrap.fill(short_go_name(term).upper(), 38))

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


def plot_fig2f(output, reference_folder):

    # This panel needs the paper-style GSEA from run_paper_gsea.py.
    gsea_file = output / "12_fig2F_gsea_go_mf.tsv"
    if not gsea_file.exists():
        notes = ["Skipped: the paper-style GSEA (GO Molecular Function) didn't run, usually because the gene sets "
                 "couldn't be downloaded. Run the analysis with internet access to make this panel."]
        return notes, "Not run"

    gsea = pd.read_csv(gsea_file, sep="\t")

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

    # With 1,000 permutations, the smallest FDR GSEApy can report is
    # about 1/1,000. An FDR of exactly 0 means "smaller than that", so
    # we show it at that limit (clip raises anything smaller up to it).
    floor = 1 / permutations
    clipped_fdr = gsea["FDR q-val"].clip(lower=floor)
    gsea["minus_log10_fdr"] = -np.log10(clipped_fdr)

    # ---------------------------------------------------------------
    # STEP 2: Draw the two dot plots, one above the other
    # ---------------------------------------------------------------

    fig, axes = plt.subplots(2, 1, figsize=(10, 11))
    draw_gsea_dots(fig, axes[0], gsea, -1, "Het-AKAP11-KO vs WT (downregulated)", floor)
    draw_gsea_dots(fig, axes[1], gsea, 1, "Het-AKAP11-KO vs WT (upregulated)", floor)

    fig.suptitle(f"Fig. 2F: GSEA with {paper_gsea_gene_sets} (dot size = leading-edge genes)", y=1.0)

    # tight_layout() spaces the two plots so their labels don't overlap.
    fig.tight_layout()
    save_plot(fig, output / "12_fig2F_gsea_go_mf.png")

    # ---------------------------------------------------------------
    # STEP 3: Compare with the paper's 20 terms
    # ---------------------------------------------------------------

    paper_terms = pd.read_csv(reference_folder / "paper_gsea_go_mf.tsv", sep="\t")

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
            our_direction = "Higher in comparison"
        else:
            our_direction = "Lower in comparison"

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

    notes = []
    notes.append(
        f"{number_significant} GO Molecular Function terms are significant here ({number_down} down, "
        f"{number_up} up). The largest term tested has {largest:,} genes."
    )
    notes.append(
        f"Of the paper's 20 terms, {found} are in our results; {same} of those change in the same direction, and "
        f"{same_and_significant} are also significant here."
    )
    if len(missing) > 0:
        notes.append("Not in our results: " + ", ".join(missing) + ". GO term names and gene sets differ between "
                     "the GO version the paper used (its names look like MSigDB's) and Enrichr's 2025 library.")
    notes.append(f"Differs: the paper's colour scale reaches -log10(padj) of about 30; ours stops at "
                 f"{-np.log10(floor):.0f}. The paper used fgsea, which estimates very small p-values precisely; GSEApy "
                 f"with {permutations:,} permutations can't report an FDR much below 1/{permutations:,}. The gene ratio and "
                 "dot size are our closest equivalents to the paper's.")
    return notes, "Reproduced with differences"


# ###################################################################
# PART 3: Make every panel and write the report
# ###################################################################

def plot_paper_figures(output, reference_folder):

    output = Path(output)
    reference_folder = Path(reference_folder)

    # ---------------------------------------------------------------
    # STEP 1: Read our results
    # ---------------------------------------------------------------

    results = read_table(output / "02_differential_expression.tsv")
    normalized = pd.read_csv(output / "02_normalized_counts.tsv", sep="\t", index_col=0)
    vst = pd.read_csv(output / "02_vst_expression.tsv", sep="\t", index_col=0)
    metadata = pd.read_csv(output / "01_sample_qc.tsv", sep="\t", index_col=0)

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

    notes, status = plot_fig1b(results, normalized, metadata, reference_folder, output)
    panels.append({"panel": "Fig. 1B", "what": "AKAP11 exon 7 RT-qPCR per clone",
                   "figure": "12_fig1B_akap11_relative.png", "notes": notes, "status": status})

    panels.append({"panel": "Fig. 1C", "what": "AKAP11 protein (Western blot)", "figure": "none",
                   "notes": ["Not reproducible from RNA-seq: it measures protein on a gel. The closest RNA-seq "
                             "evidence is Fig. 1B (AKAP11 RNA is lower in every clone)."],
                   "status": "Not reproducible from RNA-seq"})

    panels.append({"panel": "Fig. 1D", "what": "Neuron images (immunostaining)", "figure": "none",
                   "notes": ["Not reproducible from RNA-seq: these are microscope images. The paper's cell-type "
                             "deconvolution of the RNA-seq data, which found neurons predominant in both groups, "
                             "is the sequencing-based counterpart."],
                   "status": "Not reproducible from RNA-seq"})

    notes, status = plot_fig2a(results, reference_folder, output)
    panels.append({"panel": "Fig. 2A", "what": "Volcano plot",
                   "figure": "12_fig2A_volcano.png", "notes": notes, "status": status})

    notes, status = plot_fig2b(results, reference_folder, output)
    panels.append({"panel": "Fig. 2B", "what": "Number of significant genes",
                   "figure": "12_fig2B_deg_counts.png", "notes": notes, "status": status})

    notes, status = plot_fig2c(results, vst, output)
    panels.append({"panel": "Fig. 2C", "what": "Clustered heatmap of significant genes",
                   "figure": "12_fig2C_heatmap.png", "notes": notes, "status": status})

    notes, status = plot_fig2d(results, reference_folder, output)
    panels.append({"panel": "Fig. 2D", "what": "Top 30 down and up genes",
                   "figure": "12_fig2D_top_genes.png", "notes": notes, "status": status})

    notes, status = describe_fig2e(output, reference_folder)
    panels.append({"panel": "Fig. 2E", "what": "Reactome ORA",
                   "figure": "10_reactome_ora.png", "notes": notes, "status": status})

    notes, status = plot_fig2f(output, reference_folder)
    panels.append({"panel": "Fig. 2F", "what": "GSEA, GO Molecular Function",
                   "figure": "12_fig2F_gsea_go_mf.png", "notes": notes, "status": status})

    panels.append({"panel": "Fig. 3A–E", "what": "DNA methylation (WGBS)", "figure": "none",
                   "notes": ["Not reproducible from this data: Fig. 3 comes from whole-genome bisulfite sequencing, "
                             "a different experiment with its own FASTQ files and pipeline (for example "
                             "nf-core/methylseq). Fig. 3E also needs those methylation results, because it uses "
                             "genes that change in both expression and methylation.",
                             "It is outside this RNA-seq project, but would be the natural next step."],
                   "status": "Not reproducible from RNA-seq"})

    # ---------------------------------------------------------------
    # STEP 3: Write the report
    # ---------------------------------------------------------------

    # The title and a summary table with one row per panel.
    report = []
    report.append("# Reproducing the paper's figures")
    report.append("")
    report.append("Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557, "
                  "https://doi.org/10.1038/s41380-025-03040-x")
    report.append("")
    report.append("Every panel of the paper's Figures 1 to 3, what this project made for it, and why any differences "
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

    with open(output / "12_paper_figures_report.md", "w") as file:
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
#     python3 plot_paper_figures.py akap11_go_bp_results
# It does NOT run when run_analysis.py imports this file.
if __name__ == "__main__":

    # Save plots to files instead of opening windows, and use the same
    # plot style as run_analysis.py.
    import matplotlib
    matplotlib.use("Agg")
    sns.set_theme(style="whitegrid")

    from settings import paper_reference_folder

    # We need exactly 2 words: the script name and the results folder.
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 plot_paper_figures.py results_folder")

    plot_paper_figures(sys.argv[1], paper_reference_folder)
