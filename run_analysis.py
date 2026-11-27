"""
===================================================================
run_analysis.py      <-- START HERE
===================================================================

WHAT THIS FILE DOES
    This is the main file. It runs the whole analysis from start to
    finish. It doesn't do the calculations itself: it reads what you
    typed in Terminal, then calls the other files one after another.

        1. read_counts        -> read the count table and sample
                                 sheet, and remove low-count genes
        2. fit_model          -> PyDESeq2: find the genes that differ
                                 between the two groups
        3. plot_pca           -> PCA plot and sample correlation
                                 heatmap
        4. plot_volcano       -> volcano plot (it also calls plot_ma
                                 for the MA plot)
        5. plot_expression_heatmap -> heatmap of the top genes
        6. plot_target_gene   -> bar plot for one gene (optional)
        7. run_gsea           -> gene set enrichment analysis (GSEA)
        8. plot_pathways      -> bar plot of the top GO terms, and
           plot_gsea_enrichment -> classic GSEA enrichment plots
        9. run_ora            -> Reactome over-representation analysis
       10. plot_ora           -> bar plot of the top Reactome pathways
       11. run_paper_gsea     -> a second GSEA set up like the paper's
                                 Fig. 2F
       12. compare_to_paper   -> check our results against the paper
       13. plot_paper_figures -> redraw the paper's figures from our
                                 results

    Steps 11 to 13 run only if the paper_reference folder exists, so
    the same code still works on a different dataset.

HOW TO RUN IT (type it all on one line in Terminal)

    python3 run_analysis.py counts.tsv samples.tsv WT KO results_folder AKAP11

        counts.tsv      the gene count table from nf-core/rnaseq
        samples.tsv     the sample sheet (sample name -> group)
        WT              the baseline group
        KO              the group compared with the baseline
        results_folder  where every output file is saved
        AKAP11          OPTIONAL: one gene to draw a bar plot for

OUTPUT
    Every file in the results folder, numbered by step (01_ to 12_).
    The settings for every step (cutoffs, gene set collections, and
    so on) are in settings.py.

PYTHON IDEAS USED IN THIS FILE
    sys.argv        the list of words typed in Terminal
    sys.exit(text)  print the text and stop the program
    from x import y use the function y written in the file x.py
    a, b, c = f()   a function can give back several things at once;
                    each one gets its own name
===================================================================
"""

# sys lets us read what was typed in Terminal, and stop the program.
import sys

# Path makes file and folder paths that work on every computer.
from pathlib import Path

import matplotlib

# Tell Matplotlib to save plots to files instead of opening windows.
# This line must come BEFORE we import anything that draws plots.
matplotlib.use("Agg")

import seaborn as sns

# Each line below brings in ONE function from one of our own files.
# For example, "from read_counts import read_counts" means: open
# read_counts.py and use the function called read_counts() from it.
from read_counts import read_counts
from fit_model import fit_model
from plot_pca import plot_pca
from plot_volcano import plot_volcano
from plot_expression_heatmap import plot_expression_heatmap
from plot_target_gene import plot_target_gene
from run_gsea import run_gsea
from plot_pathways import plot_pathways
from plot_gsea_enrichment import plot_gsea_enrichment
from run_ora import run_ora
from plot_ora import plot_ora
from compare_to_paper import compare_to_paper
from run_paper_gsea import run_paper_gsea
from plot_paper_figures import plot_paper_figures

# A few settings from settings.py.
from settings import heatmap_genes            # how many genes in the heatmap (50)
from settings import gsea_gene_sets           # which GO collection GSEA tests
from settings import paper_reference_folder   # where the paper's numbers are
from settings import reference_color          # the color of the reference group (WT)
from settings import comparison_color         # the color of the comparison group (KO)


# ===================================================================
# STEP 1: Read what the user typed in Terminal
# ===================================================================

# This message explains how to use the script. We show it if the
# user types the wrong number of words.
usage = (
    "Usage: python run_analysis.py "
    "counts.tsv samples.tsv reference comparison output [target_gene]"
)

# sys.argv is a list of every word typed in Terminal, for example
#   ["run_analysis.py", "counts.tsv", "samples.tsv", "WT", "KO",
#    "results", "AKAP11"]
# sys.argv[0] is the script's own name, so we skip it: [1:] means
# "from item 1 to the end".
arguments = sys.argv[1:]

number_of_arguments = len(arguments)

# If the user typed only --help, show the usage message and stop.
# sys.exit(0) stops the program and means "no error".
if arguments == ["--help"]:
    print(usage)
    sys.exit(0)

# We need 5 words, or 6 if a target gene is given.
# Anything else is a mistake, so show the usage message and stop.
if number_of_arguments != 5 and number_of_arguments != 6:
    sys.exit(usage)

# Give each word a clear name.
# Python counts from 0, so the 1st word is arguments[0].
# Path() turns text into a file path that Python can work with.
count_file = Path(arguments[0])        # the count table
sample_file = Path(arguments[1])       # the sample sheet
reference_group = arguments[2]         # the baseline group, e.g. WT
comparison_group = arguments[3]        # the group we compare, e.g. KO
output = Path(arguments[4])            # the folder for the results

# The target gene is optional. It's the 6th word, if there is one,
# which is arguments[5]. None is Python's word for "nothing".
if number_of_arguments == 6:
    target_gene = arguments[5]
else:
    target_gene = None

# The two groups must be different, or there's nothing to compare.
if reference_group == comparison_group:
    sys.exit("Choose two different conditions.")


# ===================================================================
# STEP 2: Set up the output folder and the plot colors
# ===================================================================

# Make the results folder.
#   parents=True  -> also make any missing parent folders
#   exist_ok=True -> don't complain if the folder already exists
output.mkdir(parents=True, exist_ok=True)

# Use seaborn's clean white style with light grid lines for all plots.
sns.set_theme(style="whitegrid")

# Pick one color per group, in a dictionary: group name -> color.
# Every plot uses this same dictionary, so WT is always blue and KO is
# always red. The colors themselves come from settings.py; the group
# names come from what the user typed, so nothing here is tied to "WT"
# or "KO".
colors = {}
colors[reference_group] = reference_color     # e.g. "WT" -> "steelblue"
colors[comparison_group] = comparison_color   # e.g. "KO" -> "firebrick"


# ===================================================================
# STEP 3: Read the data and fit the model
# ===================================================================

# read_counts() gives back three tables:
#   counts      the filtered counts (genes x samples)
#   metadata    which group each sample is in
#   annotation  each gene's ID and name
counts, metadata, annotation = read_counts(
    count_file,
    sample_file,
    reference_group,
    comparison_group,
)

# fit_model() runs PyDESeq2 and gives back three more tables:
#   results      one row per gene: fold change, p-value, padj, ...
#   normalized   counts corrected for sequencing depth
#   transformed  VST values, used for PCA and heatmaps
# It also saves the 02_ files.
results, normalized, transformed = fit_model(
    counts,
    metadata,
    annotation,
    reference_group,
    comparison_group,
    output,
)


# ===================================================================
# STEP 4: Make the sample plots and the gene plots
# ===================================================================

# PCA plot and sample correlation heatmap (03_ files).
plot_pca(count_file, transformed, metadata, colors, output)

# Volcano plot (05_); this also makes the MA plot (06_).
plot_volcano(results, reference_group, comparison_group, colors, output)

# Clustered heatmap of the top genes (07_). It gets the sample groups
# and colors too, to draw a WT/KO color bar above the columns.
plot_expression_heatmap(results, transformed, metadata, colors, heatmap_genes, output)

# Bar plot for one gene (08_), but only if the user asked for one
# (and didn't just type an empty word, "").
if target_gene is not None and target_gene != "":
    plot_target_gene(results, normalized, metadata, target_gene, colors, output)


# ===================================================================
# STEP 5: Run GSEA and plot the top GO terms
# ===================================================================

# gsea_gene_sets comes from settings.py. It says which GO collection
# to test, for example "GO_Biological_Process_2025".
# run_gsea() saves 09_gsea_pathways.tsv and gives back two things:
#   pathways   -> the results table
#   enrichment -> GSEApy's full result (with each term's curve)
pathways, enrichment = run_gsea(results, gsea_gene_sets, output)

# Bar plot of the top terms (09_gsea_pathways.png).
plot_pathways(pathways, reference_group, comparison_group, colors, output)

# Classic GSEA enrichment plots of the top terms in each direction
# (09_gsea_enrichment_plots.png).
plot_gsea_enrichment(pathways, enrichment, reference_group, comparison_group, output)


# ===================================================================
# STEP 6: Run Reactome ORA and plot the top pathways
# ===================================================================

# run_ora() saves 10_reactome_ora.tsv and gives back the table.
ora = run_ora(results, output)

# Only draw the plot if at least one pathway was tested.
if len(ora) > 0:
    plot_ora(ora, output)


# ===================================================================
# STEP 7: Compare our results with the published paper,
#         and redraw the paper's figures
# ===================================================================

# Path(...).exists() is True if the folder is there.
# For a different dataset there's no paper_reference folder, so this
# step is skipped and everything else still runs.
reference_folder = Path(paper_reference_folder)

if reference_folder.exists():
    # A second GSEA with the paper's GO collection, for Fig. 2F
    # (12_fig2F_gsea_go_mf.tsv). It needs internet access.
    run_paper_gsea(output)

    # Check our numbers against the paper's (11_ files).
    compare_to_paper(output, reference_folder)

    # Redraw the paper's figures and explain any differences
    # (12_ files). It needs the two group names so it knows which
    # samples are the reference (WT) and which the comparison (KO).
    plot_paper_figures(output, reference_folder, reference_group, comparison_group)
else:
    print("No paper_reference folder found, so the comparison with the paper was skipped.")

print(f"Results saved in {output}")
