"""
===================================================================
settings.py
===================================================================

WHAT THIS FILE DOES
    It holds every number, name, color, and file name that the analysis
    uses: cutoffs, column names, labels, colors, and the names of the
    files it reads and writes.

    Nothing "runs" in this file. It only creates variables
    (a variable is a name that stores a value, like  min_count = 10).

WHERE IT FITS
    Every other script reads its settings from here. For example,
    read_counts.py starts with the line:

        from settings import min_count

    which means: "open settings.py and give me the value of min_count".

WHY HAVE A SETTINGS FILE? (WHY NOTHING IS "HARD-CODED")
    "Hard-coding" means typing a value straight into the code where it
    is used, for example, writing 0.05 inside five different scripts.
    That causes two problems:
      1. To change the value you must find and edit every copy, and if
         you miss one, the scripts silently disagree with each other.
      2. A reader can't see all the choices the analysis makes in one
         place.
    Here every choice is typed ONCE, with a comment saying why it was
    chosen. Every script that needs it imports it from here.

    Facts about the PUBLISHED PAPER (its numbers, its sample order,
    its methods, its citation) are not settings: they are data copied
    from the paper, so they live in the files of the paper_reference
    folder instead (see paper_reference/README.md).

WHAT IS DELIBERATELY NOT IN THIS FILE, AND WHY
    - Column names made by other tools, such as PyDESeq2's "padj" and
      "log2FoldChange" or GSEApy's "NES" and "FDR q-val". Those tools
      decide these names, not us, so they can't be changed here.
    - Column names of the tables our own code creates, such as
      "direction" or "label". They work like variable names: they are
      part of how the code is written, not choices about the analysis.
    - Purely visual numbers (figure sizes in inches, font sizes, dot
      sizes, transparency). They change how a plot looks, never what
      it shows, so each one stays next to the plot line it controls,
      with a comment explaining it.

HOW TO READ A LINE IN THIS FILE
        fdr_cutoff = 0.05
        ^^^^^^^^^^   ^^^^
        the name     the value stored under that name

    Values can be:
        numbers         10, 0.05, 1000
        text (strings)  "gene_id"      <- text always has quote marks
        lists           ["a", "b"]     <- several values in [ ]
        None            a special Python value meaning "nothing chosen"
===================================================================
"""


# ###################################################################
# PART 1: Column names in the two input tables
# ###################################################################

# The count table (AKAP11_final_gene_counts_length_scaled.tsv) looks
# like this (only the first few columns are shown):
#   gene_id              gene_name   AKAP11_KO_clone16   ...   WT_rep3
#   ENSG00000000003.16   TSPAN6      9775.336            ...   9120.512
#   ENSG00000000005.6    TNMD        0.0                 ...   1.0

# The sample sheet (samples.tsv) looks like this:
#   sample              condition   geo_sample   sra_run
#   AKAP11_KO_clone16   KO          GSM8200086   SRR28641569

# These settings tell the scripts which columns hold what. If your own
# files use other column names, change them here and nowhere else.

# The gene ID column. A gene ID (like "ENSG00000000003.16") is unique:
# no two genes share one. That makes it the safest way to match the
# same gene across different tables.
gene_id_column = "gene_id"

# The gene name column. A gene name, also called a gene symbol
# (like "TSPAN6"), is easier for people to read, so we use it for plot
# labels. GSEA and ORA also match genes to pathways by these names.
gene_name_column = "gene_name"

# The sample sheet column with each sample's name. These names must
# match the sample column headings of the count table exactly.
sample_column = "sample"

# The sample sheet column with each sample's group (for us, WT or KO).
# The DESeq2 design formula is built from this name (fit_model.py).
condition_column = "condition"


# ###################################################################
# PART 2: Removing genes with very few reads
# ###################################################################

# Most of the 78,932 genes in the table have almost no reads in these
# neurons. A gene with 0 to 3 reads per sample can't be tested
# reliably, so we remove those genes before the statistics.

# A sample "passes" for a gene if the gene has at least this many
# reads in that sample (after rounding to whole numbers).
# 10 is the cutoff recommended in the DESeq2 / edgeR tutorials: below
# about 10 reads, random sampling noise is as big as any real change.
min_count = 10

# How many samples must pass for us to keep the gene?

# None means "not chosen here": read_counts.py then uses the size of
# the smaller group instead.
# Example with our data:
#   3 WT samples + 3 KO samples -> the smaller group has 3 samples.
#   So a gene is kept if ANY 3 of the 6 samples have at least 10 reads.
# Why "any 3" and not "all 6"? A gene that is switched ON only in the
# knockout (3 samples with reads, 3 without) is exactly the kind of
# gene we want to find. Requiring all 6 would throw it away.
min_samples = None

# Each group needs at least this many samples (replicates).
# Why 2? With only 1 sample per group, there is no way to measure how
# much samples of the same group naturally vary, so no statistical
# test is possible. 2 is the mathematical minimum (3 or more is much
# better). Don't set this below 2.
min_replicates_per_group = 2


# ###################################################################
# PART 3: Cutoffs for calling a gene "significant"
# ###################################################################

# A gene is called significant only if it passes BOTH cutoffs below.
# These are the same cutoffs the paper used, so our gene lists can be
# compared with theirs. (plot_paper_figures.py checks this against
# paper_reference/paper_summary.tsv.)

# CUTOFF 1: the adjusted p-value (padj).
# A p-value says how surprising the result would be if the gene did
# NOT really change. We test about 26,000 genes, so even with no real
# changes, about 1,300 of them (5%) would get p < 0.05 by luck alone.
# The adjusted p-value corrects for testing so many genes at once.
# A gene needs padj of 0.05 or less.
# The same cutoff is also used for pathways and GO terms (ORA, GSEA),
# so "significant" means the same thing everywhere in the project.
fdr_cutoff = 0.05

# CUTOFF 2: the size of the change, as a log2 fold change.
# "Fold change" = expression in KO divided by expression in WT.
# "log2" turns fold changes into a symmetric scale:
#     fold change     log2 fold change     meaning
#     2  (twice)       1                    doubled in KO
#     1  (same)        0                    no change
#     0.5 (half)      -1                    halved in KO
# A gene needs a log2 fold change of at least 1 in either direction,
# meaning it at least doubled or at least halved.
fold_change_cutoff = 1.0

# ---------------------------------------------------------------
# The three labels every gene gets (fit_model.py adds them in a
# column called "direction"). Many scripts look for these labels, so
# they are typed only here: a typo in one script can't break it.

# "comparison" means the comparison group (KO). The labels don't say
# "KO" on purpose, so the same code works for any two groups.

# NOTE: paper_reference/paper_gsea_go_mf.tsv uses the same two
# "Higher/Lower" labels. If you change them here, change them there.
# ---------------------------------------------------------------
higher_label = "Higher in comparison"         # significant AND higher in KO
lower_label = "Lower in comparison"           # significant AND lower in KO
not_significant_label = "Below thresholds"    # failed at least one cutoff


# ###################################################################
# PART 4: Plot settings
# ###################################################################

# How many gene names to write on EACH side of the volcano plot.
# With more labels, the names pile on top of each other.
volcano_labels = 8

# How many genes to show in the gene heatmap (07_de_heatmap.png).
# 50 rows still leaves room for every gene name to be readable.
heatmap_genes = 50

# In both heatmaps, colors stop getting darker beyond this z-score
# (in either direction). Why 2? About 95% of values in a normal
# distribution lie between -2 and +2, so this keeps the color range
# for typical values, and one extreme gene can't wash out the rest.
heatmap_z_limit = 2

# A p-value of exactly 0 can't be put on a -log10 axis (log10(0) is
# minus infinity). Before plotting, p-values smaller than this are
# raised to it. 1e-300 is close to the smallest number a computer can
# store, so it changes nothing except true zeros. The tables are not
# changed, only the plots.
smallest_plotted_p_value = 1e-300

# Image sharpness for every saved figure, in dots per inch (dpi).
# 200 gives a sharp image with a modest file size. (Screens are about
# 100 dpi; printing usually needs 300.) Used by save_plot.py.
image_resolution = 200

# ---------------------------------------------------------------
# Colors. Matplotlib understands color names ("steelblue") and hex
# codes ("#b9bfc6": red, green and blue amounts, from 00 to ff).
# Using one variable per role means the same group always gets the
# same color in every plot.
# ---------------------------------------------------------------
reference_color = "steelblue"     # the reference group (WT): blue
comparison_color = "firebrick"    # the comparison group (KO): red
not_significant_color = "#b9bfc6" # genes below the cutoffs: light gray
neutral_color = "slategray"       # plots with no group (ORA bars, paper-vs-ours dots)
cutoff_line_color = "gray"        # dashed lines that mark the cutoffs
label_line_color = "#65717a"      # thin lines from a moved label to its dot


# ###################################################################
# PART 5: GSEA (gene set enrichment analysis)
# ###################################################################

# GSEA puts ALL genes in order, from "most up in KO" to "most down in
# KO", then checks whether the genes of a GO term (a named group of
# genes, like "Translation") bunch up near the top or the bottom.

# Which collection of gene sets (GO terms) to test.

# Two kinds of value work here:
#   1. An Enrichr library name. GSEApy downloads it, so you need
#      internet. Two useful ones:
#        "GO_Biological_Process_2025"  -> processes genes take part in
#                                         (example: "Translation")
#        "GO_Molecular_Function_2025"  -> what gene products do
#                                         (example: "Actin Binding")
#   2. The name of a .gmt file in this folder. This works offline:
#        "go_biological_process.gmt"

gsea_gene_sets = "GO_Biological_Process_2025"

# A short, readable name for the collection above. It is only used in
# the GSEA plot's title. If you change gsea_gene_sets, change this too,
# for example to "GO Molecular Function".
gsea_collection_label = "GO Biological Process"

# Only test GO terms that share between 15 and 500 genes with our
# ranked list. These are GSEA's own default limits.
#   - Terms with fewer than 15 genes are unreliable: a couple of genes
#     can make them look significant.
#   - Terms with more than 500 genes are very general (like "protein
#     binding") and tell us little.
gsea_min_size = 15
gsea_max_size = 500

# How many GO terms to show in EACH direction in the GSEA bar plot.
# (The results table 09_gsea_pathways.tsv still lists every term.)
gsea_top = 8

# GSEA works out what a "random" result looks like by shuffling the
# data many times and repeating the test. Each shuffle is called a
# permutation. 1,000 is the GSEA default; more takes longer.
permutations = 1000

# The shuffles are random. Setting a "seed" makes the computer produce
# the same random numbers every run, so you get exactly the same
# results each time. The number 42 itself has no special meaning.
seed = 42


# ###################################################################
# PART 6: ORA (over-representation analysis)
# ###################################################################

# ORA asks a simpler question than GSEA: "Among our significant genes,
# are the genes of some pathway found more often than chance would
# predict?" The paper did this with Reactome pathways (its Fig. 2E).

# The pathway file to use. It must be in this folder.
ora_gene_sets = "reactome_human.gmt"

# A short, readable name for that collection, used in plot titles and
# report text. If you change ora_gene_sets, change this too.
ora_collection_label = "Reactome"

# Only test pathways with 2 to 2,000 of our genes in them. These are
# the default limits of ShinyGO, the tool the paper used.
ora_min_size = 2
ora_max_size = 2000

# How many pathways to show in the ORA bar plot (like the paper's
# Fig. 2E, which shows 10).
ora_top = 10

# Pathway and GO term names often end with a database ID, e.g.
#   "Extracellular Matrix Organization R-HSA-1474244"   (Reactome)
#   "Structural Constituent of Ribosome (GO:0003735)"   (Gene Ontology)
# The IDs make labels long and stop names from matching between
# databases, so the code cuts each name at the first of these markers.
# To use another database, add the start of its IDs to this list.
term_id_markers = [" R-HSA-", "(GO:"]


# ###################################################################
# PART 7: Comparing with the published paper (compare_to_paper.py)
# ###################################################################

# The folder with small tables copied from the paper's figures.
# compare_to_paper.py and plot_paper_figures.py check our results
# against them.

# If the folder doesn't exist (for example, when you analyze a
# different dataset), those comparison steps are skipped and
# everything else still runs.
paper_reference_folder = "paper_reference"

# ---------------------------------------------------------------
# Rules of thumb for the "verdict" of each check.
# There is no official rule for when two analyses "agree", so these
# are our choices, written down so anyone can see (and change) them.
# ---------------------------------------------------------------

# Check 1: the two totals of significant genes count as "similar" if
# they are at most 20% apart. Why 20%? Different counting tools and
# statistical models routinely move genes near the cutoffs; studies
# comparing DESeq2 and limma-voom report differences of this size.
similar_total_percent = 20

# Check 2 (and the GO term check): "agree" if at least 80% of the
# paper's genes (or terms) change in the same direction in ours.
# By pure chance, about 50% would, so 80% is far above chance.
agree_direction_percent = 80

# Check 3: "agree" if at least 70% of the paper's top Reactome
# pathways are also significant in our ORA. This is a bit lower than
# the 80% above because pathway p-values depend on the database
# version and the background genes, which differ (see the methods
# table in paper_reference/paper_methods.tsv).
agree_ora_percent = 70

# How many example terms to list in the report under each finding.
report_examples = 3

# A correlation needs at least 3 points to mean anything (with 2
# points, any two values lie on a perfect straight line).
min_points_for_correlation = 3


# ###################################################################
# PART 8: Redrawing the paper's figures (plot_paper_figures.py)
# ###################################################################

# The paper's Fig. 2F is a GSEA with GO Molecular Function terms,
# while our main GSEA (PART 5) uses GO Biological Process. To redraw
# Fig. 2F, run_paper_gsea.py runs a SECOND GSEA with this collection.
# GSEApy downloads it from Enrichr, so you need internet.
paper_gsea_gene_sets = "GO_Molecular_Function_2025"

# A short, readable name for that collection, used in plot titles and
# the report (like gsea_collection_label in PART 5).
paper_gsea_collection_label = "GO Molecular Function"

# Size limits for that second GSEA.
# The paper's Fig. 2F includes very large terms such as "RNA binding",
# which has more than 1,000 genes. Our main GSEA stops at 500 genes
# (gsea_max_size above), which would leave those terms out. So this
# second GSEA allows terms with up to 3,000 genes.
paper_gsea_min_size = 15
paper_gsea_max_size = 3000

# Colors close to the paper's own, so the redrawn figures can be put
# side by side with the published ones.
paper_down_color = "#2F5DA8"         # blue: lower in the knockout
paper_up_color = "#B2182B"           # red: higher in the knockout
paper_gray_color = "#C9CCD1"         # gray: not significant
paper_wt_color = "#2E8B57"           # Fig. 1B: green for WT
paper_ko_color = "#3F51B5"           # Fig. 1B: blue for the KO clones
paper_highlight_color = "#FFF1C2"    # Fig. 2D: pale yellow = also in the paper's list
paper_heatmap_colors = ["#4B3F9E", "#FFFFFF", "#D7301F"]   # Fig. 2C: purple, white, red


# ###################################################################
# PART 9: Names of the files in the paper_reference folder
# ###################################################################
#
# What each file holds is explained in paper_reference/README.md.

paper_info_file = "paper_info.tsv"                    # citation and a few facts about the paper
paper_summary_file = "paper_summary.tsv"              # numbers from Fig. 2B (and read off other panels)
paper_methods_file = "paper_methods.tsv"              # the paper's methods next to ours
paper_genes_file = "paper_degs.tsv"                   # Fig. 2D: the paper's top genes
paper_volcano_labels_file = "paper_volcano_labels.tsv"  # Fig. 2A: genes labeled on the volcano plot
paper_ora_file = "paper_reactome_ora.tsv"             # Fig. 2E: the paper's top Reactome pathways
paper_gsea_file = "paper_gsea_go_mf.tsv"              # Fig. 2F: the paper's GO terms and directions
paper_fig1b_file = "paper_fig1b_qpcr.tsv"             # Fig. 1B: knocked-out gene in each clone
paper_fig2c_columns_file = "paper_fig2c_columns.tsv"  # Fig. 2C: the paper's column order and labels
paper_not_reproducible_file = "paper_not_reproducible.tsv"  # panels RNA-seq can't redraw, and why


# ###################################################################
# PART 10: Names of the files the analysis writes
# ###################################################################
#
# Every output file name is typed once, here. Several scripts read
# files that another script wrote (for example, compare_to_paper.py
# reads 02_differential_expression.tsv, which fit_model.py writes),
# so the writer and the reader must use exactly the same name.
# The numbers at the start (01_, 02_, ...) are the step that makes the
# file, so the results folder lists them in the order they were made.

sample_qc_file = "01_sample_qc.tsv"
differential_expression_file = "02_differential_expression.tsv"
significant_genes_file = "02_significant_genes.tsv"
normalized_counts_file = "02_normalized_counts.tsv"
vst_expression_file = "02_vst_expression.tsv"
pca_plot_file = "03_pca.png"
pca_scores_file = "03_pca_scores.tsv"
sample_correlation_plot_file = "04_sample_correlation.png"
volcano_plot_file = "05_volcano.png"
volcano_labels_file = "05_volcano_labels.tsv"
ma_plot_file = "06_ma.png"
heatmap_plot_file = "07_de_heatmap.png"
target_gene_plot_file = "08_target_gene.png"
target_gene_results_file = "08_target_gene_results.tsv"
gsea_ranking_file = "09_gsea_gene_ranking.rnk"
gsea_results_file = "09_gsea_pathways.tsv"
gsea_plot_file = "09_gsea_pathways.png"
ora_results_file = "10_reactome_ora.tsv"
ora_plot_file = "10_reactome_ora.png"
comparison_report_file = "11_comparison_with_paper.md"
gene_comparison_file = "11_paper_gene_comparison.tsv"
term_comparison_file = "11_paper_term_comparison.tsv"
fold_change_comparison_plot_file = "11_paper_vs_our_fold_change.png"
fig1b_table_file = "12_fig1B_target_gene_relative.tsv"
fig1b_plot_file = "12_fig1B_target_gene_relative.png"
fig2a_plot_file = "12_fig2A_volcano.png"
fig2b_plot_file = "12_fig2B_deg_counts.png"
fig2c_plot_file = "12_fig2C_heatmap.png"
fig2d_table_file = "12_fig2D_top_genes.tsv"
fig2d_plot_file = "12_fig2D_top_genes.png"
paper_gsea_results_file = "12_fig2F_gsea_go_mf.tsv"
fig2f_plot_file = "12_fig2F_gsea_go_mf.png"
figures_report_file = "12_paper_figures_report.md"
