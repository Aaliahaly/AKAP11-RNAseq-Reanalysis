"""
===================================================================
settings.py
===================================================================

WHAT THIS FILE DOES
    It holds every number and name that the analysis uses:
    cutoffs, column names, plot sizes, file names.

    Nothing "runs" in this file. It only creates variables
    (a variable is a name that stores a value, like  min_count = 10).

WHERE IT FITS
    Every other script reads its settings from here. For example,
    read_counts.py starts with the line:

        from settings import min_count

    which means: "open settings.py and give me the value of min_count".

WHY HAVE A SETTINGS FILE?
    If a cutoff were typed inside five different scripts, changing it
    would mean finding and editing all five. Here you change it once,
    and every script that uses it picks up the new value.

HOW TO READ A LINE IN THIS FILE
        fdr_cutoff = 0.05
        ^^^^^^^^^^   ^^^^
        the name     the value stored under that name

    Values can be:
        numbers         10, 0.05, 1000
        text (strings)  "gene_id"      <- text always has quote marks
        None            a special Python value meaning "nothing chosen"
===================================================================
"""


# ###################################################################
# PART 1: Column names in the count table
# ###################################################################
#
# The count table (AKAP11_final_gene_counts_length_scaled.tsv) looks
# like this (only the first few columns are shown):
#
#   gene_id              gene_name   AKAP11_KO_clone16   ...   WT_rep3
#   ENSG00000000003.16   TSPAN6      9775.336            ...   9120.512
#   ENSG00000000005.6    TNMD        0.0                 ...   1.0
#
# These two settings tell the scripts which columns hold the gene
# IDs and the gene names.

# The gene ID column. A gene ID (like "ENSG00000000003.16") is unique:
# no two genes share one. That makes it the safest way to match the
# same gene across different tables.
gene_id_column = "gene_id"

# The gene name column. A gene name, also called a gene symbol
# (like "TSPAN6"), is easier for people to read, so we use it for plot
# labels. GSEA and ORA also match genes to pathways by these names.
gene_name_column = "gene_name"


# ###################################################################
# PART 2: Removing genes with very few reads
# ###################################################################
#
# Most of the 78,932 genes in the table have almost no reads in these
# neurons. A gene with 0 to 3 reads per sample can't be tested
# reliably, so we remove those genes before the statistics.

# A sample "passes" for a gene if the gene has at least this many
# reads in that sample (after rounding to whole numbers).
min_count = 10

# How many samples must pass for us to keep the gene?
#
# None means "not chosen here": read_counts.py then uses the size of
# the smaller group instead.
#
# Example with our data:
#   3 WT samples + 3 KO samples -> the smaller group has 3 samples.
#   So a gene is kept if ANY 3 of the 6 samples have at least 10 reads.
#
# Why "any 3" and not "all 6"? A gene that is switched ON only in the
# knockout (3 samples with reads, 3 without) is exactly the kind of
# gene we want to find. Requiring all 6 would throw it away.
min_samples = None


# ###################################################################
# PART 3: Cutoffs for calling a gene "significant"
# ###################################################################
#
# A gene is called significant only if it passes BOTH cutoffs below.
# These are the same cutoffs the paper used, so our gene lists can be
# compared with theirs.

# CUTOFF 1: the adjusted p-value (padj).
#
# A p-value says how surprising the result would be if the gene did
# NOT really change. We test about 26,000 genes, so even with no real
# changes, about 1,300 of them (5%) would get p < 0.05 by luck alone.
# The adjusted p-value corrects for testing so many genes at once.
# A gene needs padj of 0.05 or less.
fdr_cutoff = 0.05

# CUTOFF 2: the size of the change, as a log2 fold change.
#
# "Fold change" = expression in KO divided by expression in WT.
# "log2" turns fold changes into a symmetric scale:
#
#     fold change     log2 fold change     meaning
#     2  (twice)       1                    doubled in KO
#     1  (same)        0                    no change
#     0.5 (half)      -1                    halved in KO
#
# A gene needs a log2 fold change of at least 1 in either direction,
# meaning it at least doubled or at least halved.
fold_change_cutoff = 1.0


# ###################################################################
# PART 4: Plot settings
# ###################################################################

# How many gene names to write on EACH side of the volcano plot.
# With more labels, the names pile on top of each other.
volcano_labels = 8

# How many genes to show in the gene heatmap (07_de_heatmap.png).
# 25 rows still leaves room for every gene name to be readable.
heatmap_genes = 25


# ###################################################################
# PART 5: GSEA (gene set enrichment analysis)
# ###################################################################
#
# GSEA puts ALL genes in order, from "most up in KO" to "most down in
# KO", then checks whether the genes of a GO term (a named group of
# genes, like "Translation") bunch up near the top or the bottom.

# Which collection of gene sets (GO terms) to test.
#
# Two kinds of value work here:
#
#   1. An Enrichr library name. GSEApy downloads it, so you need
#      internet. Two useful ones:
#        "GO_Biological_Process_2025"  -> processes genes take part in
#                                         (example: "Translation")
#        "GO_Molecular_Function_2025"  -> what gene products do
#                                         (example: "Actin Binding")
#
#   2. The name of a .gmt file in this folder. This works offline:
#        "go_biological_process.gmt"
gsea_gene_sets = "GO_Biological_Process_2025"

# A short, readable name for the collection above. It is only used in
# the GSEA plot's title. If you change gsea_gene_sets, change this too,
# for example to "GO Molecular Function".
gsea_collection_label = "GO Biological Process"

# Only test GO terms that share between 15 and 500 genes with our
# ranked list.
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
# permutation. 1,000 is a common choice; more takes longer.
permutations = 1000

# The shuffles are random. Setting a "seed" makes the computer produce
# the same random numbers every run, so you get exactly the same
# results each time. The number 42 itself has no special meaning.
seed = 42


# ###################################################################
# PART 6: ORA (over-representation analysis)
# ###################################################################
#
# ORA asks a simpler question than GSEA: "Among our significant genes,
# are the genes of some pathway found more often than chance would
# predict?" The paper did this with Reactome pathways (its Fig. 2E).

# The pathway file to use. It must be in this folder.
ora_gene_sets = "reactome_human.gmt"

# Only test pathways with 2 to 2,000 of our genes in them. These are
# the default limits of ShinyGO, the tool the paper used.
ora_min_size = 2
ora_max_size = 2000

# How many pathways to show in the ORA bar plot (like the paper's
# Fig. 2E, which shows 10).
ora_top = 10


# ###################################################################
# PART 7: Comparing with the published paper
# ###################################################################

# The folder with small tables copied from the paper's figures.
# compare_to_paper.py and plot_paper_figures.py check our results
# against them.
#
# If the folder doesn't exist (for example, when you analyze a
# different dataset), those comparison steps are skipped and
# everything else still runs.
paper_reference_folder = "paper_reference"


# ###################################################################
# PART 8: Redrawing the paper's figures (plot_paper_figures.py)
# ###################################################################

# The paper's Fig. 2F is a GSEA with GO Molecular Function terms,
# while our main GSEA (PART 5) uses GO Biological Process. To redraw
# Fig. 2F, run_paper_gsea.py runs a SECOND GSEA with this collection.
# GSEApy downloads it from Enrichr, so you need internet.
paper_gsea_gene_sets = "GO_Molecular_Function_2025"

# Size limits for that second GSEA.
#
# The paper's Fig. 2F includes very large terms such as "RNA binding",
# which has more than 1,000 genes. Our main GSEA stops at 500 genes
# (gsea_max_size above), which would leave those terms out. So this
# second GSEA allows terms with up to 3,000 genes.
paper_gsea_min_size = 15
paper_gsea_max_size = 3000
