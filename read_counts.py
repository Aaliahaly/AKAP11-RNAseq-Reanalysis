"""
===================================================================
read_counts.py
===================================================================

WHAT THIS FILE DOES
    It reads the two input tables, checks that they look right,
    removes genes with very few reads, and prepares readable gene
    names for the plots.

    The two input tables are:

    1. The sample sheet (samples.tsv): which sample is in which group.

         sample              condition   geo_sample   sra_run
         AKAP11_KO_clone16   KO          GSM8200086   SRR28641569
         ...
         WT_rep3             WT          GSM8200091   SRR28641564

    2. The count table (AKAP11_final_gene_counts_length_scaled.tsv):
       one row per gene, one column per sample. Each number is
       (roughly) how many sequencing reads came from that gene.

         gene_id              gene_name   AKAP11_KO_clone16   ...
         ENSG00000000003.16   TSPAN6      9775.336            ...

WHERE IT FITS
    Step 1 of the analysis. run_analysis.py calls read_counts() first,
    and every later step uses the three tables it returns.

INPUTS
    count_file        path to the count table
    sample_file       path to the sample sheet
    reference_group   the baseline group, here "WT"
    comparison_group  the group compared with it, here "KO"

OUTPUTS (returned to run_analysis.py, not saved as files)
    filtered    counts for the genes we keep (genes x samples)
    metadata    one row per sample: its group and two QC numbers
    annotation  one row per gene: its name and its plot label

PYTHON IDEAS USED IN THIS FILE
    DataFrame     a pandas table, with labeled rows and columns
    index         the row labels of a DataFrame
    True/False filtering
                  table[ list_of_True_False ] keeps only the rows
                  where the list says True
    raise ValueError("...")
                  stop the program with an error message. We do this
                  when an input file is wrong, so the problem is
                  caught at the start instead of causing strange
                  results later.
    list, .append()
                  a list is an ordered collection, like [1, 2, 3];
                  .append(x) adds x to the end of it
    dictionary    a lookup table of key -> value pairs,
                  like {"HOXA2": 1, "NDN": 2}
===================================================================
"""

# numpy ("np") is the basic library for numbers and arrays in Python.
# We use it here to check that every count is a normal number.
import numpy as np

# pandas ("pd") is the library for tables (DataFrames).
# Almost everything in this project uses it.
import pandas as pd

# Read four settings from settings.py (see that file for explanations).
from settings import gene_id_column     # "gene_id"
from settings import gene_name_column   # "gene_name"
from settings import min_count          # 10
from settings import min_samples        # None -> use the smaller group size


def read_counts(count_file, sample_file, reference_group, comparison_group):

    # ===============================================================
    # STEP 1: Read the sample sheet
    # ===============================================================

    # pd.read_csv() reads a text table into a pandas DataFrame.
    # Despite its name, it reads any separator, not only commas.
    #
    #   sep="\t"              -> the columns are separated by tabs
    #                            ("\t" is how Python writes a tab)
    #   dtype=str             -> read every value as text, so a
    #                            sample called "1" stays the text "1"
    #   keep_default_na=False -> leave empty cells as "" (empty text).
    #                            Without this, pandas turns empty cells
    #                            (and even the word "NA") into NaN,
    #                            its "missing value" marker.
    metadata = pd.read_csv(
        sample_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    # metadata.columns holds the column names. list() turns them into
    # a plain Python list, for example:
    #   ["sample", "condition", "geo_sample", "sra_run"]
    sample_sheet_columns = list(metadata.columns)

    # The sample sheet MUST have a "sample" column and a "condition"
    # column. If either is missing, stop with a clear message.
    if "sample" not in sample_sheet_columns:
        raise ValueError("The sample sheet needs sample and condition columns.")

    if "condition" not in sample_sheet_columns:
        raise ValueError("The sample sheet needs sample and condition columns.")

    # Every sample needs its own, non-empty name: we use the name to
    # find that sample's column in the count table.
    #
    # We go through the names one by one, and keep a list of the names
    # we have already seen, so we can spot repeats.
    names_seen_so_far = []

    # This loop runs once per row of the sample sheet.
    # Each time, sample_name holds the next name, e.g. "WT_rep1".
    for sample_name in metadata["sample"]:

        # An empty name is not allowed.
        if sample_name == "":
            raise ValueError("Sample names must be unique and nonempty.")

        # A name we have already seen is a repeat, which is not allowed.
        if sample_name in names_seen_so_far:
            raise ValueError("Sample names must be unique and nonempty.")

        # Remember this name, so a later repeat of it gets caught.
        names_seen_so_far.append(sample_name)

    # ---------------------------------------------------------------
    # Keep only the rows for the two groups we want to compare.
    # (A sample sheet could contain other groups too; we ignore them.)
    # ---------------------------------------------------------------

    # Comparing a whole column with one value gives one True/False per
    # row. With reference_group = "WT":
    #
    #   condition   == "WT" ?
    #   KO          False
    #   KO          False
    #   KO          False
    #   WT          True
    #   WT          True
    #   WT          True
    is_reference_row = metadata["condition"] == reference_group

    # Putting that True/False list inside [ ] keeps only the True rows,
    # so "reference" is a smaller table with just the WT samples.
    reference = metadata[is_reference_row]

    # The same for the comparison group (KO).
    is_comparison_row = metadata["condition"] == comparison_group
    comparison = metadata[is_comparison_row]

    # len() of a table gives its number of rows, so these are the
    # number of samples in each group (3 and 3 for our data).
    number_of_reference_samples = len(reference)
    number_of_comparison_samples = len(comparison)

    # Find the size of the smaller group. We use it twice: to check
    # there are enough samples, and (in STEP 3) for the gene filter.
    if number_of_reference_samples < number_of_comparison_samples:
        smaller_group = number_of_reference_samples
    else:
        smaller_group = number_of_comparison_samples

    # Each group needs at least 2 samples (replicates). With only one
    # sample, there is no way to measure how much samples naturally
    # vary, and the statistics in fit_model.py can't work.
    # (This also catches a typo in a group name: a misspelled group
    # matches 0 rows.)
    if smaller_group < 2:
        raise ValueError("Each selected group needs at least two samples.")

    # ---------------------------------------------------------------
    # Put the two groups back together into one table.
    # ---------------------------------------------------------------

    # pd.concat() stacks tables on top of each other. The reference
    # group (WT) goes first, then the comparison group (KO). This
    # order decides the order of the samples in every later table.
    metadata = pd.concat([reference, comparison])

    # Use the sample names as the row labels (the "index").
    # After this, we can look a sample up by its name, for example
    # metadata.loc["WT_rep1"].
    metadata = metadata.set_index("sample")

    # Keep only the "condition" column; we don't need the GEO and SRA
    # IDs for the analysis.
    # Single brackets, metadata["condition"], would give one column on
    # its own. Double brackets, metadata[["condition"]], keep it as a
    # table, which is what PyDESeq2 expects later.
    metadata = metadata[["condition"]]

    # A plain list of the sample names, in order:
    #   ["WT_rep1", "WT_rep2", "WT_rep3",
    #    "AKAP11_KO_clone16", "AKAP11_KO_clone20", "AKAP11_KO_clone21"]
    samples = list(metadata.index)

    # ---------------------------------------------------------------
    # Decide how many samples must pass the count cutoff (STEP 3).
    # ---------------------------------------------------------------

    # settings.py sets min_samples = None, which means "use the size of
    # the smaller group" (3 here). "is None" is how Python checks for
    # the special value None.
    if min_samples is None:
        required_samples = smaller_group
    else:
        required_samples = min_samples

    total_samples = len(samples)

    # A sensible value is between 1 and the number of samples (6).
    # "or" means: stop if EITHER of the two problems is true.
    if required_samples < 1 or required_samples > total_samples:
        raise ValueError("min_samples must be between 1 and the number of selected samples.")

    # ===============================================================
    # STEP 2: Read the count table
    # ===============================================================

    # Read the count table. dtype={gene_id_column: str} says: read the
    # gene ID column as text. Gene IDs are names, not numbers, and
    # reading them as text guarantees they are kept exactly as written.
    table = pd.read_csv(
        count_file,
        sep="\t",
        dtype={gene_id_column: str},
    )

    # Build a list of every column we expect in the count table:
    # the gene ID column, the gene name column, and one column per
    # sample.
    required_columns = []
    required_columns.append(gene_id_column)
    required_columns.append(gene_name_column)
    for sample_name in samples:
        required_columns.append(sample_name)

    # Check the expected columns one at a time. A missing sample column
    # usually means a sample is spelled differently in the two files.
    count_table_columns = list(table.columns)

    for column in required_columns:
        if column not in count_table_columns:
            # "+" joins two pieces of text, so the message names the
            # missing column, e.g. "Missing count table column: WT_rep4"
            raise ValueError("Missing count table column: " + column)

    # Make sure no name appears twice in our list. For example, a
    # sample must not be called "gene_id"; otherwise one column would
    # be used for two different jobs.
    columns_seen_so_far = []

    for column in required_columns:
        if column in columns_seen_so_far:
            raise ValueError("Gene ID, gene symbol, and sample columns must be separate.")
        columns_seen_so_far.append(column)

    # Use the gene IDs as the row labels (the index), so we can look a
    # gene up by its ID, e.g. table.loc["ENSG00000023516.9"].
    table = table.set_index(gene_id_column)

    # Every gene must have an ID: .hasnans is True if any are missing.
    if table.index.hasnans:
        raise ValueError("Gene IDs must be unique and nonempty.")

    # Every gene ID must be different: .has_duplicates is True if any
    # ID appears twice.
    if table.index.has_duplicates:
        raise ValueError("Gene IDs must be unique and nonempty.")

    # ---------------------------------------------------------------
    # Take the sample columns and make sure they hold numbers.
    # ---------------------------------------------------------------

    # Select only the sample columns, in the same order as the sample
    # sheet. table[samples] picks several columns by their names.
    # .copy() makes a separate table, so changing it later can't
    # accidentally change "table" too.
    counts = table[samples].copy()

    # Convert each sample column to numbers, one column at a time.
    # pd.to_numeric() stops with an error if a cell holds something
    # that isn't a number (like the word "abc").
    for sample_name in samples:
        counts[sample_name] = pd.to_numeric(counts[sample_name])

    # Check that every value is a normal number.
    #   counts.to_numpy()  -> the numbers as a plain grid (no labels)
    #   np.isfinite(...)   -> True for normal numbers; False for a
    #                         missing value (NaN) or infinity
    #   .all()             -> True only if EVERY value is True
    all_values_are_finite = np.isfinite(counts.to_numpy()).all()

    # Check for negative numbers: a count can't be below 0.
    #   counts < 0         -> a True/False table, True where negative
    #   first .any()       -> one True/False per column:
    #                         "does this column have a negative?"
    #   second .any()      -> one True/False overall:
    #                         "does ANY column have a negative?"
    is_negative = counts < 0
    any_negative_counts = is_negative.any().any()

    # "not" flips True to False and back, so this reads:
    # "if it's NOT true that all values are finite, stop".
    if not all_values_are_finite:
        raise ValueError("Sample columns must contain finite, nonnegative counts.")

    if any_negative_counts:
        raise ValueError("Sample columns must contain finite, nonnegative counts.")

    # ===============================================================
    # STEP 3: Round the counts and remove low-count genes
    # ===============================================================

    # The counts in this file are estimates from Salmon, so they have
    # decimals (like 9775.336). PyDESeq2 needs whole numbers.
    #
    # First round every value to the nearest whole number
    # (9775.336 -> 9775.0)...
    rounded_counts = counts.round()

    # ...then store them as integers, Python's whole-number type
    # (9775.0 -> 9775).
    counts = rounded_counts.astype(int)

    # For every gene in every sample, ask: "at least 10 reads?"
    # The result is a table of True/False with the same shape as
    # counts. Example for one gene:
    #
    #   counts:   WT_rep1 12   WT_rep2 15   WT_rep3 9   KO... 0  0  1
    #   >= 10 ?   True         True         False       False ...
    has_enough_counts = counts >= min_count

    # Add up each row (axis=1 means "across the columns of a row").
    # Python counts True as 1 and False as 0, so the sum is the number
    # of samples that passed. For the example above: 2.
    number_of_passing_samples = has_enough_counts.sum(axis=1)

    # One True/False per gene: did enough samples (3) pass?
    keep_this_gene = number_of_passing_samples >= required_samples

    # .loc[ True/False list ] keeps only the rows (genes) marked True.
    # For our data this keeps 26,288 of the 78,932 genes.
    filtered = counts.loc[keep_this_gene]

    # The model needs at least 2 genes. With real data this is never a
    # problem, but it catches a wrong setting (like min_count = 1000000).
    if len(filtered) < 2:
        raise ValueError("Filtering left fewer than two genes or an empty sample.")

    # Every sample must still have some reads after filtering.
    # .sum() with no axis adds DOWN each column: one total per sample.
    filtered_totals = filtered.sum()

    for sample_name in samples:
        if filtered_totals[sample_name] == 0:
            raise ValueError("Filtering left fewer than two genes or an empty sample.")

    # ---------------------------------------------------------------
    # Save two quality-control (QC) numbers for each sample.
    # These use "counts" (ALL genes, before filtering), not "filtered".
    # They are written to 01_sample_qc.tsv later by fit_model.py.
    # ---------------------------------------------------------------

    # Total reads (counts) in each sample.
    # Assigning a new column name creates a new column. pandas matches
    # the values to the right rows by sample name.
    total_counts_per_sample = counts.sum()
    metadata["total_estimated_counts_rounded"] = total_counts_per_sample

    # Number of genes with at least one read in each sample.
    #   counts > 0  -> True where the gene was seen in that sample
    #   .sum()      -> counts the Trues down each column
    gene_is_detected = counts > 0
    detected_genes_per_sample = gene_is_detected.sum()
    metadata["genes_detected_before_filter"] = detected_genes_per_sample

    # ===============================================================
    # STEP 4: Make readable gene labels for the plots
    # ===============================================================

    # We build two lists, with one entry per gene, in the same order
    # as the count table:
    #   clean_names -> the gene name, tidied up ("" if it has none)
    #   labels      -> the text to show on plots (always filled in)
    #
    # list(table.index) gives every gene ID; list(table[...]) gives
    # every gene name. Position 0 in one list matches position 0 in
    # the other.
    gene_ids = list(table.index)
    raw_names = list(table[gene_name_column])

    clean_names = []

    for name in raw_names:

        # A missing gene name is stored as NaN ("not a number"),
        # pandas' marker for an empty cell. pd.isna() is True for it.
        # Replace it with "" (empty text).
        if pd.isna(name):
            name = ""

        # Make sure it is text (str), then remove any spaces at the
        # start and end (" HOXA2 " -> "HOXA2").
        name = str(name)
        name = name.strip()

        clean_names.append(name)

    # Build the plot labels.
    # If a gene has a name, use it. If not, use its gene ID, so every
    # gene still gets a label.
    labels = []

    # range(len(gene_ids)) counts 0, 1, 2, ... up to the number of genes
    # minus 1. "position" is used to take the matching item from each
    # list.
    for position in range(len(gene_ids)):
        gene_id = gene_ids[position]
        name = clean_names[position]

        # "!=" means "is not equal to".
        if name != "":
            labels.append(name)
        else:
            labels.append(gene_id)

    # A few different genes share the same name. We want every label
    # to be unique (two genes called "PDE11A" in a plot would be
    # confusing), so first count how often each label appears.
    #
    # times_label_appears is a dictionary: label -> how many times seen.
    # For example: {"TSPAN6": 1, "PDE11A": 2, ...}
    times_label_appears = {}

    for label in labels:
        if label in times_label_appears:
            # Seen before: add 1 to its count.
            times_label_appears[label] = times_label_appears[label] + 1
        else:
            # First time: start its count at 1.
            times_label_appears[label] = 1

    # Now, for every label that appears more than once, add the gene ID
    # in brackets so it becomes unique, for example:
    #   "PDE11A"  ->  "PDE11A (ENSG00000128655.17)"
    for position in range(len(labels)):
        label = labels[position]

        if times_label_appears[label] > 1:
            gene_id = gene_ids[position]
            labels[position] = label + " (" + gene_id + ")"

    # Put the names and labels into one table, with the gene IDs as
    # row labels. The { } part says: column "gene_name" holds the list
    # clean_names, and column "label" holds the list labels.
    annotation = pd.DataFrame(
        {
            "gene_name": clean_names,
            "label": labels,
        },
        index=table.index,
    )

    # ---------------------------------------------------------------
    # Print a short summary in Terminal.
    # ---------------------------------------------------------------

    # An f-string (text starting with f") fills in the values inside
    # { }. The ":," part adds thousands separators, so 26288 prints
    # as 26,288. This prints:
    #   Kept 26,288 of 78,932 genes after filtering.
    number_kept = len(filtered)
    number_before = len(counts)
    print(f"Kept {number_kept:,} of {number_before:,} genes after filtering.")

    # Send three tables back to run_analysis.py. "return" with several
    # values separated by commas sends them all back together, and
    # run_analysis.py receives them in the same order:
    #   filtered   -> counts for the genes we kept (used by the model)
    #   metadata   -> each sample's group and QC numbers
    #   annotation -> each gene's name and plot label
    return filtered, metadata, annotation
