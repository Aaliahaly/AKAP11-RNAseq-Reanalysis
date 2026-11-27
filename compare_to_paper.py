"""
===================================================================
compare_to_paper.py
===================================================================

WHAT THIS FILE DOES
    It checks our results against the published paper, automatically.
    (Which paper, and all its numbers, are read from the files in the
    paper_reference folder; for this project that is
    Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557.)

    It reads small tables copied from the paper's Figure 2 (in the
    paper_reference folder) and answers four questions:

      1. Did we find a similar NUMBER of significant genes?    (Fig. 2B)
      2. Do the paper's top genes change in the SAME
         DIRECTION in our results, and by a similar amount?    (Fig. 2D)
      3. Did our ORA find the paper's top pathways?            (Fig. 2E)
      4. Did our GSEA find the paper's GO terms, going the
         same way?                                             (Fig. 2F)

WHERE IT FITS
    Step 12. run_analysis.py calls it near the end, if the
    paper_reference folder exists. You can also run it on its own, on
    a results folder you already have:

        python3 compare_to_paper.py akap11_go_bp_results

INPUTS
    output            the results folder (it reads the 01_, 02_, 09_,
                      10_ and, if present, 12_fig2F_ files from it)
    reference_folder  the paper_reference folder

OUTPUT FILES (names set in settings.py, PART 10)
    11_comparison_with_paper.md      a short report you can read
    11_paper_gene_comparison.tsv     question 2, gene by gene
    11_paper_term_comparison.tsv     questions 3 and 4, term by term
    11_paper_vs_our_fold_change.png  question 2 as a scatter plot

HOW THE "VERDICTS" ARE DECIDED (rules of thumb)
    The cutoffs are in settings.py (PART 7), each with its reason:
    question 1: totals within similar_total_percent (20%)  -> "similar"
    question 2: at least agree_direction_percent (80%) of
                the paper's genes going the same way       -> "agree"
    question 3: at least agree_ora_percent (70%) of the
                paper's pathways significant in ours       -> "agree"
    question 4: the paper's headline GSEA terms (for this
                paper: ribosome / protein synthesis) mostly
                going the paper's way                      -> "agree"
                (and, if tested, at least 80% of the paper's GO
                terms going the same way                   -> "agree")

WHY NO NUMBER FROM THE PAPER IS TYPED IN THIS FILE
    Every number from the paper comes from a file in paper_reference,
    and every sentence in the report that mentions a number fills it
    in from the data. If the data or a setting changes, the report
    changes with it; it can never show an old number.

PYTHON IDEAS USED IN THIS FILE
    Path(...)        a file or folder path; "/" joins a folder and a
                     file name, e.g.  output / "02_results.tsv"
    np.nan           NumPy's "missing value" marker
    "x" in text      True if the text contains "x"
    ", ".join(list)  glue a list of texts into one, with ", " between
    text.format(**dictionary)
                     fill in the {names} in a text from a dictionary:
                     "{a} and {b}".format(**{"a": 1, "b": 2}) -> "1 and 2"
    if __name__ == "__main__":
                     code under this line runs only when you start
                     this file directly from Terminal
===================================================================
"""

# sys lets us read the words typed in Terminal (sys.argv).
import sys

# Path makes file and folder paths that work on every computer.
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

from save_plot import save_plot

# Cutoffs and labels shared with the rest of the analysis.
from settings import fdr_cutoff                  # 0.05
from settings import min_count                   # 10
from settings import min_samples                 # None -> smaller group size
from settings import condition_column            # "condition"
from settings import higher_label                # "Higher in comparison"
from settings import lower_label                 # "Lower in comparison"
from settings import gsea_gene_sets              # our main GSEA collection
from settings import paper_gsea_gene_sets        # the paper-style GSEA collection
from settings import ora_gene_sets               # our ORA pathway file
from settings import ora_collection_label        # "Reactome"
from settings import term_id_markers             # [" R-HSA-", "(GO:"]

# The verdict rules (PART 7 of settings.py).
from settings import similar_total_percent       # 20
from settings import agree_direction_percent     # 80
from settings import agree_ora_percent           # 70
from settings import report_examples             # 3
from settings import min_points_for_correlation  # 3

# Colors for the scatter plot.
from settings import neutral_color               # "slategray"
from settings import cutoff_line_color           # gray
from settings import label_line_color            # dark gray

# The files in the paper_reference folder (PART 9 of settings.py).
from settings import paper_info_file
from settings import paper_summary_file
from settings import paper_methods_file
from settings import paper_genes_file
from settings import paper_ora_file
from settings import paper_gsea_file

# Our own results files, read and written here (PART 10).
from settings import sample_qc_file
from settings import differential_expression_file
from settings import gsea_results_file
from settings import ora_results_file
from settings import paper_gsea_results_file
from settings import comparison_report_file
from settings import gene_comparison_file
from settings import term_comparison_file
from settings import fold_change_comparison_plot_file


# ###################################################################
# Helper functions for reading the paper_reference files
# (plot_paper_figures.py uses them too, so they live here once)
# ###################################################################

def read_paper_info(reference_folder):
    """Read paper_info.tsv into a dictionary of field -> text.

    The file has two columns, "field" and "value", for example:
        short_citation   Farhangdoost et al. (2025)
        fig1b_gene       AKAP11
    and we return {"short_citation": "Farhangdoost et al. (2025)", ...}
    """

    # dtype=str reads every value as text, so nothing is turned into a
    # number by accident.
    table = pd.read_csv(Path(reference_folder) / paper_info_file, sep="\t", dtype=str)

    info = {}
    for row_number, row in table.iterrows():
        info[row["field"]] = row["value"]
    return info


def read_paper_numbers(reference_folder):
    """Read paper_summary.tsv into a dictionary of measure -> number.

    For example {"significant_genes": 849.0, "fdr_cutoff": 0.05, ...}.
    Every value is read as a decimal number (float); a count like 849
    becomes 849.0, and int() turns it back into 849 where needed.
    """

    table = pd.read_csv(Path(reference_folder) / paper_summary_file, sep="\t")

    numbers = {}
    for row_number, row in table.iterrows():
        numbers[row["measure"]] = float(row["value"])
    return numbers


def paper_citation(info):
    """Build the citation line used at the top of both reports."""

    # For example:
    #   "Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557,
    #    https://doi.org/10.1038/s41380-025-03040-x"
    return info["short_citation"] + " " + info["full_citation"] + ", " + info["doi_link"]


def read_methods(output, reference_folder):
    """Read paper_methods.tsv and fill in our side from settings.py.

    The "this_project" column can contain names in curly brackets,
    such as {min_count}. These are filled in here from settings.py and
    from our results, so the table always shows the values that were
    really used. Returns a dictionary:
        step -> {"paper": ..., "this_project": ..., "affects_gene_list": "yes" or "no"}
    """

    # ---------------------------------------------------------------
    # The values that can be filled in.
    # ---------------------------------------------------------------

    # How many genes went into the model = rows of the results table.
    results = pd.read_csv(Path(output) / differential_expression_file, sep="\t", index_col=0)
    genes_tested = len(results)

    # How many samples a gene had to pass in. If min_samples is None
    # in settings.py, read_counts.py used the size of the smaller
    # group, so we work that out again from the sample table.
    if min_samples is None:
        samples = pd.read_csv(Path(output) / sample_qc_file, sep="\t", index_col=0)
        # .value_counts() counts how many samples are in each group,
        # e.g. WT 3, KO 3; .min() takes the smaller count.
        samples_needed = samples[condition_column].value_counts().min()
    else:
        samples_needed = min_samples

    values = {
        "min_count": min_count,
        "min_samples": samples_needed,
        "genes_tested": genes_tested,
        "gsea_gene_sets": gsea_gene_sets,
        "ora_gene_sets": ora_gene_sets,
    }

    # ---------------------------------------------------------------
    # Read the table and fill in the {names}.
    # ---------------------------------------------------------------

    table = pd.read_csv(Path(reference_folder) / paper_methods_file, sep="\t", dtype=str)

    methods = {}
    for row_number, row in table.iterrows():
        # "...{min_count}...".format(**values) swaps each {name} for
        # its value from the dictionary. Text without { } is unchanged.
        ours = row["this_project"].format(**values)
        methods[row["step"]] = {
            "paper": row["paper"],
            "this_project": ours,
            "affects_gene_list": row["affects_gene_list"],
        }
    return methods


# ###################################################################
# Helper function: make term names comparable
# ###################################################################

def simplify_name(name):
    """Make a pathway or GO term name easy to compare.

    Different databases write the same term in slightly different ways:

        paper: "GO STRUCTURAL CONSTITUENT OF RIBOSOME"
        ours:  "Structural Constituent of Ribosome (GO:0003735)"

        paper: "G alpha (q) signalling events"
        ours:  "G Alpha (Q) Signaling Events R-HSA-416476"

    This function turns both versions of each into the same simple
    text, so they can be matched:

        "structural constituent of ribosome"
        "g alpha q signaling events"
    """

    # Make sure it's text, then make it all lowercase.
    name = str(name).lower()

    # Remove the database ID at the end. The ID markers (" R-HSA-",
    # "(GO:") are listed in settings.py. We lowercase each marker too,
    # because the name is now lowercase.
    # text.split(x) cuts the text at x and gives a list of the pieces;
    # [0] takes the first piece, i.e. everything BEFORE x.
    #   "structural constituent of ribosome (go:0003735)"
    #   -> "structural constituent of ribosome "
    for marker in term_id_markers:
        name = name.split(marker.lower())[0]

    # Some collections (MSigDB, used by many papers) start GO names
    # with "GO_" or "GO ". Remove those 3 characters: name[3:] means
    # "the text from character 3 to the end".
    if name.startswith("go ") or name.startswith("go_"):
        name = name[3:]

    # Make British and American spellings the same.
    name = name.replace("signalling", "signaling")

    # Replace every character that isn't a letter or a number with a
    # space. This removes brackets, dashes, slashes and commas.
    # .isalnum() is True for letters (a-z) and digits (0-9).
    cleaned = ""
    for character in name:
        if character.isalnum():
            cleaned = cleaned + character
        else:
            cleaned = cleaned + " "

    # .split() with nothing in the brackets cuts at every run of spaces
    # and drops empty pieces; " ".join() glues the words back with
    # single spaces. Together they turn "g alpha  q  signaling" into
    # "g alpha q signaling".
    words = cleaned.split()
    return " ".join(words)


# ###################################################################
# The main comparison
# ###################################################################

def compare_to_paper(output, reference_folder):

    # Turn the two folder names into Path objects, so "/" can join
    # them with file names. (They may arrive as plain text when this
    # file is run on its own from Terminal.)
    output = Path(output)
    reference_folder = Path(reference_folder)

    # Facts about the paper (citation, headline finding) and its
    # numbers, from the paper_reference folder.
    info = read_paper_info(reference_folder)
    paper_numbers = read_paper_numbers(reference_folder)

    # The report is built as lists of text lines, and written to a
    # Markdown file at the end. It has three parts:
    #
    #   title_lines    -> the title and the paper's reference
    #   summary_lines  -> a summary table (made at the end, from
    #                     summary_rows)
    #   body_lines     -> one section per question
    #
    # summary_rows holds one entry per check:
    #   [what was checked, the result in words, the verdict]
    title_lines = []
    title_lines.append("# Comparison with " + info["short_citation"])
    title_lines.append("")
    title_lines.append("Paper: " + paper_citation(info))
    title_lines.append("")

    body_lines = []
    summary_rows = []

    # term_rows collects one row per pathway or GO term for questions
    # 3 and 4. It's saved as 11_paper_term_comparison.tsv at the end.
    term_rows = []

    # ===============================================================
    # STEP 0: Read our results
    # ===============================================================

    # Read the full results table from fit_model.py.
    #   index_col=0          -> use the first column (gene IDs) as
    #                           row labels
    #   keep_default_na=False, na_values=[""]
    #                        -> treat ONLY truly empty cells as missing,
    #                           so a gene name like "NA" stays a name
    results = pd.read_csv(
        output / differential_expression_file,
        sep="\t",
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )

    # Genes without a name are read as missing (NaN). fillna("")
    # replaces those with empty text, as in the rest of the project.
    results["gene_name"] = results["gene_name"].fillna("")

    # ===============================================================
    # QUESTION 1: How many significant genes?          (paper Fig. 2B)
    # ===============================================================

    # The paper's numbers, from paper_summary.tsv. int() turns 849.0
    # into 849, because these are counts of genes.
    paper_total = int(paper_numbers["significant_genes"])
    paper_lower = int(paper_numbers["lower_in_comparison"])
    paper_higher = int(paper_numbers["higher_in_comparison"])

    # Our numbers. Comparing the direction column with a label gives
    # True/False per gene; .sum() counts the Trues.
    is_lower = results["direction"] == lower_label
    our_lower = is_lower.sum()

    is_higher = results["direction"] == higher_label
    our_higher = is_higher.sum()

    our_total = our_lower + our_higher

    # How far apart are the two totals, as a percentage of the paper's?
    # abs() removes a minus sign, so the order of the subtraction
    # doesn't matter.  |848 - 849| / 849 x 100 = 0.1%
    difference_percent = abs(our_total - paper_total) / paper_total * 100

    # Rule of thumb (settings.py): within similar_total_percent (20%)
    # counts as "similar".
    if difference_percent <= similar_total_percent:
        verdict = "similar"
    else:
        verdict = "different"

    # The report section for question 1: a small Markdown table.
    # In Markdown, "|" separates table columns, and "|---|" is the line
    # under the header row.
    body_lines.append("## 1. Number of significant genes (paper Fig. 2B)")
    body_lines.append("")
    body_lines.append("| | Paper | Ours |")
    body_lines.append("|---|---|---|")
    body_lines.append(f"| Total | {paper_total} | {our_total} |")
    body_lines.append(f"| Lower in comparison group | {paper_lower} | {our_lower} |")
    body_lines.append(f"| Higher in comparison group | {paper_higher} | {our_higher} |")
    body_lines.append("")
    body_lines.append(f"The totals differ by {difference_percent:.1f}% "
                      f"(\"similar\" means at most {similar_total_percent}% apart).")
    body_lines.append("")

    summary_rows.append([
        "Number of significant genes",
        f"paper {paper_total}, ours {our_total} ({difference_percent:.1f}% apart)",
        verdict,
    ])

    # ===============================================================
    # QUESTION 2: Do the paper's top genes agree with ours? (Fig. 2D)
    # ===============================================================

    # The paper's top down and top up genes, with their log2 fold
    # change and adjusted p-value.
    paper_genes = pd.read_csv(reference_folder / paper_genes_file, sep="\t")

    # ---------------------------------------------------------------
    # Prepare our results for looking genes up by name.
    # ---------------------------------------------------------------

    # Keep only genes that have a name.
    has_name = results["gene_name"] != ""
    our_genes = results[has_name].copy()

    # We'll compare names in UPPERCASE, so a difference in letter case
    # can't cause a missed match.
    our_genes["symbol_upper"] = our_genes["gene_name"].str.upper()

    # A few names appear under more than one gene ID. For each name,
    # keep the most highly expressed copy: sort by baseMean, highest
    # first, then keep the first row for each name.
    our_genes = our_genes.sort_values("baseMean", ascending=False, kind="stable")
    our_genes = our_genes.drop_duplicates("symbol_upper")

    # Use the uppercase names as row labels, so we can do
    # our_genes.loc["NDN"].
    our_genes = our_genes.set_index("symbol_upper")

    # ---------------------------------------------------------------
    # Compare the paper's genes one at a time.
    # ---------------------------------------------------------------

    # One dictionary per paper gene; becomes a table at the end.
    gene_rows = []

    for row_number, paper_gene in paper_genes.iterrows():

        symbol = paper_gene["gene_name"]
        paper_fold_change = paper_gene["log2FoldChange"]

        # Which way did the gene change in the paper?
        if paper_fold_change > 0:
            paper_direction = higher_label
        else:
            paper_direction = lower_label

        # Look the gene up in our results.
        if symbol.upper() in our_genes.index:
            ours = our_genes.loc[symbol.upper()]
            found = True
            our_fold_change = ours["log2FoldChange"]
            our_padj = ours["padj"]
            our_direction = ours["direction"]

            # Same direction = both positive, or both negative.
            # (our_fold_change > 0) is True or False, and so is
            # (paper_fold_change > 0). If the two match, the
            # directions match.
            same_direction = (our_fold_change > 0) == (paper_fold_change > 0)

            # Is it also SIGNIFICANT in our results, in the same
            # direction?
            significant_same_way = our_direction == paper_direction
        else:
            # Not found. The gene may have a different name in our
            # annotation, or it was removed by our low-count filter.
            # np.nan marks the missing numbers.
            found = False
            our_fold_change = np.nan
            our_padj = np.nan
            our_direction = "not in our results"
            same_direction = False
            significant_same_way = False

        # Store this gene's comparison as one row.
        gene_rows.append({
            "gene_name": symbol,
            "paper_log2FoldChange": paper_fold_change,
            "paper_padj": paper_gene["padj"],
            "our_log2FoldChange": our_fold_change,
            "our_padj": our_padj,
            "our_direction": our_direction,
            "found_in_our_results": found,
            "same_direction": same_direction,
            "significant_in_ours_same_direction": significant_same_way,
        })

    # Turn the rows into a table and save it.
    gene_table = pd.DataFrame(gene_rows)
    gene_table.to_csv(output / gene_comparison_file, sep="\t", index=False)

    # ---------------------------------------------------------------
    # Summarize question 2.
    # ---------------------------------------------------------------

    # .sum() of a True/False column counts the Trues.
    number_of_paper_genes = len(gene_table)
    number_found = gene_table["found_in_our_results"].sum()
    number_same_direction = gene_table["same_direction"].sum()
    number_significant = gene_table["significant_in_ours_same_direction"].sum()

    # Keep only the genes we found (True in the found column).
    found_genes = gene_table[gene_table["found_in_our_results"]]

    # The correlation of the paper's fold changes with ours (see
    # fold_change_correlation() below).
    correlation = fold_change_correlation(found_genes)

    # The share of found genes that change in the same direction.
    if number_found > 0:
        percent_same_direction = number_same_direction / number_found * 100
    else:
        percent_same_direction = 0

    # Rule of thumb (settings.py): at least agree_direction_percent
    # (80%) in the same direction = agreement.
    if percent_same_direction >= agree_direction_percent:
        verdict = "agree"
    else:
        verdict = "disagree"

    # The names of the paper's genes we didn't find.
    not_found = []
    for row_number, row in gene_table.iterrows():
        if row["found_in_our_results"] == False:
            not_found.append(row["gene_name"])

    body_lines.append(f"## 2. The paper's top {number_of_paper_genes} genes (paper Fig. 2D)")
    body_lines.append("")
    body_lines.append(f"- Found in our results: {number_found} of {number_of_paper_genes}")
    body_lines.append(f"- Same direction of change: {number_same_direction} of "
                      f"{number_found} found ({percent_same_direction:.0f}%)")
    body_lines.append(f"- Also significant in our results, same direction: "
                      f"{number_significant} of {number_of_paper_genes}")
    body_lines.append(f"- Correlation of log2 fold changes: {correlation:.2f}")
    if len(not_found) > 0:
        not_found_text = ", ".join(not_found)
        body_lines.append(f"- Not found (renamed in our annotation, or removed by "
                          f"our low-count filter): {not_found_text}")
    body_lines.append("")
    body_lines.append(f"Gene-by-gene details: `{gene_comparison_file}`. "
                      f"Plot: `{fold_change_comparison_plot_file}`.")
    body_lines.append("")

    summary_rows.append([
        f"Paper's top {number_of_paper_genes} genes",
        f"{number_same_direction} of {number_found} found go the same way; "
        f"{number_significant} significant in ours; r = {correlation:.2f}",
        verdict,
    ])

    # Draw the scatter plot (function at the bottom of this file).
    plot_gene_comparison(found_genes, correlation, number_of_paper_genes, output)

    # ===============================================================
    # QUESTION 3: ORA pathways                         (paper Fig. 2E)
    # ===============================================================

    # The paper's top pathways (one column, "term").
    paper_ora = pd.read_csv(reference_folder / paper_ora_file, sep="\t")

    # The label for these rows in the term table, e.g.
    # "Reactome ORA (Fig. 2E)". plot_paper_figures.py looks for rows
    # starting with ora_collection_label.
    ora_analysis_name = ora_collection_label + " ORA (Fig. 2E)"

    ora_file = output / ora_results_file
    section_title = f"## 3. {ora_collection_label} pathways (paper Fig. 2E)"

    # .exists() is True if the file is there. It might not be, if ORA
    # was skipped.
    if ora_file.exists():

        our_ora = pd.read_csv(ora_file, sep="\t")

        # Number our rows 1, 2, 3, ... (the table is already sorted by
        # FDR), so we can report each pathway's rank.
        # range(1, n + 1) counts 1 up to n.
        our_ora["rank"] = range(1, len(our_ora) + 1)

        # Make a simplified name for every one of our pathways, so it
        # can be matched with the paper's spelling.
        simple_names = []
        for term in our_ora["Term"]:
            simple_names.append(simplify_name(term))
        our_ora["simple_name"] = simple_names

        # Keep one row per simplified name, and use it as row labels.
        our_ora = our_ora.drop_duplicates("simple_name")
        our_ora = our_ora.set_index("simple_name")

        number_ora_found = 0
        number_ora_significant = 0

        # Check the paper's pathways one at a time.
        for term in paper_ora["term"]:
            key = simplify_name(term)

            if key in our_ora.index:
                ours = our_ora.loc[key]
                number_ora_found = number_ora_found + 1

                # Is it significant in our ORA?
                is_significant = ours["fdr"] < fdr_cutoff
                if is_significant:
                    number_ora_significant = number_ora_significant + 1

                term_rows.append({
                    "analysis": ora_analysis_name,
                    "paper_term": term,
                    "our_term": ours["Term"],
                    "our_rank": int(ours["rank"]),
                    "our_fdr": ours["fdr"],
                    "paper_direction": "",
                    "our_direction": "",
                    "agrees": is_significant,
                })
            else:
                term_rows.append({
                    "analysis": ora_analysis_name,
                    "paper_term": term,
                    "our_term": "not tested",
                    "our_rank": np.nan,
                    "our_fdr": np.nan,
                    "paper_direction": "",
                    "our_direction": "",
                    "agrees": False,
                })

        number_paper_ora = len(paper_ora)

        # Rule of thumb (settings.py): at least agree_ora_percent (70%)
        # of the paper's pathways significant in ours = agreement.
        # "/ 100" turns the percentage into a fraction (70 -> 0.7).
        if number_ora_significant >= number_paper_ora * agree_ora_percent / 100:
            verdict = "agree"
        else:
            verdict = "partly agree"

        body_lines.append(section_title)
        body_lines.append("")
        body_lines.append(f"Of the paper's top {number_paper_ora} pathways, "
                          f"{number_ora_found} were tested in our ORA and "
                          f"{number_ora_significant} have FDR < {fdr_cutoff} in ours.")
        body_lines.append("")

        summary_rows.append([
            f"{ora_collection_label} ORA pathways",
            f"{number_ora_significant} of the paper's {number_paper_ora} "
            f"top pathways are significant in ours",
            verdict,
        ])
    else:
        body_lines.append(section_title)
        body_lines.append("")
        body_lines.append(f"Skipped: {ora_results_file} was not found.")
        body_lines.append("")

    # ===============================================================
    # QUESTION 4: GSEA GO terms                        (paper Fig. 2F)
    # ===============================================================

    # Our main GSEA results (GO Biological Process by default).
    gsea = pd.read_csv(output / gsea_results_file, sep="\t")

    # ---------------------------------------------------------------
    # Part A: the paper's GO terms, one by one.
    # ---------------------------------------------------------------

    # Which of our GSEA tables should we compare them with?
    # If the paper-style GSEA (from run_paper_gsea.py) exists, use
    # that: it's the same GO collection as the paper's. Otherwise use
    # our main GSEA.
    paper_style_file = output / paper_gsea_results_file
    if paper_style_file.exists():
        lookup_table = pd.read_csv(paper_style_file, sep="\t")
        compared_collection = paper_gsea_gene_sets
    else:
        lookup_table = gsea.copy()
        compared_collection = gsea_gene_sets

    # Simplified names, so the paper's spelling can be matched.
    simple_names = []
    for term in lookup_table["Term"]:
        simple_names.append(simplify_name(term))
    lookup_table["simple_name"] = simple_names

    lookup_table = lookup_table.drop_duplicates("simple_name")
    gsea_lookup = lookup_table.set_index("simple_name")

    # The paper's terms, each with its direction.
    paper_go = pd.read_csv(reference_folder / paper_gsea_file, sep="\t")
    number_paper_go = len(paper_go)

    # The label for these rows in the term table.
    gsea_analysis_name = "GSEA " + compared_collection + " (Fig. 2F)"

    number_go_found = 0
    number_go_same_direction = 0

    for row_number, paper_term in paper_go.iterrows():
        key = simplify_name(paper_term["term"])
        paper_direction = paper_term["direction"]

        if key in gsea_lookup.index:
            ours = gsea_lookup.loc[key]
            number_go_found = number_go_found + 1

            # A positive NES means the term leans toward KO (higher in
            # the comparison group).
            if ours["NES"] > 0:
                our_direction = higher_label
            else:
                our_direction = lower_label

            same_direction = our_direction == paper_direction
            if same_direction:
                number_go_same_direction = number_go_same_direction + 1

            term_rows.append({
                "analysis": gsea_analysis_name,
                "paper_term": paper_term["term"],
                "our_term": ours["Term"],
                "our_rank": np.nan,
                "our_fdr": ours["FDR q-val"],
                "paper_direction": paper_direction,
                "our_direction": our_direction,
                "agrees": same_direction,
            })
        else:
            term_rows.append({
                "analysis": gsea_analysis_name,
                "paper_term": paper_term["term"],
                "our_term": "not tested",
                "our_rank": np.nan,
                "our_fdr": np.nan,
                "paper_direction": paper_direction,
                "our_direction": "",
                "agrees": False,
            })

    # ---------------------------------------------------------------
    # Part B: the paper's headline GSEA finding.
    # ---------------------------------------------------------------

    # paper_info.tsv describes the paper's main GSEA result in three
    # fields. For this paper:
    #   gsea_headline           "ribosome and protein-synthesis terms
    #                            go down in the knockout"
    #   gsea_headline_keywords  "ribosom;translation;rrna;trna"
    #   gsea_headline_direction "Lower in comparison"
    #
    # We check this in our MAIN GSEA, whichever GO collection it used.
    # A term counts if its name contains any keyword. (Short keywords
    # catch more spellings: "ribosom" matches "ribosome" and
    # "ribosomal".)
    headline_text = info["gsea_headline"]
    headline_direction = info["gsea_headline_direction"]

    # .split(";") cuts the text at every ";" -> a list of keywords.
    keywords = info["gsea_headline_keywords"].split(";")

    # The NES sign that matches the paper's direction: negative (-1)
    # for "lower", positive (+1) for "higher".
    if headline_direction == lower_label:
        expected_sign = -1
    else:
        expected_sign = 1

    headline_terms_expected_way = []     # significant terms going the paper's way
    number_headline_other_way = 0        # significant terms going the other way

    for row_number, row in gsea.iterrows():
        name = row["Term"].lower()

        # Does the term's name contain any of the keywords?
        is_headline_term = False
        for keyword in keywords:
            if keyword.lower() in name:
                is_headline_term = True

        # Only count significant terms.
        if is_headline_term and row["FDR q-val"] < fdr_cutoff:
            # np.sign() gives -1 for a negative number, +1 for a
            # positive one, and 0 for exactly zero.
            if np.sign(row["NES"]) == expected_sign:
                headline_terms_expected_way.append(row)
            elif np.sign(row["NES"]) == -expected_sign:
                number_headline_other_way = number_headline_other_way + 1

    number_headline_expected_way = len(headline_terms_expected_way)

    # The finding is reproduced if some of these terms go the paper's
    # way, and more go that way than the other way.
    if (number_headline_expected_way > 0
            and number_headline_expected_way > number_headline_other_way):
        headline_verdict = "agree"
    else:
        headline_verdict = "disagree"

    # Words for the two directions, for the sentence below.
    if expected_sign == -1:
        expected_word = "lower"
        other_word = "higher"
    else:
        expected_word = "higher"
        other_word = "lower"

    keyword_text = ", ".join(keywords)

    body_lines.append("## 4. GSEA (paper Fig. 2F)")
    body_lines.append("")
    body_lines.append(f"The paper's main GSEA finding is that {headline_text}. In our GSEA "
                      f"({gsea_gene_sets}), {number_headline_expected_way} matching terms "
                      f"(names containing: {keyword_text}) are significantly {expected_word} "
                      f"and {number_headline_other_way} are significantly {other_word}.")
    body_lines.append("")
    body_lines.append(f"Top {report_examples} of them in our results:")
    body_lines.append("")

    # List the first few terms going the paper's way (report_examples
    # = 3 in settings.py). The GSEA table is sorted by FDR, so these
    # are the most significant.
    number_listed = 0
    for row in headline_terms_expected_way:
        if number_listed == report_examples:
            break   # "break" leaves the loop early
        number_listed = number_listed + 1

        # GSEApy reports an FDR of exactly 0 when it's smaller than its
        # permutations can measure, so we label it an estimate.
        if row["FDR q-val"] == 0:
            fdr_text = "0 (estimated)"
        else:
            fdr_text = f"{row['FDR q-val']:.3g}"

        term = row["Term"]
        nes = row["NES"]
        body_lines.append(f"- {term}: NES {nes:.2f}, FDR {fdr_text}")
    body_lines.append("")

    # A sentence about part A.
    if number_go_found == 0:
        body_lines.append(f"None of the paper's {number_paper_go} GO terms were in our "
                          f"GSEA table ({compared_collection}), probably because it used a "
                          f"different collection. To compare term by term, run the analysis "
                          f"with internet access (so the {paper_gsea_gene_sets} GSEA can run), "
                          f"or set `gsea_gene_sets = \"{paper_gsea_gene_sets}\"` in settings.py.")
    else:
        body_lines.append(f"Of the paper's {number_paper_go} GO terms, "
                          f"{number_go_found} were tested in our GSEA ({compared_collection}) and "
                          f"{number_go_same_direction} change in the same direction.")
    body_lines.append("")

    summary_rows.append([
        "GSEA: " + headline_text,
        f"{number_headline_expected_way} significant terms {expected_word}, "
        f"{number_headline_other_way} {other_word}",
        headline_verdict,
    ])

    # The term-by-term check only gets a summary row if any terms were
    # found. Rule of thumb: at least agree_direction_percent (80%) the
    # same way = agreement.
    if number_go_found > 0:
        if number_go_same_direction >= number_go_found * agree_direction_percent / 100:
            verdict = "agree"
        else:
            verdict = "partly agree"

        summary_rows.append([
            f"GSEA: the paper's {number_paper_go} GO terms",
            f"{number_go_same_direction} of {number_go_found} tested go the same way",
            verdict,
        ])

    # Save the term-by-term table for questions 3 and 4.
    term_table = pd.DataFrame(term_rows)
    term_table.to_csv(output / term_comparison_file, sep="\t", index=False)

    # ===============================================================
    # Why our results can differ from the paper's
    # ===============================================================

    # The methods table comes from paper_methods.tsv, with our side
    # filled in from settings.py (see read_methods() above).
    methods = read_methods(output, reference_folder)

    body_lines.append("## Why our numbers are not identical to the paper's")
    body_lines.append("")
    body_lines.append("| Step | Paper | This project |")
    body_lines.append("|---|---|---|")
    # Looping over a dictionary gives its keys, here the step names,
    # in the order of the file.
    for step in methods:
        paper_method = methods[step]["paper"]
        our_method = methods[step]["this_project"]
        body_lines.append(f"| {step} | {paper_method} | {our_method} |")
    body_lines.append("")

    # ===============================================================
    # Build the summary table, write the report, print the summary
    # ===============================================================

    summary_lines = []
    summary_lines.append("## Summary")
    summary_lines.append("")
    summary_lines.append("| Check | Result | Verdict |")
    summary_lines.append("|---|---|---|")
    for summary_row in summary_rows:
        check = summary_row[0]
        result = summary_row[1]
        verdict = summary_row[2]
        summary_lines.append(f"| {check} | {result} | {verdict} |")
    summary_lines.append("")

    # The report, in reading order: title, then summary, then the
    # sections. "+" joins lists one after the other.
    report = title_lines + summary_lines + body_lines

    # Write every line to the file. "w" opens it for writing (replacing
    # any old version); "\n" ends each line.
    # encoding="utf-8" makes sure symbols such as "≥" are saved
    # correctly on every computer (Windows uses another default).
    with open(output / comparison_report_file, "w", encoding="utf-8") as file:
        for line in report:
            file.write(line + "\n")

    # Print the summary in Terminal too, e.g.
    #   [agree] Paper's top 60 genes: 55 of 55 found go the same way; ...
    print("Comparison with the paper:")
    for summary_row in summary_rows:
        check = summary_row[0]
        result = summary_row[1]
        verdict = summary_row[2]
        print(f"  [{verdict}] {check}: {result}")


# ###################################################################
# Helper: the correlation of the paper's fold changes with ours
# (plot_paper_figures.py uses it too)
# ###################################################################

def fold_change_correlation(found_genes):
    """Pearson correlation of the paper's log2 fold changes with ours.

    1 = perfect agreement, 0 = no relationship. found_genes must have
    the columns "paper_log2FoldChange" and "our_log2FoldChange".
    """

    # A correlation needs at least min_points_for_correlation (3)
    # points to mean anything (see settings.py). With fewer, return
    # np.nan ("missing").
    if len(found_genes) < min_points_for_correlation:
        return np.nan

    # np.corrcoef() returns a small 2 x 2 grid of correlations:
    #   [[ paper vs paper,  paper vs ours ],
    #    [ ours vs paper,   ours vs ours  ]]
    # The number we want is in row 0, column 1.
    correlation_grid = np.corrcoef(
        found_genes["paper_log2FoldChange"],
        found_genes["our_log2FoldChange"],
    )
    return correlation_grid[0, 1]


# ###################################################################
# The scatter plot for question 2
# ###################################################################

def plot_gene_comparison(found_genes, correlation, number_of_paper_genes, output):
    """Scatter plot: the paper's fold change vs ours, one dot per gene.

    If the two analyses agreed perfectly, every dot would sit on the
    dashed diagonal line.
    """

    # Nothing to draw if no genes were found.
    if len(found_genes) == 0:
        return

    # A square figure, 6 x 6 inches (visual only).
    fig, ax = plt.subplots(figsize=(6, 6))

    # x = the paper's fold change, y = ours.
    x = found_genes["paper_log2FoldChange"]
    y = found_genes["our_log2FoldChange"]

    ax.scatter(x, y, color=neutral_color, s=30)

    # Write each gene's name next to its dot, and keep the labels in a
    # list so adjust_text() can move them apart afterwards.
    labels = []
    for row_number, gene in found_genes.iterrows():
        label = ax.text(gene["paper_log2FoldChange"], gene["our_log2FoldChange"],
                        gene["gene_name"], fontsize=6)
        labels.append(label)

    # The dashed diagonal: where the two fold changes would be equal.
    # It runs from the smallest value on either axis to the largest.
    lowest = min(x.min(), y.min())
    highest = max(x.max(), y.max())
    ax.plot([lowest, highest], [lowest, highest], linestyle="--",
            color=cutoff_line_color, linewidth=0.8)

    # Light lines at zero (= no change) on both axes: dots in the
    # top-right and bottom-left quarters change in the same direction
    # in both analyses.
    ax.axhline(0, color="lightgray", linewidth=0.8)
    ax.axvline(0, color="lightgray", linewidth=0.8)

    ax.set_xlabel(r"Paper $\log_2$ fold change")
    ax.set_ylabel(r"Our $\log_2$ fold change")
    ax.set_title(f"Paper's top {number_of_paper_genes} genes (r = {correlation:.2f})")

    # Move overlapping labels apart; the dots stay where they are.
    adjust_text(
        labels,
        x=x,
        y=y,
        ax=ax,
        arrowprops={"arrowstyle": "-", "color": label_line_color, "linewidth": 0.5},
        iter_lim=300,
    )

    save_plot(fig, output / fold_change_comparison_plot_file)


# ###################################################################
# Run this file on its own
# ###################################################################

# Python sets __name__ to "__main__" only for the file you start from
# Terminal. So this part runs when you type
#     python3 compare_to_paper.py akap11_go_bp_results
# but NOT when run_analysis.py imports this file.
if __name__ == "__main__":

    # Save plots to files instead of opening windows.
    import matplotlib
    matplotlib.use("Agg")

    # Use the same plot style as run_analysis.py.
    import seaborn as sns
    sns.set_theme(style="whitegrid")

    from settings import paper_reference_folder

    # sys.argv is the list of words typed in Terminal:
    #   ["compare_to_paper.py", "akap11_go_bp_results"]
    # We need exactly 2 words (the script and one folder).
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 compare_to_paper.py results_folder")

    results_folder = sys.argv[1]
    compare_to_paper(results_folder, paper_reference_folder)
