"""
===================================================================
compare_to_paper.py
===================================================================

WHAT THIS FILE DOES
    It checks our results against the published paper, automatically:

        Farhangdoost et al. (2025) Molecular Psychiatry 30:4543-4557
        https://doi.org/10.1038/s41380-025-03040-x

    It reads small tables copied from the paper's Figure 2 (in the
    paper_reference folder) and answers four questions:

      1. Did we find a similar NUMBER of significant genes?    (Fig. 2B)
      2. Do the paper's top 60 genes change in the SAME
         DIRECTION in our results, and by a similar amount?    (Fig. 2D)
      3. Did our Reactome ORA find the paper's top pathways?   (Fig. 2E)
      4. Did our GSEA find the paper's GO terms, going the
         same way?                                             (Fig. 2F)

WHERE IT FITS
    Step 11. run_analysis.py calls it near the end, if the
    paper_reference folder exists. You can also run it on its own, on
    a results folder you already have:

        python3 compare_to_paper.py akap11_go_bp_results

INPUTS
    output            the results folder (it reads 02_, 09_, 10_ and,
                      if present, 12_fig2F_ files from it)
    reference_folder  the paper_reference folder

OUTPUT FILES
    11_comparison_with_paper.md      a short report you can read
    11_paper_gene_comparison.tsv     question 2, gene by gene
    11_paper_term_comparison.tsv     questions 3 and 4, term by term
    11_paper_vs_our_fold_change.png  question 2 as a scatter plot

HOW THE "VERDICTS" ARE DECIDED (rules of thumb)
    question 1: totals within 20% of each other -> "similar"
    question 2: at least 80% same direction     -> "agree"
    question 3: at least 70% of the paper's pathways
                significant in ours             -> "agree"
    question 4: ribosome/protein-synthesis terms mostly down
                in the knockout                 -> "agree"
                (and, if tested, at least 80% of the paper's GO
                terms going the same way        -> "agree")

PYTHON IDEAS USED IN THIS FILE
    Path(...)        a file or folder path; "/" joins a folder and a
                     file name, e.g.  output / "02_results.tsv"
    np.nan           NumPy's "missing value" marker
    "x" in text      True if the text contains "x"
    ", ".join(list)  glue a list of texts into one, with ", " between
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
from settings import fdr_cutoff    # 0.05


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

    # Remove the database ID at the end.
    # text.split(x) cuts the text at x and gives a list of the pieces;
    # [0] takes the first piece, i.e. everything BEFORE x.
    #   "structural constituent of ribosome (go:0003735)"
    #   -> "structural constituent of ribosome "
    if "(go:" in name:
        name = name.split("(go:")[0]
    if " r-hsa-" in name:
        name = name.split(" r-hsa-")[0]

    # The paper's GO names start with "go ". Remove those 3 characters.
    # name[3:] means "the text from character 3 to the end".
    if name.startswith("go "):
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
    title_lines.append("# Comparison with Farhangdoost et al. (2025)")
    title_lines.append("")
    title_lines.append("Paper: Molecular Psychiatry 30:4543-4557, "
                       "https://doi.org/10.1038/s41380-025-03040-x")
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
        output / "02_differential_expression.tsv",
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

    # The paper's numbers. The file has two columns, "measure" and
    # "value", one row per number:
    #   significant_genes      849
    #   lower_in_comparison    421
    #   higher_in_comparison   428
    paper_summary = pd.read_csv(reference_folder / "paper_summary.tsv", sep="\t")

    # Turn the table into a dictionary, so we can look a number up by
    # its name: {"significant_genes": 849, "lower_in_comparison": 421, ...}
    paper_numbers = {}
    for row_number, row in paper_summary.iterrows():
        paper_numbers[row["measure"]] = int(row["value"])

    paper_total = paper_numbers["significant_genes"]
    paper_lower = paper_numbers["lower_in_comparison"]
    paper_higher = paper_numbers["higher_in_comparison"]

    # Our numbers. Comparing the direction column with a label gives
    # True/False per gene; .sum() counts the Trues.
    is_lower = results["direction"] == "Lower in comparison"
    our_lower = is_lower.sum()

    is_higher = results["direction"] == "Higher in comparison"
    our_higher = is_higher.sum()

    our_total = our_lower + our_higher

    # How far apart are the two totals, as a percentage of the paper's?
    # abs() removes a minus sign, so the order of the subtraction
    # doesn't matter.  |848 - 849| / 849 x 100 = 0.1%
    difference_percent = abs(our_total - paper_total) / paper_total * 100

    # Rule of thumb: within 20% counts as "similar".
    if difference_percent <= 20:
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
    body_lines.append(f"The totals differ by {difference_percent:.1f}%.")
    body_lines.append("")

    summary_rows.append([
        "Number of significant genes",
        f"paper {paper_total}, ours {our_total} ({difference_percent:.1f}% apart)",
        verdict,
    ])

    # ===============================================================
    # QUESTION 2: Do the paper's top genes agree with ours? (Fig. 2D)
    # ===============================================================

    # The paper's top 30 down and top 30 up genes, with their log2
    # fold change and adjusted p-value.
    paper_genes = pd.read_csv(reference_folder / "paper_degs.tsv", sep="\t")

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
            paper_direction = "Higher in comparison"
        else:
            paper_direction = "Lower in comparison"

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
            # annotation (GENCODE v47 vs the paper's Ensembl 104), or
            # it was removed by our low-count filter.
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
    gene_table.to_csv(output / "11_paper_gene_comparison.tsv", sep="\t", index=False)

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

    # The correlation of the paper's fold changes with ours, for the
    # genes we found: 1 = perfect agreement, 0 = no relationship.
    # np.corrcoef() returns a small 2 x 2 grid of correlations:
    #   [[ paper vs paper,  paper vs ours ],
    #    [ ours vs paper,   ours vs ours  ]]
    # The number we want is in row 0, column 1.
    # A correlation needs at least 3 points to mean anything.
    if len(found_genes) >= 3:
        correlation_grid = np.corrcoef(
            found_genes["paper_log2FoldChange"],
            found_genes["our_log2FoldChange"],
        )
        correlation = correlation_grid[0, 1]
    else:
        correlation = np.nan

    # The share of found genes that change in the same direction.
    if number_found > 0:
        percent_same_direction = number_same_direction / number_found * 100
    else:
        percent_same_direction = 0

    # Rule of thumb: at least 80% in the same direction = agreement.
    if percent_same_direction >= 80:
        verdict = "agree"
    else:
        verdict = "disagree"

    # The names of the paper's genes we didn't find.
    not_found = []
    for row_number, row in gene_table.iterrows():
        if row["found_in_our_results"] == False:
            not_found.append(row["gene_name"])

    body_lines.append("## 2. The paper's top 60 genes (paper Fig. 2D)")
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
    body_lines.append("Gene-by-gene details: `11_paper_gene_comparison.tsv`. "
                      "Plot: `11_paper_vs_our_fold_change.png`.")
    body_lines.append("")

    summary_rows.append([
        "Paper's top 60 genes",
        f"{number_same_direction} of {number_found} found go the same way; "
        f"{number_significant} significant in ours; r = {correlation:.2f}",
        verdict,
    ])

    # Draw the scatter plot (function at the bottom of this file).
    plot_gene_comparison(found_genes, correlation, output)

    # ===============================================================
    # QUESTION 3: Reactome pathways                    (paper Fig. 2E)
    # ===============================================================

    # The paper's top 10 Reactome pathways (one column, "term").
    paper_ora = pd.read_csv(reference_folder / "paper_reactome_ora.tsv", sep="\t")

    ora_file = output / "10_reactome_ora.tsv"

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
                    "analysis": "Reactome ORA (Fig. 2E)",
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
                    "analysis": "Reactome ORA (Fig. 2E)",
                    "paper_term": term,
                    "our_term": "not tested",
                    "our_rank": np.nan,
                    "our_fdr": np.nan,
                    "paper_direction": "",
                    "our_direction": "",
                    "agrees": False,
                })

        number_paper_ora = len(paper_ora)

        # Rule of thumb: at least 70% of the paper's pathways
        # significant in ours = agreement.
        if number_ora_significant >= number_paper_ora * 0.7:
            verdict = "agree"
        else:
            verdict = "partly agree"

        body_lines.append("## 3. Reactome pathways (paper Fig. 2E)")
        body_lines.append("")
        body_lines.append(f"Of the paper's top {number_paper_ora} pathways, "
                          f"{number_ora_found} were tested in our ORA and "
                          f"{number_ora_significant} have FDR < {fdr_cutoff} in ours.")
        body_lines.append("")

        summary_rows.append([
            "Reactome ORA pathways",
            f"{number_ora_significant} of the paper's {number_paper_ora} "
            f"top pathways are significant in ours",
            verdict,
        ])
    else:
        body_lines.append("## 3. Reactome pathways (paper Fig. 2E)")
        body_lines.append("")
        body_lines.append("Skipped: 10_reactome_ora.tsv was not found.")
        body_lines.append("")

    # ===============================================================
    # QUESTION 4: GSEA GO terms                        (paper Fig. 2F)
    # ===============================================================

    # Our main GSEA results (GO Biological Process by default).
    gsea = pd.read_csv(output / "09_gsea_pathways.tsv", sep="\t")

    # ---------------------------------------------------------------
    # Part A: the paper's 20 GO Molecular Function terms, one by one.
    # ---------------------------------------------------------------

    # Which of our GSEA tables should we compare them with?
    # If the paper-style GSEA (GO Molecular Function, from
    # run_paper_gsea.py) exists, use that: it's the same GO collection
    # as the paper's. Otherwise use our main GSEA.
    paper_style_file = output / "12_fig2F_gsea_go_mf.tsv"
    if paper_style_file.exists():
        lookup_table = pd.read_csv(paper_style_file, sep="\t")
    else:
        lookup_table = gsea.copy()

    # Simplified names, so the paper's spelling can be matched.
    simple_names = []
    for term in lookup_table["Term"]:
        simple_names.append(simplify_name(term))
    lookup_table["simple_name"] = simple_names

    lookup_table = lookup_table.drop_duplicates("simple_name")
    gsea_lookup = lookup_table.set_index("simple_name")

    # The paper's 20 terms, each with its direction.
    paper_go = pd.read_csv(reference_folder / "paper_gsea_go_mf.tsv", sep="\t")

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
                our_direction = "Higher in comparison"
            else:
                our_direction = "Lower in comparison"

            same_direction = our_direction == paper_direction
            if same_direction:
                number_go_same_direction = number_go_same_direction + 1

            term_rows.append({
                "analysis": "GSEA GO Molecular Function (Fig. 2F)",
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
                "analysis": "GSEA GO Molecular Function (Fig. 2F)",
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

    # The paper's main GSEA result: genes for ribosomes and protein
    # synthesis go DOWN in the knockout. We check this in our MAIN
    # GSEA, whichever GO collection it used.
    #
    # GO Biological Process names these terms with words like
    # "ribosome" and "translation"; GO Molecular Function uses words
    # like "rRNA binding" and "tRNA binding", so we look for all four.
    # "ribosom" (without the final e) also matches "ribosomal".
    ribosome_down = []       # the significant terms going DOWN
    number_ribosome_up = 0   # how many significant terms go UP

    for row_number, row in gsea.iterrows():
        name = row["Term"].lower()

        is_ribosome_term = False
        if "ribosom" in name:
            is_ribosome_term = True
        if "translation" in name:
            is_ribosome_term = True
        if "rrna" in name:
            is_ribosome_term = True
        if "trna" in name:
            is_ribosome_term = True

        # Only count significant terms.
        if is_ribosome_term and row["FDR q-val"] < fdr_cutoff:
            if row["NES"] < 0:
                ribosome_down.append(row)
            elif row["NES"] > 0:
                number_ribosome_up = number_ribosome_up + 1

    number_ribosome_down = len(ribosome_down)

    # The finding is reproduced if some of these terms go down, and
    # more go down than up.
    if number_ribosome_down > 0 and number_ribosome_down > number_ribosome_up:
        ribosome_verdict = "agree"
    else:
        ribosome_verdict = "disagree"

    body_lines.append("## 4. GSEA (paper Fig. 2F)")
    body_lines.append("")
    body_lines.append(f"The paper's main GSEA finding is that ribosome and protein-synthesis "
                      f"terms go down in the knockout. In our GSEA, "
                      f"{number_ribosome_down} such terms (ribosome, translation, rRNA, tRNA) are "
                      f"significantly lower and {number_ribosome_up} are significantly higher.")
    body_lines.append("")
    body_lines.append("Top three of them in our results:")
    body_lines.append("")

    # List the first three of the down terms. The GSEA table is sorted
    # by FDR, so these are the three most significant.
    number_listed = 0
    for row in ribosome_down:
        if number_listed == 3:
            break   # "break" leaves the loop early
        number_listed = number_listed + 1

        # GSEApy reports an FDR of exactly 0 when it's smaller than the
        # 1,000 permutations can measure, so we label it an estimate.
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
        body_lines.append("None of the paper's 20 GO Molecular Function terms were in our "
                          "GSEA table, because our GSEA used a different collection. "
                          "To compare term by term, set `gsea_gene_sets = "
                          "\"GO_Molecular_Function_2025\"` in settings.py.")
    else:
        body_lines.append(f"Of the paper's 20 GO Molecular Function terms, "
                          f"{number_go_found} were tested in our GSEA and "
                          f"{number_go_same_direction} change in the same direction.")
    body_lines.append("")

    summary_rows.append([
        "GSEA: ribosome/protein-synthesis terms lower in knockout",
        f"{number_ribosome_down} significant terms lower, {number_ribosome_up} higher",
        ribosome_verdict,
    ])

    # The term-by-term check only gets a summary row if any terms were
    # found. Rule of thumb: at least 80% the same way = agreement.
    if number_go_found > 0:
        if number_go_same_direction >= number_go_found * 0.8:
            verdict = "agree"
        else:
            verdict = "partly agree"

        summary_rows.append([
            "GSEA: the paper's 20 GO Molecular Function terms",
            f"{number_go_same_direction} of {number_go_found} tested go the same way",
            verdict,
        ])

    # Save the term-by-term table for questions 3 and 4.
    term_table = pd.DataFrame(term_rows)
    term_table.to_csv(output / "11_paper_term_comparison.tsv", sep="\t", index=False)

    # ===============================================================
    # Why our results can differ from the paper's
    # ===============================================================

    body_lines.append("## Why our numbers are not identical to the paper's")
    body_lines.append("")
    body_lines.append("| Step | Paper | This project |")
    body_lines.append("|---|---|---|")
    body_lines.append("| Trimming | Trimmomatic | Trim Galore (nf-core/rnaseq) |")
    body_lines.append("| Alignment | STAR | STAR |")
    body_lines.append("| Counting | HTSeq, Ensembl 104 | Salmon + tximport, GENCODE v47 |")
    body_lines.append("| Gene filter | mean count >= 20 (23,590 genes) | "
                      ">= 10 counts in >= 3 samples |")
    body_lines.append("| Statistics | limma-voom, robust fit | PyDESeq2 |")
    body_lines.append("| GSEA | fgsea, GO Molecular Function | GSEApy, setting in settings.py |")
    body_lines.append("| ORA | ShinyGO, curated Reactome | hypergeometric test, Reactome (.gmt) |")
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
    with open(output / "11_comparison_with_paper.md", "w") as file:
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
# The scatter plot for question 2
# ###################################################################

def plot_gene_comparison(found_genes, correlation, output):
    """Scatter plot: the paper's fold change vs ours, one dot per gene.

    If the two analyses agreed perfectly, every dot would sit on the
    dashed diagonal line.
    """

    # Nothing to draw if no genes were found.
    if len(found_genes) == 0:
        return

    # A square figure, 6 x 6 inches.
    fig, ax = plt.subplots(figsize=(6, 6))

    # x = the paper's fold change, y = ours.
    x = found_genes["paper_log2FoldChange"]
    y = found_genes["our_log2FoldChange"]

    ax.scatter(x, y, color="slategray", s=30)

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
    ax.plot([lowest, highest], [lowest, highest], linestyle="--", color="gray", linewidth=0.8)

    # Light lines at zero on both axes: dots in the top-right and
    # bottom-left quarters change in the same direction in both
    # analyses.
    ax.axhline(0, color="lightgray", linewidth=0.8)
    ax.axvline(0, color="lightgray", linewidth=0.8)

    ax.set_xlabel(r"Paper $\log_2$ fold change")
    ax.set_ylabel(r"Our $\log_2$ fold change")
    ax.set_title(f"Paper's top 60 genes (r = {correlation:.2f})")

    # Move overlapping labels apart; the dots stay where they are.
    adjust_text(
        labels,
        x=x,
        y=y,
        ax=ax,
        arrowprops={"arrowstyle": "-", "color": "#65717a", "linewidth": 0.5},
        iter_lim=300,
    )

    save_plot(fig, output / "11_paper_vs_our_fold_change.png")


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
