# AKAP11 RNA Sequencing Workflow on Expanse

**Part 1 of the project: from FASTQ files to gene counts.** Part 2, the Python analysis of the counts, is described in the [main README](../README.md) and summarized in section 14.

*Bulk RNA-seq reanalysis of GSE263847 with nf-core/rnaseq on SDSC Expanse  - November 2026*

I analyzed six paired-end RNA-seq samples from GEO series GSE263847 (BioProject PRJNA1099597) to compare three heterozygous AKAP11 knockout clones with three wild-type replicates. All six are iPSC-derived neurons from the SBP009 cell line. I used Nextflow 25.04.8 and nf-core/rnaseq 3.19.0 on SDSC Expanse, with STAR for alignment, Salmon for quantification, and the full set of QC tools, including RSeQC read duplication and biotype QC.

This document walks through the whole analysis, from logging in to Expanse to opening the final QC report on my Mac, and then summarizes the Python analysis of the gene counts and how my results compare with the paper (section 14). For each step, I explain what the commands do, why the step is needed, and what could go wrong, so that someone else (or future me) can repeat the analysis from scratch.

### Project Details

| Item | Details |
|---|---|
| Dataset | GEO GSE263847 / BioProject PRJNA1099597 |
| Samples | 6 iPSC-derived neuron samples (SBP009 line): 3 heterozygous AKAP11 KO clones (16, 20, 21) and 3 wild-type replicates |
| Sequencing | Illumina NovaSeq 6000, paired-end, about 145 to 294 million read pairs per sample |
| Reference | GRCh38 primary assembly + GENCODE v47 annotation |
| Pipeline | nf-core/rnaseq 3.19.0 (STAR + Salmon), run with Nextflow 25.04.8 and SingularityPRO 4.1.2 |
| Computing | SDSC Expanse, Slurm shared partition, allocation gtw111 |
| Main outputs | Gene count table (salmon.merged.gene_counts_length_scaled.tsv) and MultiQC report |

## 0 The big picture: what I'm doing and why

### The biological question

AKAP11 (A-kinase anchoring protein 11) is a scaffold protein that holds protein kinase A (PKA) in particular places inside the cell. Large exome-sequencing studies have linked rare loss-of-function variants in AKAP11 to bipolar disorder and schizophrenia. The study I'm reanalyzing made neurons from induced pluripotent stem cells (iPSCs) and used CRISPR to make clones that lose one copy of AKAP11 (heterozygous knockouts), which mimics what the patient variants do. My question is: which genes change their expression in neurons when one copy of AKAP11 is lost?

### What RNA-seq measures

RNA-seq measures gene expression by sequencing the RNA in a sample. The RNA is copied into DNA (cDNA), cut into short fragments, and each fragment is read by an Illumina sequencer. This dataset is paired-end, so each fragment was read from both ends, giving me two files per sample: read 1 (\_1.fastq.gz) and read 2 (\_2.fastq.gz). The basic logic: if a gene is more active, more of its RNA is in the cell, so more reads come from it. Counting reads per gene tells me how active each gene is, and comparing those counts between knockout and wild-type tells me what AKAP11 affects.

### The workflow

| Step | What I do | Main tools | What comes out |
|---|---|---|---|
| 1 | Log in and set up folders | ssh, module, lfs | A project folder on Expanse |
| 2 | Get the human reference | curl, gzip, sha256sum | Genome FASTA + gene annotation GTF |
| 3 | Find the samples and their labels | ENA API, GEO SOFT | Run table and list of SRA runs |
| 4 | Download the raw reads | SRA Toolkit, pigz, Slurm array | 12 FASTQ files (6 samples x 2 reads) |
| 5 | Write the sample sheet | Python | Checked sample sheet + metadata |
| 6 | Tell Nextflow how to use Expanse | Nextflow config, Singularity | configs/expanse.config |
| 7 | Test with small data | nf-core test profile | Proof the setup works |
| 8 | Run the real analysis | nf-core/rnaseq: STAR, Salmon, QC tools | Gene counts + QC |
| 9 | Monitor the jobs | squeue, sacct, sstat | Confirmation the run finished |
| 10 - 12 | Find the results and copy them to my Mac | find, scp, open | Counts table and MultiQC report on my laptop |
| 14 | Analyze the counts on my Mac | Python: PyDESeq2, GSEApy, SciPy | Differential expression, plots, GSEA, ORA, comparison with the paper |

### Where each command runs

There are three different computers involved, and it matters which one I'm typing into. My Mac is where I open Terminal, and where scp and open run in steps 11 and 12. The Expanse login node is where I land after ssh. It's shared by everyone, so it's only for light work like editing files, submitting jobs and checking on them. The Expanse compute nodes are where the real work happens. I never log in to them directly: I write a batch script, submit it with sbatch, and the scheduler (Slurm) runs it on a compute node when one is free.

### How to read this document

Commands are in code blocks. Lines that start with # are my comments explaining the command underneath. The boxes mean:

> [!TIP]
> **Expected outcome:** what the step produces.

> [!WARNING]
> **What could go wrong:** common problems at this step, how to recognize them, and how to fix them.

> [!NOTE]
> **Key idea:** a concept I had to learn to understand the step.

> [!IMPORTANT]
> **What I learned:** the main takeaways from each step.

### Using this guide with your own account

Everything in this guide ran under my own Expanse account. To repeat it under yours, change these values wherever they appear:

| What | In this guide | Where it appears |
|---|---|---|
| Expanse username | aaly1 | ssh and scp commands, the project path, and configs/expanse.config. The batch scripts use $USER, so they need no change. |
| Allocation (account) | gtw111 | Every #SBATCH --account line, clusterOptions in configs/expanse.config, and the project path |
| Project folder | /expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq | configs/expanse.config (TMPDIR and the Singularity bind), and remote= in steps 11 and 12 |
| Job IDs | for example 55597317 and 55746106 | sacct, squeue and log-file commands. Use the job IDs that sbatch prints for you. |

### Why nf-core/rnaseq with STAR + Salmon suits this data

Before choosing a pipeline, I wrote down what kind of data this is, because the data type decides which tools make sense:

| Property of the data | What it means for the pipeline |
|---|---|
| Bulk RNA-seq from a cell culture (not single-cell) | I need a bulk pipeline that ends in one count per gene per sample |
| Human, with a good reference genome (GRCh38) | I can align reads to the genome instead of assembling transcripts from scratch |
| Illumina short reads, paired-end | The aligner must handle read pairs |
| Reads come from spliced mRNA | The aligner must be splice-aware: one read can span two exons with an intron removed in between |
| Library strandedness not stated on GEO | The pipeline should detect it from the data |
| Very deep (145 to 294 million read pairs per sample) | The pipeline must run on a cluster and restart after failures without redoing finished work |
| Goal: compare KO vs WT gene expression | The output must be count-like values that DESeq2 or limma can use |

nf-core/rnaseq fits every row of that table:

**It's built for exactly this job.** It turns bulk, short-read FASTQ files into gene counts plus quality control. It is a community pipeline that is peer-reviewed, versioned (I used 3.19.0) and tested on every release, so other people can rerun the same version and get the same steps.

**STAR is splice-aware,** and it's the same aligner the original authors used. That keeps the alignment step directly comparable with the paper.

**Salmon handles ambiguous reads better than simple counting.** Many reads match more than one transcript or overlapping genes. The paper counted with HTSeq (-m intersection-nonempty), which throws those reads away. In star_salmon mode, Salmon takes STAR's alignments and estimates how the ambiguous reads should be shared between transcripts, so more of the data is used. tximport then sums transcripts into genes and gives the length-scaled counts that DESeq2 expects.

**Strandedness is detected automatically.** I set it to auto, so the pipeline tests a sample of reads and uses the right setting, and RSeQC infer_experiment confirms it for every sample. All six libraries came out reverse-stranded, which matches the paper's HTSeq setting (stranded = reverse). See step 12.

**Quality control is built in.** FastQC, Trim Galore, STAR alignment rates, Qualimap, RSeQC, dupRadar, featureCounts biotypes and a DESeq2 PCA are all collected into one MultiQC report, so I can check that the data is good before trusting the counts.

**It runs reliably on Expanse.** Nextflow sends each step to Slurm as its own job and runs each tool in a Singularity container, so tool versions are fixed. -resume restarts from the last finished step after a timeout, which mattered for these deep samples.

Alternatives I considered:

| Option | Why I didn't choose it |
|---|---|
| Rebuild the paper's exact chain (Trimmomatic, STAR, HTSeq) | Possible, but I would have to write and maintain every step myself, with no built-in QC or restart. nf-core keeps the same core (trim, then STAR) and adds reproducibility. |
| Salmon only, no genome alignment (--pseudo_aligner salmon) | Faster, but it gives no genome alignments (BAM files), so most of the alignment QC (RSeQC, Qualimap, dupRadar) would be missing. |
| HISAT2 instead of STAR | Uses less memory, but STAR is what the paper used, and Expanse has enough memory for it. |

How my pipeline differs from the paper's, and so where small differences in results can come from:

| Step | Paper | This project |
|---|---|---|
| Trimming | Trimmomatic (Phred &lt; 30) | Trim Galore (nf-core default) |
| Alignment | STAR, GRCh38 | STAR, GRCh38 primary assembly |
| Annotation | Ensembl 104 | GENCODE v47 (newer, so some genes have been renamed) |
| Counting | HTSeq, ambiguous reads dropped | Salmon + tximport, ambiguous reads shared out |
| Gene filter | Mean count ≥ 20 (23,590 genes) | ≥ 10 counts in ≥ 3 samples (26,288 genes) |
| Differential expression | limma-voom, robust fit | PyDESeq2 |
| Significance cutoffs | Adjusted p &lt; 0.05 and \|log2FC\| &gt; 1 | Adjusted p ≤ 0.05 and \|log2FC\| ≥ 1 |

Despite these differences, my results agree closely with the paper (section 14).

## 1 Connect to Expanse and set up my project

I opened Terminal on my Mac and logged in over SSH. Expanse asks for my password and then a six-digit TOTP code from my authenticator app. My nextflow_env conda environment turns on by itself when I log in, which is why my prompt starts with (nextflow_env).

```bash
# Log in to Expanse (password first, then the TOTP code)
ssh aaly1@login.expanse.sdsc.edu

# Make sure the tools I need are installed and see which versions I have
nextflow -version
java -version
singularity --version
command -v singularity

# Check how many compute hours (SUs) my allocations have left
module load sdsc
expanse-client user -r expanse

# See which SRA Toolkit module Expanse offers
module avail sra
```

> [!NOTE]
> **Key idea: Modules and conda environments.** Expanse has hundreds of programs installed, but most aren't switched on by default, because different users need different versions. module load &lt;name&gt; switches one on for my current session (for example module load sratoolkit/2.10.9/rn4humf), and module avail &lt;word&gt; lists the versions available. Conda is a second way to get software: I made my own environment called nextflow_env that contains Nextflow and Java, and it turns on whenever I log in.

Next, I made one project folder under my gtw111 project storage, with a subfolder for each file type, and pointed Nextflow's work and container cache folders inside it.

```bash
# Save my project path in a variable so I don't have to keep typing it
export AKAP11_ROOT=/expanse/lustre/projects/gtw111/"$USER"/akap11_rnaseq

# Make all the subfolders I'll use in one go
mkdir -p "$AKAP11_ROOT"/{metadata,results,logs,scripts,configs,tmp}
mkdir -p "$AKAP11_ROOT"/{work,singularity_cache,sra_cache,sra_tmp,fastq_sra}

# Tell Nextflow to keep its work files and containers inside this project
export NXF_WORK="$AKAP11_ROOT/work"
export NXF_SINGULARITY_CACHEDIR="$AKAP11_ROOT/singularity_cache"

cd "$AKAP11_ROOT"
pwd

# Check my storage quota (user and group) and free space on the filesystem
lfs quota -u "$USER" /expanse/lustre/projects
lfs quota -g gtw111 /expanse/lustre/projects
df -h .
```

> [!NOTE]
> **Key idea: Why project storage and not my home folder.** My home folder (/home/aaly1) is small. Project storage on /expanse/lustre/projects is a Lustre parallel filesystem made for big data that many compute nodes read at the same time. This project needs far more space than home can give: about 73 GB of compressed reads, plus the SRA archives, the unzipped reference, the STAR index and Nextflow's work folder. That's why lfs quota (the Lustre version of a disk-quota check) is part of my setup.

> [!TIP]
> **Expected outcome:** Nextflow 25.04.8, Java 17 and SingularityPRO 4.1.2 are available. The gtw111 allocation has 100,000 SUs available, and pwd prints /expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq.

> [!WARNING]
> **What could go wrong:** Expanse SSH sessions often drop ("Broken pipe" or "Connection reset by peer"), especially when idle. That's why every long step in this workflow runs as a Slurm batch job: the job keeps going on the cluster even if my laptop disconnects. After a drop, I log back in, run cd /expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq, and check my jobs with squeue. Exported variables like $AKAP11_ROOT disappear when the session ends, which is why every script below uses the full project path.

> [!IMPORTANT]
> **What I learned**
> - The login node is for setting up, submitting and checking jobs. Anything heavy goes into an sbatch script.
> - An SU (service unit) is roughly one CPU core used for one hour on Expanse, so checking my allocation tells me how much computing I can afford.
> - Making all the folders at the start keeps everything organized and avoids typos in long paths later.

## 2 Download the reference genome and annotation

> [!NOTE]
> **Key idea: FASTA vs GTF.** To know which gene a read came from, I need two reference files. The FASTA file (GRCh38.primary_assembly.genome.fa) is the DNA sequence of every human chromosome, the "map". The GTF file (gencode.v47...gtf) is the annotation: it lists where every gene, transcript and exon sits on that map. They have to agree on chromosome names and coordinates, which is why I got both from the same GENCODE release (v47). I used the primary assembly, which leaves out alternative versions (haplotypes) of some regions. Otherwise the same read could match two near-identical copies of a region, and the aligner couldn't tell where it really came from.

I saved this script as scripts/prepare_reference.sbatch.

```bash
#!/bin/bash
#SBATCH --job-name=akap11_reference
#SBATCH --account=gtw111
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --constraint=lustre
#SBATCH --output=logs/reference_%j.out
#SBATCH --error=logs/reference_%j.err

# Stop right away if any command fails
set -euo pipefail

cd "/expanse/lustre/projects/gtw111/$USER/akap11_rnaseq"

mkdir -p reference/gencode_v47
cd reference/gencode_v47

base="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47"

for filename in \
    GRCh38.primary_assembly.genome.fa.gz \
    gencode.v47.primary_assembly.annotation.gtf.gz
do
    # Download to a .part file first so a half-finished download never looks complete
    if [ ! -f "$filename" ]; then
        curl --fail --location \
            --retry 5 --retry-delay 15 \
            --connect-timeout 60 \
            --continue-at - \
            --output "$filename.part" "$base/$filename"

        gzip -t "$filename.part"
        mv "$filename.part" "$filename"
    fi

    # Make sure the gzip file isn't corrupted
    gzip -t "$filename"

    # Unzip it (again through a .part file)
    output="${filename%.gz}"
    if [ ! -f "$output" ]; then
        gzip -dc "$filename" > "$output.part"
        mv "$output.part" "$output"
    fi
done

# Record checksums so I can tell later if these files ever change
sha256sum \
    GRCh38.primary_assembly.genome.fa \
    gencode.v47.primary_assembly.annotation.gtf \
    > reference.sha256

ls -lh
echo "Reference files prepared successfully."
```

Every batch script starts with #SBATCH lines. They look like comments to bash, but Slurm reads them as the job's resource request. These are the ones I use in every script:

| #SBATCH option | What it means for my job |
|---|---|
| --account=gtw111 | Which allocation pays for the job (charged in SUs) |
| --partition=shared | Share a node with other users, so I only pay for the cores I ask for |
| --nodes=1 --ntasks=1 | One computer, one task. Expanse requires ntasks to be set |
| --cpus-per-task / --mem | How many CPU cores and how much RAM to reserve |
| --time=HH:MM:SS | Wall-clock limit. Slurm kills the job when time runs out (TIMEOUT) |
| --constraint=lustre | Only use nodes that can reach the Lustre project storage |
| --output / --error | Where stdout and stderr go. %j becomes the job ID |

```bash
# Submit the reference job (mine was job 55597317)
sbatch scripts/prepare_reference.sbatch

# Check whether it finished and how long it took
sacct -j 55597317 --format=JobID,JobName,State,Elapsed,ExitCode --allocations

# Read the end of the job log and list the files it made
tail -n 10 logs/reference_55597317.out
ls -lh reference/gencode_v47/
```

> [!TIP]
> **Expected outcome:** The job COMPLETED in 12 min 40 s. The folder holds the 3.0 GB genome FASTA, the 1.7 GB GTF, both original .gz files and reference.sha256, and the log ends with "Reference files prepared successfully."

> [!WARNING]
> **What could go wrong:** Long downloads from public servers can stop partway through. The --retry and --continue-at - options make curl try again and pick up where it stopped, instead of starting over. If the genome and GTF came from different releases or sources (for example, Ensembl chromosome names like "1" with a GENCODE GTF that uses "chr1"), the pipeline would find almost no reads in genes. Using one release for both avoids that.

> [!IMPORTANT]
> **What I learned**
> - set -euo pipefail makes a script stop at the first error, instead of carrying on and making broken files. -e stops on errors, -u stops on unset variables, and pipefail catches errors in the middle of a pipe.
> - Downloading to a .part file and renaming it only after it passes a check means a half-finished file can never be mistaken for a finished one.
> - ${filename%.gz} is bash for "filename with .gz cut off the end", so it turns the .fa.gz name into the .fa name.
> - A sha256 checksum is a fingerprint of a file. If even one letter of the file changed, the fingerprint would be completely different.

## 3 Get the study metadata and list the runs

> [!NOTE]
> **Key idea: GEO, SRA, ENA and accession numbers.** When scientists publish sequencing data, it gets IDs at several levels. GEO (GSE263847 for the study, GSM8200086 and so on for each sample) stores the description: sample names, genotypes, cell type. SRA at NCBI stores the raw reads, organised as BioProject (PRJNA1099597) → BioSample (SAMN...) → Experiment (SRX...) → Run (SRR...). ENA in Europe mirrors the same data. The SRR run accession is what I actually download, so the most important job of this step is matching each SRR to its GEO sample title, because that title is the only place that says which sample is a knockout and which is wild type.

I pulled the run table from ENA to get the six SRA run accessions and file sizes, then downloaded the GEO SOFT file to find out which run belongs to which sample.

```bash
# Download the run table for the BioProject from ENA
curl --fail --location --retry 3 \
  'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1099597&result=read_run&fields=run_accession,experiment_accession,sample_accession,library_strategy,library_layout,instrument_model,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv' \
  --output metadata/ena_runs.tsv

# Look at the table and count the runs (should be 6)
head -n 12 metadata/ena_runs.tsv
awk 'END {print NR-1}' metadata/ena_runs.tsv

# Add up the FASTQ sizes so I know how much space the reads need
awk -F '\t' 'NR>1 {
    n=split($11,sizes,";")
    for(i=1;i<=n;i++) total+=sizes[i]
}
END {
    printf "%.2f GB compressed\n",total/1000000000
}' metadata/ena_runs.tsv

# Download the GEO SOFT file, which has the sample titles and genotypes
curl --fail --location --retry 3 \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263847/soft/GSE263847_family.soft.gz' \
  --output metadata/GSE263847_family.soft.gz

# Print only the lines I care about for each sample
gzip -dc metadata/GSE263847_family.soft.gz |
awk '/^\^SAMPLE =/ ||
     /^!Sample_title =/ ||
     /^!Sample_characteristics_ch1 =/ ||
     /^!Sample_relation =/'

# Save just the six SRR accessions, one per line, for the download array in step 4
awk -F '\t' 'NR > 1 && $1 ~ /^SRR[0-9]+$/ { print $1 }' metadata/ena_runs.tsv > metadata/sra_runs.txt
cat metadata/sra_runs.txt
```

> [!TIP]
> **Expected outcome:** Six paired-end Illumina NovaSeq 6000 RNA-seq runs, about 72.90 GB of compressed FASTQ in total, and sra_runs.txt lists the six SRR accessions. Following the GEO SOFT file from sample to BioSample to SRA experiment gives this mapping:

| My label | SRA run | GEO sample | GEO title | Genotype |
|---|---|---|---|---|
| AKAP11_KO_clone16 | SRR28641569 | GSM8200086 | Het_AKAP11-KO-Clone 16 | Het AKAP11 KO |
| AKAP11_KO_clone20 | SRR28641568 | GSM8200087 | Het_AKAP11-KO-Clone 20 | Het AKAP11 KO |
| AKAP11_KO_clone21 | SRR28641567 | GSM8200088 | Het_AKAP11-KO-Clone 21 | Het AKAP11 KO |
| WT_rep1 | SRR28641566 | GSM8200089 | WT_SBP009_Rep1 | WT |
| WT_rep2 | SRR28641565 | GSM8200090 | WT_SBP009_Rep2 | WT |
| WT_rep3 | SRR28641564 | GSM8200091 | WT_SBP009_Rep3 | WT |

A note on the study design: the three knockouts are three independent CRISPR clones, and the three wild-type samples are replicates of the parental line, all from the same donor (SBP009). Using three separate clones is good, because any change seen in all three is less likely to be a side effect of one clone. But because everything comes from one person, my results can't tell me how much these effects would vary between different people. Also, AKAP11_KO_clone20 is sequenced about twice as deeply as the others (about 294 million read pairs vs 145 to 180 million), which matters for memory in step 8.

> [!WARNING]
> **What could go wrong:** The biggest risk in this whole project is mislabeling a sample (calling a knockout wild type). The run accessions appear in a different order in the ENA table than in GEO, so I matched each SRR to its GEO title through the SRX experiment ID, instead of assuming the order lines up.

## 4 Download and convert the six SRA runs

> [!NOTE]
> **Key idea: What a FASTQ file is.** Each read takes four lines: a name line starting with @, the DNA sequence, a + line, and a quality line with one character per base that encodes how sure the sequencer was about that base. In paired-end data, read number N in the \_1 file and read number N in the \_2 file come from the two ends of the same RNA fragment, so the two files must stay in exactly the same order.

I used the NCBI SRA Toolkit instead of downloading the FASTQ files directly. prefetch downloads one compressed .sra archive per run, vdb-validate checks the archive's internal checksums, fasterq-dump splits it into read 1 and read 2 FASTQ files, and pigz compresses them using several CPU cores at once.

> [!NOTE]
> **Key idea: Slurm job arrays.** Instead of writing six nearly identical scripts, a job array (#SBATCH --array=1-6%2) runs the same script six times. Each copy gets its own $SLURM_ARRAY_TASK_ID (1 to 6), and the script uses it to read the matching line of sra_runs.txt with sed -n "${SLURM_ARRAY_TASK_ID}p". The %2 means at most two tasks run at the same time. %A in the log name is the array's job ID and %a is the task number.

I saved the array script as scripts/download_sra.sbatch.

```bash
#!/bin/bash
#SBATCH --job-name=akap11_sra
#SBATCH --account=gtw111
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --constraint=lustre
#SBATCH --array=1-6%2
#SBATCH --output=logs/sra_%A_%a.out
#SBATCH --error=logs/sra_%A_%a.err

set -euo pipefail
module load sratoolkit/2.10.9/rn4humf

cd "/expanse/lustre/projects/gtw111/$USER/akap11_rnaseq"

root="$PWD"
# Pick this task's run from the list, and stop if it doesn't look like an SRR ID
run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metadata/sra_runs.txt)
[[ "$run" =~ ^SRR[0-9]+$ ]] || exit 1

# Work in a hidden staging folder and only rename it to the real name at the end,
# so a half-finished run never looks finished
out="$root/fastq_sra/$run"
stage="$root/fastq_sra/.${run}_${SLURM_JOB_ID}"
tmp="$root/sra_tmp/${run}_${SLURM_JOB_ID}"

mkdir -p "$root/sra_cache" "$root/fastq_sra" "$tmp"

# Skip runs that already finished (makes the script safe to rerun)
if [ -f "$out/COMPLETE" ]; then
    echo "$run already completed"
    exit 0
fi

# Don't overwrite a folder that exists but never finished
if [ -e "$out" ]; then
    echo "Output folder already exists without COMPLETE marker: $out"
    exit 1
fi

mkdir "$stage"

echo "$(date -Is) Starting $run"

# Search for the archive, because prefetch saves it one folder deeper than you'd expect
# (sra_cache/SRR.../SRR.../SRR....sra)
archive=$(find "$root/sra_cache" -maxdepth 3 \
    -type f -name "$run.sra" -print -quit)

if [ -z "$archive" ]; then
    echo "$(date -Is) Downloading $run"
    prefetch "$run" --max-size 100G \
        -O "$root/sra_cache/$run"

    archive=$(find "$root/sra_cache" -maxdepth 3 \
        -type f -name "$run.sra" -print -quit)
fi

if [ -z "$archive" ]; then
    echo "SRA archive not found for $run"
    exit 1
fi

# Check the archive is complete and not corrupted
echo "$(date -Is) Validating $archive"
vdb-validate "$archive"

# Split into read 1 and read 2 FASTQ files, using the job's 6 CPUs
echo "$(date -Is) Converting $run"
fasterq-dump "$archive" \
    --split-files \
    -e "$SLURM_CPUS_PER_TASK" \
    -t "$tmp" \
    -O "$stage" \
    -p

# Make sure both FASTQ files exist and aren't empty
test -s "$stage/${run}_1.fastq"
test -s "$stage/${run}_2.fastq"

# Compress with pigz, then test the gzip files
echo "$(date -Is) Compressing $run"
pigz -p "$SLURM_CPUS_PER_TASK" \
    "$stage/${run}_1.fastq" \
    "$stage/${run}_2.fastq"

gzip -t "$stage/${run}_1.fastq.gz"
gzip -t "$stage/${run}_2.fastq.gz"

# Mark the run as done and move it into place
touch "$stage/COMPLETE"
mv "$stage" "$out"

echo "$(date -Is) Completed $run"
ls -lh "$out"
# Check the tools exist, then submit the array
module load sratoolkit/2.10.9/rn4humf
command -v prefetch fasterq-dump vdb-validate pigz
bash -n scripts/download_sra.sbatch

# --parsable makes sbatch print only the job ID, so I can save it in a variable
array_id=$(sbatch --parsable scripts/download_sra.sbatch)
echo "$array_id"

# Watch it, then check every task's final state
squeue -j "$array_id"
sacct -j "$array_id" -X --format=JobID,State,Elapsed,ExitCode

# When it's done: six COMPLETE markers and 12 FASTQ files
find fastq_sra -maxdepth 2 -name COMPLETE -print
find fastq_sra -maxdepth 2 -name '*.fastq.gz' | wc -l
ls -lh fastq_sra/SRR28641567/
```

> [!TIP]
> **Expected outcome:** Every array task finishes as COMPLETED, in roughly 20 to 55 minutes each. There are six COMPLETE markers and 12 FASTQ files. For example, fastq_sra/SRR28641567/ has a 4.4 GB R1 and a 4.9 GB R2 file, and fasterq-dump reports 155,374,783 spots (310,749,566 reads) for that run.

> [!WARNING]
> **What could go wrong:** Downloading the FASTQ files straight from ENA with curl is possible, but very large transfers (4 to 15 GB per file) can drop partway through, or finish with a file that fails its MD5 check. The SRA route is more reliable because vdb-validate checks the archive before anything is converted. Two other traps: fasterq-dump writes uncompressed FASTQ that is about three to four times bigger than gzipped, so running all six at once could fill the quota (that's what %2 prevents), and pointing vdb-validate at sra_cache/SRR.../SRR....sra fails because prefetch nests the file one folder deeper (that's why the script uses find).

> [!IMPORTANT]
> **What I learned**
> - An MD5 or internal checksum is the only real proof a downloaded file is complete. The right file size doesn't prove anything.
> - The staging-folder-plus-COMPLETE-marker pattern makes a script safe to rerun: finished samples are skipped, and half-finished ones are never confused with real outputs.
> - A job array is the cleanest way to run the same step on many samples.

## 5 Write the sample sheet

> [!NOTE]
> **Key idea: The nf-core sample sheet.** nf-core/rnaseq reads a CSV with four columns: sample (the name I want in all outputs), fastq_1 and fastq_2 (full paths to read 1 and read 2) and strandedness. Strandedness describes whether the library kept track of which DNA strand the RNA came from. Setting it to "auto" lets the pipeline work it out from the data, instead of me guessing. I kept the condition (KO or WT) in a separate sample_metadata.csv, because the pipeline doesn't need it, but my differential expression analysis will.

This script writes both files from the mapping in step 3. It refuses to write anything unless all 12 FASTQ files exist and aren't empty.

```bash
python - <<'PY'
import csv
from pathlib import Path

root = Path.cwd()

# (my label, condition, SRA run) for each sample, from the GEO mapping in step 3
samples = [
    ("AKAP11_KO_clone16", "KO", "SRR28641569"),
    ("AKAP11_KO_clone20", "KO", "SRR28641568"),
    ("AKAP11_KO_clone21", "KO", "SRR28641567"),
    ("WT_rep1", "WT", "SRR28641566"),
    ("WT_rep2", "WT", "SRR28641565"),
    ("WT_rep3", "WT", "SRR28641564"),
]

# Safety check: six different runs, no duplicates
if len({run for _, _, run in samples}) != 6:
    raise SystemExit("Expected six different SRA runs.")

rows = []
for sample, group, run in samples:
    read1 = root / "fastq_sra" / run / f"{run}_1.fastq.gz"
    read2 = root / "fastq_sra" / run / f"{run}_2.fastq.gz"
    # Stop if any file is missing or empty
    for path in (read1, read2):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty file: {path}")
    rows.append([sample, str(read1), str(read2), "auto"])

# nf-core sample sheet
with open(root / "metadata/samplesheet.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
    writer.writerows(rows)

# Metadata for my own downstream analysis (sample -> condition)
with open(root / "metadata/sample_metadata.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["sample", "condition", "run"])
    writer.writerows(samples)

print("Wrote samplesheet.csv and sample_metadata.csv; all 12 FASTQ files exist and are nonempty.")
PY

# Check both files
cat metadata/samplesheet.csv
cat metadata/sample_metadata.csv
```

> [!TIP]
> **Expected outcome:** Both files have a header plus six rows: three KO clones and three WT replicates. Every row of the sample sheet points to fastq_sra/&lt;run&gt;/&lt;run&gt;\_1.fastq.gz and \_2.fastq.gz.

> [!WARNING]
> **What could go wrong:** Typing 12 long paths by hand is an easy way to swap two samples without noticing, and the pipeline would run happily on the wrong labels. Writing the sheet from one small table in a script, with checks, avoids that. The pipeline also stops at the start if a path in the sheet doesn't exist, which is why the script checks the files first.

## 6 Configure Nextflow for Slurm and Singularity

> [!NOTE]
> **Key idea: Nextflow, nf-core and containers.** Nextflow is a workflow manager: it reads a pipeline, works out which tasks can run and in what order, and submits them. nf-core is a community that writes and tests standard pipelines, and nf-core/rnaseq is the one for bulk RNA-seq. Singularity runs each tool inside a container, a sealed package holding that exact tool version and everything it needs, so I don't have to install STAR, Salmon and the rest myself. The config file is where I tell Nextflow how to do all this on Expanse. When the pipeline runs, there are two kinds of jobs: one long "launcher" job running Nextflow, and many short "child" jobs (named nf-NFCORE_RNASEQ\_...) that the launcher submits for each task.

This config file makes Nextflow submit every pipeline task as its own Slurm job on the shared partition, run each tool inside a Singularity container, and give the tools a temporary folder they can actually use.

```bash
cat > configs/expanse.config <<'EOF'
// Submit every pipeline task to Slurm on the shared partition.
// Expanse rejects jobs that don't set ntasks, so --nodes=1 --ntasks=1 are required.
process {
    executor = 'slurm'
    queue = 'shared'
    clusterOptions = '--account=gtw111 --constraint=lustre --nodes=1 --ntasks=1'
}

// Keep at most two pipeline jobs in the queue, and submit at most one per second
executor {
    queueSize = 2
    submitRateLimit = '1 sec'
}

// Give tools a temporary folder in my project instead of the node's local /scratch
env {
    TMPDIR = '/expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq/tmp'
    SINGULARITYENV_TMPDIR = '/expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq/tmp'
}

// Run tools in Singularity containers, and mount the tmp folder inside them
singularity {
    enabled = true
    autoMounts = true
    runOptions = '--bind /expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq/tmp:/expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq/tmp'
}
EOF
```

queueSize = 2 means Nextflow keeps at most two pipeline jobs in Slurm at any time. That limits how many cores, how much memory and how much temporary disk space the pipeline uses at once, and it avoids flooding the shared queue. The cost is speed: tasks for different samples wait for each other. With enough allocation and quota, raising it (for example to 20) lets several samples run in parallel.

Before trusting the config, I checked by hand that a container (here the bedtools image the pipeline uses) can see and write to the mounted temporary folder.

```bash
image="$PWD/singularity_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-63-6397750e9730a3fbcc5b4c43f14bd141c64c723fd7dad80e47921a68a7c3cd21-data.img"

# Inside the container: check TMPDIR exists, make a test file there, then delete it
SINGULARITYENV_TMPDIR="$PWD/tmp" singularity exec \
    --bind "$PWD/tmp:$PWD/tmp" \
    "$image" sh -c '
        set -eu
        test -d "$TMPDIR"
        probe=$(mktemp "$TMPDIR/check.XXXXXX")
        rm "$probe"
        echo "Temporary folder is accessible and writable inside the container."
    '
```

> [!TIP]
> **Expected outcome:** The container prints "Temporary folder is accessible and writable inside the container."

> [!WARNING]
> **What could go wrong:** Each line of this config prevents a specific failure. Without --ntasks=1, every child job is rejected with "bank_limit plugin: Please set the ntasks or ntasks-per-node parameter". Without TMPDIR, steps that pipe into sort (like BEDTools genomecov) fail with "sort: cannot create temporary file in '/scratch/aaly1/job\_...'", because the node's local scratch folder doesn't exist inside the container. And setting TMPDIR alone isn't enough: a container only sees folders that are mounted into it, so without --bind the same error comes back with "No such file or directory" for my project tmp folder.

> [!IMPORTANT]
> **What I learned**
> - The same pipeline that runs fine on a laptop can fail on a cluster because of site-specific rules, like Expanse requiring ntasks.
> - Containers are isolated. If a tool needs a folder, the folder has to exist and be mounted inside the container.
> - Testing a container by hand with singularity exec takes seconds. Finding the same problem inside the pipeline can cost hours.

## 7 Test the pipeline with the nf-core test profile

Before running my real data, I ran nf-core's built-in test dataset. This checks that Slurm submission, the container downloads and my config all work, without spending hours on six full samples. I saved this launcher as scripts/test_rnaseq.sbatch. The launcher itself only needs 1 CPU and 6 GB, because it just runs Nextflow and Nextflow submits the real work as separate jobs.

```bash
#!/bin/bash
#SBATCH --job-name=rnaseq_test
#SBATCH --account=gtw111
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --constraint=lustre
#SBATCH --output=logs/rnaseq_test_%j.out
#SBATCH --error=logs/rnaseq_test_%j.err

set -euo pipefail

cd "/expanse/lustre/projects/gtw111/$USER/akap11_rnaseq"

# Pin the Nextflow version, keep containers in my project, limit Java memory
export NXF_VER=25.04.8
export NXF_SINGULARITY_CACHEDIR="$PWD/singularity_cache"
export NXF_OPTS="-Xms512m -Xmx4g"

# Separate work and results folders, so test files never mix with my real data
mkdir -p "$NXF_SINGULARITY_CACHEDIR" work_test results/test

nextflow run nf-core/rnaseq \
    -r 3.19.0 \
    -profile test,singularity \
    -c configs/expanse.config \
    -work-dir "$PWD/work_test" \
    --outdir "$PWD/results/test" \
    --max_cpus 6 \
    --max_memory '24.GB' \
    --max_time '4.h' \
    -resume
# Submit the test (job 55623214) and follow it
sbatch scripts/test_rnaseq.sbatch
squeue -u "$USER"
sacct -j 55623214 --format=JobID,State,Elapsed,ExitCode
tail -n 30 logs/rnaseq_test_55623214.out
```

> [!TIP]
> **Expected outcome:** The log ends with "Pipeline completed successfully": duration 58 min 40 s, 209 tasks succeeded.

> [!IMPORTANT]
> **What I learned**
> - The test profile checks the whole system (scheduler, containers, config) in about an hour, before I spend many hours and SUs on real data.
> - Using a separate work_test/ and results/test/ folder means test outputs never mix with my real results.
> - A passing test doesn't guarantee the real run will pass: real samples are thousands of times bigger, so time and memory limits only show up with real data.

## 8 Run the pipeline on my six samples

This is the main step, where nf-core/rnaseq takes my 12 FASTQ files and turns them into a gene count table and a full set of quality-control reports. I didn't write the analysis tools myself: the pipeline chains together about 20 published tools and runs each one as its own Slurm job in its own container.

> [!NOTE]
> **Key idea: What happens inside nf-core/rnaseq.** The pipeline goes through these stages for every sample. (1) Check and clean the reads: FQ_LINT checks the FASTQ files are well formed, FastQC measures read quality, and Trim Galore cuts off adapter sequence and low-quality ends. (2) Work out the strandedness: because my sample sheet says "auto", the pipeline takes a small random subset of reads (FQ_SUBSAMPLE) and runs Salmon on them to see which strand the reads come from. (3) Build the indexes once: STAR_GENOMEGENERATE builds a genome index and SALMON_INDEX builds a transcript index, which are like a book's index that lets the tools find matches fast. (4) Align: STAR maps every read to the genome. It is "splice-aware", which matters for RNA because a read can span two exons with an intron cut out in between. (5) Quantify: Salmon uses STAR's alignments to the transcripts to estimate how many reads came from each transcript, and tximport sums transcripts into genes. (6) Run quality control on the aligned reads: Picard MarkDuplicates, Qualimap, dupRadar, RSeQC, featureCounts biotype QC, StringTie and BEDTools/bigWig coverage tracks. (7) MultiQC collects every QC result into one HTML report.

### 8a Give read-duplication QC enough memory

RSeQC's read_duplication.py keeps a record of every read position in memory, so its memory use grows with sequencing depth. The pipeline's default for this process (36 GB) isn't enough for these deep samples, especially AKAP11_KO_clone20, which has about 294 million read pairs and can need more than 64 GB. This small config file gives that one process 96 GB and leaves every other process at its default.

```bash
# withName picks out pipeline processes by name. The pattern '.*:RSEQC_READDUPLICATION'
# matches the read-duplication step wherever it sits in the pipeline,
# so only that process gets the extra memory.
cat > configs/readdup_72g.config <<'EOF'
process {
    withName: '.*:RSEQC_READDUPLICATION' {
        memory = 96.GB
    }
}
EOF
```

Why the file is called readdup_72g.config: I first set this value to 72 GB. Clone 20 ran out of memory twice at 72 GB (jobs 55736969 and 55742427), so I raised the value in the same file to 96 GB and resubmitted with -resume. The other samples had already finished and were reused, and clone 20 then finished at 96 GB (job 55746476). If you repeat this workflow, you can name the file readdup_96g.config, as long as the -c option in the launcher uses the same name.

### 8b Make biotype QC work with the GENCODE annotation

The pipeline's biotype QC (featureCounts) looks for an attribute called gene_biotype in the GTF, which is the name Ensembl uses. GENCODE calls the same thing gene_type. I checked my GTF to confirm:

```bash
# Print the attributes of the first three exons that have a gene_type field
awk -F '\t' '$3 == "exon" && $9 ~ /(^|;[[:space:]]*)gene_type "/ {
    print $9
    if (++n == 3) exit
}' reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf
```

> [!TIP]
> **Expected outcome:** Each line includes gene_type "lncRNA" (and no gene_biotype), so I tell the pipeline to use gene_type with --featurecounts_group_type gene_type.

### 8c The launcher

I saved the launcher for my real data as scripts/run_rnaseq.sbatch. It uses the GENCODE v47 reference from step 2, STAR plus Salmon (star_salmon), all seven RSeQC modules, and the two config files.

```bash
#!/bin/bash
# The #SBATCH lines are for the "launcher" job, which only runs Nextflow itself.
# It needs 1 CPU and 6 GB because Nextflow submits the heavy work as separate jobs,
# but it has to stay alive until the very last task finishes, so I give it 48 hours.
#SBATCH --job-name=akap11_rnaseq
#SBATCH --account=gtw111
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=48:00:00
#SBATCH --constraint=lustre
#SBATCH --output=logs/rnaseq_%j.out
#SBATCH --error=logs/rnaseq_%j.err

set -euo pipefail

cd "/expanse/lustre/projects/gtw111/$USER/akap11_rnaseq"

# Pin the Nextflow version, keep containers in my project, limit Java memory
export NXF_VER=25.04.8
export NXF_SINGULARITY_CACHEDIR="$PWD/singularity_cache"
export NXF_OPTS="-Xms512m -Xmx4g"

mkdir -p "$NXF_SINGULARITY_CACHEDIR" work results

# -r 3.19.0 pins the pipeline version so my results are reproducible.
# Two -c files: Expanse settings, plus extra memory for read duplication.
# --save_reference keeps the STAR and Salmon indexes so I can reuse them later.
# --rseqc_modules lists all seven RSeQC checks, including read_duplication.
# --featurecounts_group_type gene_type makes biotype QC work with GENCODE.
# -resume is a safety net: if the launcher is ever interrupted, rerunning it
#   reuses every task that already finished instead of starting over.
nextflow run nf-core/rnaseq \
    -r 3.19.0 \
    -profile singularity \
    -c configs/expanse.config \
    -c configs/readdup_72g.config \
    -work-dir "$PWD/work" \
    --input "$PWD/metadata/samplesheet.csv" \
    --outdir "$PWD/results/rnaseq_with_readdup" \
    --fasta "$PWD/reference/gencode_v47/GRCh38.primary_assembly.genome.fa" \
    --gtf "$PWD/reference/gencode_v47/gencode.v47.primary_assembly.annotation.gtf" \
    --aligner star_salmon \
    --save_reference \
    --max_cpus 16 \
    --max_memory '120.GB' \
    --max_time '12.h' \
    --rseqc_modules bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication \
    --featurecounts_group_type gene_type \
    -resume
# bash -n checks the script for syntax mistakes without running it;
# && means sbatch only runs if that check passes.
bash -n scripts/run_rnaseq.sbatch && sbatch scripts/run_rnaseq.sbatch
```

> [!TIP]
> **Expected outcome:** Slurm returns a job ID (mine was 55746106), and within a few minutes child jobs named nf-NFCORE_RNASEQ\_... start appearing in squeue.

> [!NOTE]
> **Key idea: How -resume works.** Every task Nextflow runs gets its own folder in work/ (like work/70/376871...), named by a hash of the task's inputs, script and settings. With -resume, Nextflow checks each task before it runs: if a finished task with the same hash already exists, it reuses that output ("cached") instead of running it again. That's why I never delete work/ until I'm completely done with the project.

> [!WARNING]
> **What could go wrong:** Three things can stop a run like this. (1) Memory: if a task uses more RAM than it asked for, Slurm kills it with exit code 137 and sacct shows OUT_OF_MEMORY. With the default 36 GB, read_duplication would be killed on these deep samples, and even 72 GB isn't enough for clone20. That's what readdup_72g.config prevents. (2) Time: Qualimap alone can take more than 3 hours per sample on these BAM files, so a launcher with a 24-hour limit would hit TIMEOUT before the pipeline finishes. That's why the launcher asks for 48 hours. (3) Biotype QC: without --featurecounts_group_type gene_type, the log shows "Biotype attribute 'gene_biotype' not found in the last column of the GTF file!" and that QC is silently skipped. In every case the fix is the same: change the setting and resubmit with -resume, so the finished work is kept.

> [!WARNING]
> **What could go wrong:** --max_cpus, --max_memory and --max_time had no effect in this run (and in the test run in step 7). These options come from older nf-core pipelines. The pipeline started without an error, but the workflow summary at the end of my MultiQC report, which lists every option that differs from the pipeline's defaults, doesn't include them, so nf-core/rnaseq 3.19.0 ignored them. Newer nf-core pipelines set these limits with process.resourceLimits in a config file instead, as shown below. The extra memory for read duplication still worked, because it is set for that one process in its own config file.

```bash
// To cap resources in a future run, add this line inside the
// process { } block of configs/expanse.config
resourceLimits = [cpus: 16, memory: 120.GB, time: 12.h]
```

> [!IMPORTANT]
> **What I learned**
> - Real data needs much more time and memory than the test data, and the biggest sample sets the requirement.
> - Exit code 137 means the process was killed (128 + signal 9). On a cluster, that almost always means it ran out of memory.
> - Give extra resources only to the process that needs them (withName), instead of making every job bigger and waiting longer in the queue.
> - Check that annotation attribute names match what the pipeline expects. GENCODE and Ensembl use different names for the same thing.

## 9 Monitor the run and confirm it finished

I used these commands to keep an eye on the pipeline while it ran. A wide squeue format shows the full Nextflow task names, so I can tell which step and which sample is running, and the %m field shows how much memory each job requested.

```bash
# My jobs with full names, state, time used, memory requested and node
squeue -u "$USER" -o "%.18i %.110j %.2t %.10M %.12m %.30R"

# Time used and time left for the launcher
squeue -j 55746106 -o "%.18i %.2t %.10M %.12L %.30R"

# List the read-duplication tasks: how much memory each asked for, and whether it finished
sacct -X -S 2026-11-15 -u "$USER" -P \
  --format=JobID,JobName%150,State,ReqMem |
  grep -E 'JobID|RSEQC_READDUPLICATION'

# Memory actually used (MaxRSS) by one finished task, e.g. clone20's read duplication
sacct -j 55746476 --format=JobID,State,Elapsed,ReqMem,MaxRSS,ExitCode

# Final state of the launcher, and the end of its log
sacct -X -j 55746106 --format=JobID,State%20,Elapsed,Timelimit,ExitCode
tail -n 50 logs/rnaseq_55746106.out
```

> [!TIP]
> **Expected outcome:** sacct shows the launcher as COMPLETED with exit code 0:0, every read-duplication task as COMPLETED (clone20's with 96G requested, the others with 72G), and the log ends with "Pipeline completed successfully" (15 Nov 2026).

If a task does fail, Nextflow prints its work folder. sacct can also find it from the Slurm job ID, and the files inside show exactly what happened:

```bash
# Find a task's work folder from its Slurm job ID, then read its log and exit code
task_dir=$(sacct -X -n -P -j 55746476 --format=WorkDir)
tail -n 30 "$task_dir/.command.log"
cat "$task_dir/.exitcode"

# Live memory use of a task that is still running
sstat -j 55746476.batch --format=JobID,AveCPU,AveRSS,MaxRSS,AveDiskRead,AveDiskWrite
```

Learning to read Slurm's short codes makes monitoring much easier:

| Code or state | Where I see it | What it means |
|---|---|---|
| PD (Priority) | squeue | Pending: waiting in the queue for a free node |
| R | squeue | Running right now |
| COMPLETED, exit 0:0 | sacct | Finished without errors |
| FAILED | sacct | The script ended with an error (non-zero exit code) |
| TIMEOUT | sacct | Hit its --time limit and was stopped |
| OUT_OF_MEMORY | sacct | Used more RAM than it asked for and was killed |
| Exit 137 | Nextflow log | Killed with signal 9 (128 + 9), usually memory |
| Exit 143 | Nextflow log | Stopped with signal 15 (128 + 15), a polite stop request |

> [!WARNING]
> **What could go wrong:** Two messages in the log look scary but aren't fatal. "WARN: [SLURM] queue (shared) status cannot be fetched ... exit status : 143" means one squeue check by Nextflow timed out, and Nextflow simply tries again. "sbatch: error: Batch job submission failed: Unexpected message received" means Slurm was briefly overloaded. The job often gets submitted anyway, so I check squeue before resubmitting anything, otherwise I'd end up with duplicate jobs.

> [!IMPORTANT]
> **What I learned**
> - squeue shows what's happening now. sacct shows what happened, including jobs that already ended. sstat shows live memory use for a running job.
> - MaxRSS in sacct is the peak memory a job really used, which is how I can tell whether a memory request was too low or wastefully high.
> - Every Nextflow task folder has .command.sh (what ran), .command.log (what it printed) and .exitcode (how it ended). Those three files explain almost any failure.

## 10 Find the results

```bash
# Find the gene-level count table I'll use for differential expression,
# and the MultiQC report with all the QC
find results/rnaseq_with_readdup -name 'salmon.merged.gene_counts_length_scaled.tsv' -print
find results/rnaseq_with_readdup -name 'multiqc_report.html' -print
```

> [!TIP]
> **Expected outcome:** The counts table is at results/rnaseq_with_readdup/star_salmon/salmon.merged.gene_counts_length_scaled.tsv, and the MultiQC report is at results/rnaseq_with_readdup/multiqc/star_salmon/multiqc_report.html.

> [!NOTE]
> **Key idea: Which counts file to use.** The star_salmon folder has several gene-level tables. gene_tpm.tsv holds TPM values (transcripts per million). They're good for comparing expression within a sample or for plots, but not for statistical testing. gene_counts.tsv holds Salmon's estimated read counts per gene. gene_counts_length_scaled.tsv holds those counts corrected for differences in average transcript length between samples, which matters when a gene's dominant isoform changes between conditions. I chose the length-scaled table for differential expression with DESeq2, because DESeq2 expects count-like values and does its own normalization.

## 11 Copy the gene counts and sample sheet to my Mac

I opened a local Terminal window on my Mac and used scp to copy the files from Expanse into my Downloads folder. I copied the gene counts for downstream analysis, and sample_metadata.csv from step 5, which says which sample is a knockout and which is wild type. I saved the counts table under a name that starts with AKAP11\_, so it can't get mixed up with files from other projects.

> [!NOTE]
> **Key idea: How scp works.** scp (secure copy) uses the same login as ssh. Its pattern is scp FROM TO, and the remote side is written as user@host:/path. Because I run it on my Mac, it logs in to Expanse (password plus TOTP code again), reads the file there, and writes it to ~/Downloads on my Mac. I put quotes around "$remote/..." so the path is treated as one piece.

```bash
# Run these on my Mac, NOT on Expanse (the prompt should show my Mac's name, not login.expanse)
# Save the Expanse project path in a variable so the scp lines stay short
remote='aaly1@login.expanse.sdsc.edu:/expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq'

# Copy the gene counts table into Downloads, and give it a clearer name on my Mac
# (if the destination is a file name instead of a folder, scp saves it under that name)
scp "$remote/results/rnaseq_with_readdup/star_salmon/salmon.merged.gene_counts_length_scaled.tsv" \
  ~/Downloads/AKAP11_final_gene_counts_length_scaled.tsv

# Copy the sample-to-group table from step 5 (samplesheet.csv only has file paths, not groups)
scp "$remote/metadata/sample_metadata.csv" ~/Downloads/
```

> [!TIP]
> **Expected outcome:** The counts transfer reaches 100% (6,741 KB), and AKAP11_final_gene_counts_length_scaled.tsv and sample_metadata.csv are in my Mac's Downloads folder.

> [!WARNING]
> **What could go wrong:** If scp is run while still logged in to Expanse, it tries to save to a Downloads folder on Expanse instead of on my Mac ("open local ... Is a directory" or "No such file"). Always check the prompt first: scp runs from the Mac. Also, Expanse locks you out briefly after several wrong password or TOTP attempts ("Permission denied"), so if a login fails twice, wait a few minutes before trying again.

The Python analysis reads a tab-separated sample sheet called samples.tsv. I made it on my Mac from sample_metadata.csv, adding the GEO sample ID from step 3 so every row can be traced back to GEO:

```bash
# Run on my Mac, in the folder with sample_metadata.csv
cd ~/Downloads

python3 - <<'PY'
import csv

# GEO sample ID for each sample, from the mapping table in step 3
geo_ids = {
    "AKAP11_KO_clone16": "GSM8200086",
    "AKAP11_KO_clone20": "GSM8200087",
    "AKAP11_KO_clone21": "GSM8200088",
    "WT_rep1": "GSM8200089",
    "WT_rep2": "GSM8200090",
    "WT_rep3": "GSM8200091",
}

# Read sample_metadata.csv (columns: sample, condition, run)
with open("sample_metadata.csv") as handle:
    rows = list(csv.DictReader(handle))

# Write samples.tsv with tabs between the columns
with open("samples.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample", "condition", "geo_sample", "sra_run"])
    for row in rows:
        writer.writerow([
            row["sample"],
            row["condition"],
            geo_ids[row["sample"]],
            row["run"],
        ])

print("Wrote samples.tsv with", len(rows), "samples.")
PY

cat samples.tsv
```

> [!TIP]
> **Expected outcome:** The script prints "Wrote samples.tsv with 6 samples." and samples.tsv has a header (sample, condition, geo_sample, sra_run) plus one row per sample, with the same pairings as the table in step 3.

## 12 Copy and open the MultiQC report

Next I copied the MultiQC report to my Mac so I could review the sequencing and alignment quality before analyzing the gene counts. Like the counts table, I saved it with an AKAP11\_ name.

```bash
# On my Mac: same remote path variable as in step 11
remote='aaly1@login.expanse.sdsc.edu:/expanse/lustre/projects/gtw111/aaly1/akap11_rnaseq'

# Copy the MultiQC report to Downloads under a clearer name,
# then open it in my default browser
scp "$remote/results/rnaseq_with_readdup/multiqc/star_salmon/multiqc_report.html" \
  ~/Downloads/AKAP11_MultiQC_Final.html
open ~/Downloads/AKAP11_MultiQC_Final.html
```

> [!TIP]
> **Expected outcome:** The transfer reaches 100% (6,188 KB), and AKAP11_MultiQC_Final.html opens in my default browser.

When I open the report, these are the things I check before trusting the counts:

| Section of the report | What I look for | Why it matters |
|---|---|---|
| General statistics | Read numbers per sample (clone20 should stand out as about twice as deep) | A sample with far fewer reads gives noisier results |
| STAR alignment | A high percentage of uniquely mapped reads | A low rate can mean contamination or the wrong reference |
| FastQC / Trim Galore | Good base quality, adapters removed | Bad reads can misalign |
| RSeQC infer_experiment | All samples have the same strandedness | A mismatch would mean mixed library types |
| RSeQC read duplication / Picard / dupRadar | Duplication levels look similar between groups | Very different levels could mean a technical problem |
| featureCounts biotypes | Most reads assigned to protein-coding genes, low rRNA | High rRNA means the RNA prep didn't remove it well |
| DESeq2 QC (PCA, sample clustering) | KO samples group together, WT samples group together | First sign of whether genotype is the main source of variation |

This is what my report showed:

| Measure | KO clone16 | KO clone20 | KO clone21 | WT rep1 | WT rep2 | WT rep3 |
|---|---|---|---|---|---|---|
| Read pairs sequenced (millions) | 165.1 | 293.8 | 155.4 | 179.6 | 153.7 | 145.0 |
| Uniquely mapped (STAR) | 96.1% | 95.7% | 95.5% | 94.1% | 95.5% | 95.7% |
| Strandedness (Salmon and RSeQC) | reverse | reverse | reverse | reverse | reverse | reverse |
| Antisense reads (RSeQC) | 71.5% | 76.7% | 72.9% | 76.0% | 80.5% | 77.8% |
| Exonic reads (Qualimap) | 66.3% | 40.0% | 60.9% | 59.0% | 65.9% | 66.9% |
| Intronic reads (Qualimap) | 30.9% | 56.4% | 36.2% | 34.4% | 31.4% | 30.2% |
| Duplicate reads, read 1 (FastQC) | 38.0% | 67.1% | 37.5% | 36.1% | 38.0% | 36.3% |
| Duplication (Picard) | 36.9% | 42.9% | 35.7% | 38.4% | 37.7% | 35.0% |
| GC content, read 1 | 45% | 43% | 46% | 46% | 46% | 46% |
| Protein-coding reads (featureCounts) | 95.3% | 93.9% | 94.7% | 95.5% | 95.2% | 94.7% |
| rRNA reads | &lt; 0.001% | &lt; 0.001% | &lt; 0.001% | &lt; 0.001% | &lt; 0.001% | &lt; 0.001% |

**Overall, the data is good.** More than 94% of reads map uniquely in every sample, rRNA is almost absent, and both Salmon and RSeQC call every library reverse-stranded, which matches the paper's HTSeq setting (stranded = reverse). Every FASTQ file fails FastQC's "per base sequence content" check, but that is normal for RNA-seq: library kits prime cDNA with random hexamers, which aren't perfectly random, so the first bases of every read have uneven base composition.

**Clone 20 differs from the other five samples in more than depth.** 56% of its aligned reads fall in introns (30 to 36% in the others), its FastQC duplication is 67% (36 to 38% in the others), its GC content is lower, and it is the only sample that also fails FastQC's GC-content and duplication checks. Some extra duplication is expected from its deeper sequencing, but a high intronic fraction usually means more unspliced pre-mRNA in the library, for example from a larger share of nuclear RNA or a difference in RNA extraction.

**This shows up in the PCA.** In the DESeq2 PCA in the report, clone 20 sits alone along PC1, while PC2 separates the WT replicates from clones 16 and 21. My own PCA in Python shows the same pattern. So the largest source of variation in this dataset is a technical difference in one sample, not genotype. I kept clone 20, because it passes the basic checks (mapping, strandedness, rRNA) and the paper also used all three clones. But it is the main caveat when I interpret the results: a gene that changes only in clone 20 is weaker evidence than one that changes in all three clones.

## 13 Troubleshooting guide

This table summarizes the problems that can come up in this workflow, how to recognize each one, and how the workflow prevents or fixes it.

| Symptom | Likely cause | Fix used in this workflow |
|---|---|---|
| SSH "Broken pipe" / "Connection reset" | Idle or network timeouts | Run all long work with sbatch, then log back in and check squeue |
| curl exit 18 or an MD5 mismatch on FASTQ downloads | Unstable very large transfers | Download with SRA prefetch and check with vdb-validate (step 4) |
| vdb-validate: archive "could not be found" | prefetch nests the .sra one folder deeper | Locate the archive with find (step 4) |
| "Please set the ntasks or ntasks-per-node parameter" | Expanse rule for every job | --nodes=1 --ntasks=1 in clusterOptions (step 6) |
| "sort: cannot create temporary file" | No usable temp folder inside the container | TMPDIR plus --bind (step 6) |
| Launcher shows TIMEOUT | Real data runs much longer than the test | 48-hour launcher; resubmit with -resume (step 8) |
| Exit 137 / OUT_OF_MEMORY in read_duplication | Deep samples need more than 36 GB (clone20 more than 72 GB) | 96 GB for that process only, with withName (step 8) |
| "gene_biotype not found" warning | GENCODE uses gene_type | --featurecounts_group_type gene_type (step 8) |
| "Unexpected message received" from sbatch | Slurm briefly overloaded | Check squeue before resubmitting (step 9) |
| scp saves to the wrong place | Command run on Expanse instead of the Mac | Run scp from the Mac's own Terminal (step 11) |

## 14 From counts to results: the Python analysis

The Expanse part of the project ends with the gene count table and the sample information on my Mac. Everything after that runs on my laptop, with the Python code in this repository. The [main README](../README.md) explains how to install and run it.

### What goes in

| File | Where it comes from |
|---|---|
| AKAP11_final_gene_counts_length_scaled.tsv | salmon.merged.gene_counts_length_scaled.tsv from step 11, renamed |
| samples.tsv | Made from sample_metadata.csv in step 11 |

### How to run it

```bash
# Run in the project folder on my Mac
python3 run_analysis.py AKAP11_final_gene_counts_length_scaled.tsv samples.tsv WT KO akap11_go_bp_results AKAP11
```

WT is the reference group and KO is the comparison group, so a positive log2 fold change means higher in the knockout. The last word, AKAP11, is an optional gene to plot on its own.

### What the analysis does

| Step | File | What it does | Main outputs |
|---|---|---|---|
| 1 | read_counts.py | Checks the inputs, rounds the counts, and keeps genes with ≥ 10 counts in ≥ 3 samples | 26,288 of 78,932 genes kept |
| 2 | fit_model.py | PyDESeq2 with the design ~condition. Significant = adjusted p ≤ 0.05 and \|log2FC\| ≥ 1, the same cutoffs as the paper | 01_sample_qc.tsv, 02_differential_expression.tsv, 02_significant_genes.tsv |
| 3 | plot_pca.py | PCA of all genes (log2 CPM) and sample correlations | 03_pca.png, 04_sample_correlation.png |
| 4 | plot_volcano.py, plot_ma.py | Volcano and MA plots | 05_volcano.png, 06_ma.png |
| 5 | plot_expression_heatmap.py | Clustered heatmap of the top 50 significant genes | 07_de_heatmap.png |
| 6 | plot_target_gene.py | AKAP11 expression in each sample | 08_target_gene.png |
| 7 | run_gsea.py, plot_pathways.py | GSEA on all genes, ranked by the Wald statistic | 09_gsea_pathways.tsv, 09_gsea_pathways.png |
| 8 | run_ora.py, plot_ora.py | Reactome over-representation analysis of the significant genes (paper Fig. 2E) | 10_reactome_ora.tsv, 10_reactome_ora.png |
| 10 | run_paper_gsea.py | A second GSEA with GO Molecular Function and gene sets of up to 3,000 genes, as in the paper's Fig. 2F | 12_fig2F_gsea_go_mf.tsv |
| 11 | compare_to_paper.py | Automatic check against the paper's Figure 2 | 11_comparison_with_paper.md |
| 12 | plot_paper_figures.py | Redraws each panel of the paper's Figures 1 and 2 and explains every difference | 12_paper_figures_report.md, 12_fig\*.png |

I used PyDESeq2 rather than R, so the whole downstream analysis is one Python program that anyone can run with a single command. DESeq2 is designed for count data like Salmon's length-scaled counts.

### How my results compare with the paper

compare_to_paper.py runs automatically at the end of the analysis and checks my results against tables copied from the paper's Figure 2:

| Check | Paper | Mine | Verdict |
|---|---|---|---|
| Significant genes (Fig. 2B) | 849 (421 down, 428 up) | 848 (475 down, 373 up) | Similar |
| Paper's top 60 genes (Fig. 2D) | 30 down, 30 up | 55 of them are in my data. All 55 change in the same direction, 37 are significant, and the fold changes correlate at r = 0.98 | Agree |
| Reactome pathways (Fig. 2E) | GPCR signalling, extracellular matrix | 7 of their top 10 pathways are significant; GPCR pathways are my top 5 | Agree |
| GSEA (Fig. 2F) | Ribosome and protein-synthesis terms down in KO | 16 ribosome/protein-synthesis terms significantly down, none up | Agree |

I also ran GSEA with GO Molecular Function, the collection the paper used (gsea_gene_sets = "GO_Molecular_Function_2025" in settings.py). 7 of the paper's 20 Fig. 2F terms are in that run's results, and all 7 change in the same direction as in the paper. rRNA binding, tRNA binding and aminoacyl-tRNA ligase activity are significantly lower in the knockout, but most of the paper's other terms don't reach FDR &lt; 0.05.

The 5 genes missing from my data are C5orf17, ZBED6CL, C9orf64, LINC02302 and FRG1DP. The first four were renamed between Ensembl 104 (paper) and GENCODE v47 (mine); for example, C9orf64 is now QNG1. FRG1DP has zero reads in all six of my samples.

The total number of significant genes almost matches, but the up/down split doesn't. Likely reasons are the differences listed in section 0, especially the different counting method (Salmon vs HTSeq) and statistical model (DESeq2 vs limma-voom).

AKAP11 itself is about 35% lower in the knockout clones (log2 fold change -0.61), but it doesn't pass the cutoffs (adjusted p = 0.15). That fits a heterozygous knockout, where one working copy of the gene remains.

### Reproducing the paper's figures

plot_paper_figures.py redraws every panel of the paper's figures that can be made from RNA-seq data, and 12_paper_figures_report.md explains each difference. Figures 1C, 1D and 3 come from Western blots, microscopy and DNA methylation sequencing, so they can't be made from this data.

| Panel | Status | Main reason for any difference |
|---|---|---|
| Fig. 1B, AKAP11 per clone | Reproduced with differences | RNA-seq counts the whole gene, including the edited copy; the paper's qPCR measures only the targeted exon 7. Both show AKAP11 lower in every clone and highest in clone 20 |
| Fig. 1C, Western blot | Not reproducible from RNA-seq | Measures protein |
| Fig. 1D, immunostaining | Not reproducible from RNA-seq | Microscope images |
| Fig. 2A, volcano plot | Reproduced with differences | DESeq2 gives much smaller p-values for the strongest genes than limma-voom (y-axis up to 72 vs about 7) |
| Fig. 2B, number of significant genes | Reproduced with differences | Same total (848 vs 849); genes near the cutoffs split differently |
| Fig. 2C, heatmap | Reproduced | Clone 20 follows the shared pattern least closely, as in the paper |
| Fig. 2D, top 30 genes | Reproduced with differences | The two methods rank the strongest genes differently; 8 genes overlap |
| Fig. 2E, Reactome ORA | Reproduced with differences | Same leading GPCR pathways; 7 of 10 significant; different Reactome version and background |
| Fig. 2F, GSEA GO Molecular Function | Reproduced with differences (needs internet) | GSEApy with 1,000 permutations can't report p-values as small as fgsea, and GO versions differ |
| Fig. 3, DNA methylation | Not reproducible from RNA-seq | Different experiment (WGBS) with its own data and pipeline |

With only three samples per group, any method can find strong effects but will miss small ones. Clone 20 is also technically different from the other samples (step 12), so I treat a gene that changes the same way in all three clones as stronger evidence than one driven by a single clone.

## Glossary

| Term | Meaning |
|---|---|
| Allocation / SU | My group's computing budget on Expanse. One SU is about one CPU core for one hour |
| BAM | Binary file of aligned reads: where each read landed in the genome |
| Container (Singularity) | A sealed package with one tool and everything it needs, so it runs the same everywhere |
| FASTA | Text file of DNA sequences. My reference genome is one |
| FASTQ | Text file of sequencing reads with a quality score for every base |
| GTF | Annotation file listing the positions of genes, transcripts and exons |
| Heterozygous knockout | One of the two copies of a gene is disrupted |
| Job array | One Slurm submission that runs the same script many times, once per task ID |
| Login node / compute node | The shared machine I ssh into vs. the machines that run batch jobs |
| Lustre | Expanse's large parallel filesystem for project data |
| MaxRSS | The peak memory a job actually used, as reported by sacct or sstat |
| MD5 / SHA-256 | Checksums: fingerprints that prove a file hasn't changed |
| Nextflow / nf-core | Workflow manager / community library of tested pipelines |
| Paired-end | Each RNA fragment is read from both ends, giving two FASTQ files |
| Salmon | Fast tool that estimates how many reads came from each transcript |
| Slurm / sbatch | Expanse's job scheduler / the command that submits a batch script |
| STAR | Splice-aware aligner that maps RNA reads to the genome |
| Strandedness | Whether the library records which DNA strand the RNA was copied from |
| TPM | Transcripts per million: expression normalized for gene length and sequencing depth |
| work/ folder | Where Nextflow keeps every task's files, and what makes -resume possible |

