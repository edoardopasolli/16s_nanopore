# Full-length 16S rRNA Nanopore preprocessing and mapping workflow

This repository provides a reproducible command-line workflow for preprocessing full-length 16S rRNA Nanopore reads and mapping them to a reference database.

The workflow performs:

1. Barcode-based demultiplexing;
2. Sequential barcode/primer trimming;
3. Read length and quality filtering;
4. Quality control with FastQC;
5. Mapping with minimap2;
6. Extraction of primary alignments with samtools;
7. Filtering of high-quality primary alignments based on mapping quality, aligned length, and sequence identity.

The workflow is designed to be run directly from the command line.  
All user-editable parameters are collected in a single configuration block, so the analysis can be adapted to different datasets without modifying the individual commands.

---

## Requirements

The workflow was tested on Linux.

Required software:

- `cutadapt`
- `FastQC`
- `minimap2`
- `samtools`
- `awk`
- `conda`

A single conda environment can be created as follows:

```bash
conda create -y -n 16s_nanopore -c conda-forge -c bioconda \
  cutadapt fastqc minimap2 samtools
```

Activate the environment before running the workflow:

```bash
conda activate 16s_nanopore
```

If you are running on an HPC cluster where these programs are already available as modules, load the equivalent modules instead.

---

## Input files

Place all required input files in the working directory, or use absolute paths in the configuration block below.

Required files:

| File | Description |
|---|---|
| `bamboo22.fastq.gz` | Raw Nanopore reads before demultiplexing |
| `forward_barcodes.fasta` | FASTA file containing barcode sequences |
| `silva_nr99_v138.2_toSpecies_trainset_uq.mmi` | Pre-built minimap2 reference index |

The barcode FASTA headers are used as sample identifiers by `cutadapt`.

For example, if the FASTA file contains headers such as:

```text
>S4
ACGT...
>S5
TGCA...
```

the demultiplexing step will generate files named:

```text
bamboo22.S4.fastq
bamboo22.S5.fastq
```

---

## Edit variables

Before running the workflow, edit the variables in this block.

```bash
# =========================
# User configuration
# =========================

# Project/sample prefix
PROJECT="bamboo22"

# Input files
RAW_FASTQ="bamboo22.fastq.gz"
BARCODE_FASTA="forward_barcodes.fasta"
REF_INDEX="silva_nr99_v138.2_toSpecies_trainset_uq.mmi"

# Sample names expected after demultiplexing
SAMPLES=("S4" "S5" "S6" "S10" "S11" "S12")

# Adapter / primer sequences
BARCODE_SEQ="GGTAGTATATACAGAGAG"
PRIMER_FWD="AGRGTTYGATYMTGGCTCAG"
PRIMER_REV="RGYTACCTTGTTACGACTT"

# Cutadapt parameters
ERROR_RATE="0.1"
MIN_LENGTH="1350"
MAX_LENGTH="1650"
QUALITY_CUTOFF="20,20"
CUTADAPT_THREADS="64"

# Mapping parameters
MINIMAP2_THREADS="32"
MINIMAP2_PRESET="map-ont"

# High-quality alignment filtering parameters
MIN_MAPQ="20"
MIN_ALIGNED_LENGTH="500"
MIN_PID="0.85"

# Output directories
FASTQC_DIR="fastqc_output"
MAPPING_DIR="mapping_output"

# =========================
# End of configuration
# =========================
```

---

## Complete workflow

The following block runs the complete workflow using the parameters defined above.

```bash
echo "Starting 16S Nanopore preprocessing and mapping workflow"
echo "Project: ${PROJECT}"
echo "Samples: ${SAMPLES[*]}"

# -------------------------
# Check input files
# -------------------------

for file in "${RAW_FASTQ}" "${BARCODE_FASTA}" "${REF_INDEX}"; do
  if [[ ! -f "${file}" ]]; then
    echo "ERROR: required file not found: ${file}" >&2
    exit 1
  fi
done

# -------------------------
# Create output directories
# -------------------------

mkdir -p "${FASTQC_DIR}"
mkdir -p "${MAPPING_DIR}"

# -------------------------
# Step 1: Demultiplex reads
# -------------------------

echo "Step 1: demultiplexing reads"

cutadapt \
  -g "file:${BARCODE_FASTA}" \
  -e "${ERROR_RATE}" \
  --rc \
  -j "${CUTADAPT_THREADS}" \
  -o "${PROJECT}.{name}.fastq" \
  "${RAW_FASTQ}"

# -------------------------
# Step 2: Sequential trimming
# -------------------------

echo "Step 2: trimming barcode and primer sequences"

for SAMPLE in "${SAMPLES[@]}"; do

  INPUT_FASTQ="${PROJECT}.${SAMPLE}.fastq"

  if [[ ! -f "${INPUT_FASTQ}" ]]; then
    echo "WARNING: expected demultiplexed file not found: ${INPUT_FASTQ}" >&2
    echo "Skipping sample ${SAMPLE}" >&2
    continue
  fi

  echo "Processing sample ${SAMPLE}"

  cutadapt \
    -g "${BARCODE_SEQ}" \
    -e "${ERROR_RATE}" \
    --rc \
    -j "${CUTADAPT_THREADS}" \
    -o "${PROJECT}.${SAMPLE}.b.fastq" \
    "${PROJECT}.${SAMPLE}.fastq"

  cutadapt \
    -g "${PRIMER_FWD}" \
    -e "${ERROR_RATE}" \
    --rc \
    -j "${CUTADAPT_THREADS}" \
    -o "${PROJECT}.${SAMPLE}.c.fastq" \
    "${PROJECT}.${SAMPLE}.b.fastq"

  cutadapt \
    -g "${PRIMER_REV}" \
    -e "${ERROR_RATE}" \
    --rc \
    -j "${CUTADAPT_THREADS}" \
    -o "${PROJECT}.${SAMPLE}.d.fastq" \
    "${PROJECT}.${SAMPLE}.c.fastq"

  cutadapt \
    --minimum-length "${MIN_LENGTH}" \
    -l "${MAX_LENGTH}" \
    --quality-cutoff "${QUALITY_CUTOFF}" \
    -j "${CUTADAPT_THREADS}" \
    -o "${PROJECT}.${SAMPLE}.e.fastq" \
    "${PROJECT}.${SAMPLE}.d.fastq"

done

# -------------------------
# Step 3: Quality control
# -------------------------

echo "Step 3: running FastQC"

fastqc \
  "${PROJECT}".*.e.fastq \
  -o "${FASTQC_DIR}"

# -------------------------
# Step 4: Mapping
# -------------------------

echo "Step 4: mapping reads with minimap2"

for SAMPLE in "${SAMPLES[@]}"; do

  FINAL_FASTQ="${PROJECT}.${SAMPLE}.e.fastq"

  if [[ ! -f "${FINAL_FASTQ}" ]]; then
    echo "WARNING: final trimmed file not found: ${FINAL_FASTQ}" >&2
    echo "Skipping mapping for sample ${SAMPLE}" >&2
    continue
  fi

  SAMPLE_PREFIX="${MAPPING_DIR}/${PROJECT}.${SAMPLE}.e"

  echo "Mapping sample ${SAMPLE}"

  minimap2 \
    -t "${MINIMAP2_THREADS}" \
    -ax "${MINIMAP2_PRESET}" \
    "${REF_INDEX}" \
    "${FINAL_FASTQ}" \
    > "${SAMPLE_PREFIX}.sam"

  samtools view \
    -b \
    -o "${SAMPLE_PREFIX}.bam" \
    "${SAMPLE_PREFIX}.sam"

  samtools sort \
    -o "${SAMPLE_PREFIX}.sorted.bam" \
    "${SAMPLE_PREFIX}.bam"

  samtools index \
    "${SAMPLE_PREFIX}.sorted.bam"

  # Keep only primary alignments.
  # This removes secondary alignments, flag 0x100,
  # and supplementary alignments, flag 0x800.
  samtools view \
    -h \
    -F 0x900 \
    "${SAMPLE_PREFIX}.sorted.bam" \
    > "${SAMPLE_PREFIX}.sorted.primary.sam"

  samtools view \
    -b \
    -o "${SAMPLE_PREFIX}.sorted.primary.bam" \
    "${SAMPLE_PREFIX}.sorted.primary.sam"

  samtools view \
    "${SAMPLE_PREFIX}.sorted.primary.bam" \
    | awk '{print $3}' \
    > "${SAMPLE_PREFIX}.sorted.primary.alignedseqs.txt"

done

# -------------------------
# Step 5: Filter high-quality primary alignments
# -------------------------
  
echo "Step 5: filtering high-quality primary alignments for sample ${SAMPLE}"

for SAMPLE in "${SAMPLES[@]}"; do

  samtools view \
    -h \
    -q "${MIN_MAPQ}" \
    -F 0x904 \
    "${SAMPLE_PREFIX}.sorted.bam" \
    | awk -v min_len="${MIN_ALIGNED_LENGTH}" -v min_pid="${MIN_PID}" '
      BEGIN { OFS="\t" }

      /^@/ {
        print
        next
      }

      {
        cigar = $6

        # Compute aligned length from CIGAR.
        # Sum M, =, X, I, and D operations.
        # Ignore S, H, P, and N.
        aln = 0
        tmp = cigar

        while (match(tmp, /([0-9]+)([MIDNSHP=X])/, a)) {
          n = a[1]
          op = a[2]

          if (op == "M" || op == "=" || op == "X" || op == "I" || op == "D") {
            aln += n
          }

          tmp = substr(tmp, RSTART + RLENGTH)
        }

        # Parse NM tag, edit distance.
        nm = 0
        if (match($0, /NM:i:([0-9]+)/, m)) {
          nm = m[1]
        }

        pid = (aln > 0) ? (1 - nm / aln) : 0

        # Keep only alignments passing length and identity filters.
        if (aln >= min_len && pid >= min_pid) {
          print
        }
      }
    ' \
    | samtools view \
        -b \
        -o "${SAMPLE_PREFIX}.filtered.q${MIN_MAPQ}.pid85.len${MIN_ALIGNED_LENGTH}.bam" \
        -

  samtools index \
    "${SAMPLE_PREFIX}.filtered.q${MIN_MAPQ}.pid85.len${MIN_ALIGNED_LENGTH}.bam"

  samtools view \
    "${SAMPLE_PREFIX}.filtered.q${MIN_MAPQ}.pid85.len${MIN_ALIGNED_LENGTH}.bam" \
    | awk '{print $3}' \
    > "${SAMPLE_PREFIX}.filtered.q${MIN_MAPQ}.pid85.txt"

done

echo "Workflow completed successfully"
```

---

## Running a single sample

If you want to process only one sample, edit the `SAMPLES` variable in the configuration block:

```bash
SAMPLES=("S4")
```

Then run the complete workflow block.

---

## Running a different project

For a new sequencing run, change the following variables in the configuration block:

```bash
PROJECT="my_project"
RAW_FASTQ="my_project.fastq.gz"
BARCODE_FASTA="my_barcodes.fasta"
REF_INDEX="my_reference_index.mmi"
SAMPLES=("Sample1" "Sample2" "Sample3")
```

The workflow will then generate files named according to the new project and sample names.

---

## Output files

For each sample, the workflow generates intermediate and final files.

For sample `S4`, the trimming outputs are:

| File | Description |
|---|---|
| `bamboo22.S4.fastq` | Demultiplexed reads |
| `bamboo22.S4.b.fastq` | After barcode trimming |
| `bamboo22.S4.c.fastq` | After forward primer trimming |
| `bamboo22.S4.d.fastq` | After reverse primer trimming |
| `bamboo22.S4.e.fastq` | Final reads after length and quality filtering |

FastQC results are written to:

```text
fastqc_output/
```

Mapping results are written to:

```text
mapping_output/
```

For sample `S4`, mapping outputs include:

```text
mapping_output/bamboo22.S4.e.sam
mapping_output/bamboo22.S4.e.bam
mapping_output/bamboo22.S4.e.sorted.bam
mapping_output/bamboo22.S4.e.sorted.bam.bai
mapping_output/bamboo22.S4.e.sorted.primary.sam
mapping_output/bamboo22.S4.e.sorted.primary.bam
mapping_output/bamboo22.S4.e.sorted.primary.alignedseqs.txt
mapping_output/bamboo22.S4.e.filtered.q20.pid85.len500.bam
mapping_output/bamboo22.S4.e.filtered.q20.pid85.len500.bam.bai
mapping_output/bamboo22.S4.e.filtered.q20.pid85.txt
```

The file:

```text
mapping_output/bamboo22.S4.e.sorted.primary.alignedseqs.txt
```

contains the reference identifiers associated with all primary alignments.

The file:

```text
mapping_output/bamboo22.S4.e.filtered.q20.pid85.txt
```

contains the reference identifiers associated with high-quality primary alignments passing the following filters:

- mapping quality ≥ 20;
- aligned length ≥ 500 bp;
- percentage identity ≥ 0.85;
- exclusion of unmapped, secondary, and supplementary alignments.

---

## Parameter description

| Parameter | Meaning | Default value |
|---|---|---|
| `PROJECT` | Prefix used for input and output files | `bamboo22` |
| `RAW_FASTQ` | Raw input FASTQ file | `bamboo22.fastq.gz` |
| `BARCODE_FASTA` | Barcode FASTA file | `forward_barcodes.fasta` |
| `REF_INDEX` | Minimap2 reference index | `silva_nr99_v138.2_toSpecies_trainset_uq.mmi` |
| `SAMPLES` | Sample names expected after demultiplexing | `S4 S5 S6 S10 S11 S12` |
| `BARCODE_SEQ` | Barcode sequence trimmed after demultiplexing | `GGTAGTATATACAGAGAG` |
| `PRIMER_FWD` | Forward primer sequence | `AGRGTTYGATYMTGGCTCAG` |
| `PRIMER_REV` | Reverse primer sequence | `RGYTACCTTGTTACGACTT` |
| `ERROR_RATE` | Maximum adapter/primer matching error rate | `0.1` |
| `MIN_LENGTH` | Minimum retained read length | `500` |
| `MAX_LENGTH` | Maximum retained read length after cropping | `1550` |
| `QUALITY_CUTOFF` | 5' and 3' quality trimming cutoff | `20,20` |
| `CUTADAPT_THREADS` | Number of threads for cutadapt | `64` |
| `MINIMAP2_THREADS` | Number of threads for minimap2 | `32` |
| `MINIMAP2_PRESET` | Minimap2 preset for Nanopore reads | `map-ont` |
| `MIN_MAPQ` | Minimum mapping quality for high-quality alignments | `20` |
| `MIN_ALIGNED_LENGTH` | Minimum aligned length for high-quality alignments | `500` |
| `MIN_PID` | Minimum percentage identity for high-quality alignments | `0.85` |
| `FASTQC_DIR` | Directory for FastQC output | `fastqc_output` |
| `MAPPING_DIR` | Directory for mapping output | `mapping_output` |

---

## Citation

If you use this workflow, please cite:

> Ida Romano, Edoardo Pasolli, Jean-Claude Walser, Valeria Ventorino, Sonja Reinhard, Giuseppina Magaraci, Olimpia Pepe, Natacha Bodenhausen.  
> *A hybrid and cost-efficient barcoding strategy for full-length 16S rRNA Nanopore sequencing of environmental samples.*  
> Under review.

---

## License

This repository is licensed under the GNU General Public License v3.0 or later.
