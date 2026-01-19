# Long‑read Amplicon Preprocessing & Mapping Pipeline  
**(Cutadapt → FastQC → Minimap2 → Samtools)**

This repository provides a reproducible workflow to preprocess long‑read amplicon data (demultiplexing and multi‑step adapter/primer trimming), run QC, and map reads to a reference database to extract primary alignments.

> **Highlights**
> - Demultiplex with barcode FASTA  
> - Multi‑step 5’ trimming (barcodes/primers) with reverse‑complement matching  
> - Length & quality trimming  
> - QC with FastQC  
> - Mapping with minimap2 (`map-ont`) and primary-alignment extraction with samtools

---

## Contents
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Quick Start](#quick-start)
- [Step-by-Step Pipeline](#step-by-step-pipeline)
- [Run All Samples in a Loop](#run-all-samples-in-a-loop)
- [Parameter Notes](#parameter-notes)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)
- [License](#license)

---

## Requirements

Tested on Linux (x86_64). Recommended ≥32 CPU threads and sufficient RAM for indexing/mapping.

**Software (tested versions):**
- `cutadapt` ≥ 4.6
- `FastQC` ≥ 0.12.1
- `minimap2` = 2.28
- `samtools` = 1.20
- `conda` (for environment management)

> **Conda install (recommended):**
```bash
# Create and activate environments (or use your existing ones)
conda create -y -n cutadapt cutadapt=4.6
conda create -y -n fastqc fastqc=0.12.1
# minimap2 and samtools can be installed in a separate env or globally:
conda create -y -n mapping minimap2=2.28 samtools=1.20
```

If you already have site-specific binaries (e.g., on an HPC), ensure they’re on your `PATH` or use absolute paths.

---

## Input Files

Place (or symlink) the following in the working directory:

- **Raw reads (demultiplexing input)**  
  `bamboo22.fastq.gz`

- **Barcode file for demultiplexing**  
  `forward_barcodes.fasta`  
  > FASTA headers should contain the logical sample names; cutadapt will substitute `{name}` with this header for per-barcode outputs.

- **Per‑sample intermediate files** (after demultiplexing), expected names:  
  `bamboo22.S4.fastq`, `bamboo22.S5.fastq`, `bamboo22.S6.fastq`, `bamboo22.S10.fastq`, `bamboo22.S11.fastq`, `bamboo22.S12.fastq`  
  > These are the demultiplexed files that subsequent steps will trim.

- **Reference index for mapping**  
  `silva_nr99_v138.2_toSpecies_trainset_uq.mmi`  
  (minimap2 prebuilt index; see notes if you need to build it from a FASTA)

---

## Output Files

Per sample (e.g., `S4`), after each trimming step:

- `bamboo22.S4.b.fastq` – after removing barcode sequence 1  
- `bamboo22.S4.c.fastq` – after removing primer sequence 2  
- `bamboo22.S4.d.fastq` – after removing primer sequence 3  
- `bamboo22.S4.e.fastq` – after quality & length trimming (final pre‑mapping reads)

QC:
- `fastqc_output/` – FastQC HTML and ZIP for all `*.e.fastq`

Mapping & postprocessing (example for `S4`):
- `bamboo22.S4.e.sam`, `bamboo22.S4.e.bam`, `bamboo22.S4.e.sorted.bam(.bai)`  
- `bamboo22.S4.e.sorted.primary.sam`, `bamboo22.S4.e.sorted.primary.bam`  
- `bamboo22.S4.e.sorted.primary.alignedseqs.txt` (RNAME list of primary alignments)

---

## Quick Start

If you want to run a **single sample** end‑to‑end (e.g., `S4`), do:

```bash
# 1) Demultiplex raw reads by barcode file
conda activate cutadapt
cutadapt -g file:forward_barcodes.fasta -e 0.1 --rc -j 64 \
  -o bamboo22.{name}.fastq bamboo22.fastq.gz

# 2) 5' adapter/primer trimming (three rounds) + length/quality trimming
cutadapt -g GGTAGTATATACAGAGAG -e 0.1 --rc -j 64 -o bamboo22.S4.b.fastq bamboo22.S4.fastq
cutadapt -g AGRGTTYGATYMTGGCTCAG -e 0.1 --rc -j 64 -o bamboo22.S4.c.fastq bamboo22.S4.b.fastq
cutadapt -g RGYTACCTTGTTACGACTT   -e 0.1 --rc -j 64 -o bamboo22.S4.d.fastq bamboo22.S4.c.fastq
cutadapt --minimum-length 500 -l 1550 --quality-cutoff 20,20 -j 64 \
  -o bamboo22.S4.e.fastq bamboo22.S4.d.fastq

# 3) QC on final reads
conda activate fastqc
mkdir -p fastqc_output
fastqc bamboo22.S4.e.fastq -o fastqc_output/

# 4) Mapping + primary alignments extraction (adjust paths or ensure tools in PATH)
sample=bamboo22.S4.e

# If minimap2 is not on PATH, use its full path; otherwise just "minimap2"
minimap2 -t 32 -ax map-ont silva_nr99_v138.2_toSpecies_trainset_uq.mmi $sample.fastq > $sample.sam

samtools view -b -o $sample.bam $sample.sam
samtools sort -o $sample.sorted.bam $sample.bam
samtools index $sample.sorted.bam

# Keep only primary alignments (drop secondary 0x100 and supplementary 0x800)
samtools view -h -F 0x900 ${sample}.sorted.bam > ${sample}.sorted.primary.sam
samtools view -b -o ${sample}.sorted.primary.bam ${sample}.sorted.primary.sam
samtools view ${sample}.sorted.primary.bam | awk '{print $3}' > ${sample}.sorted.primary.alignedseqs.txt
```

---

## Step-by-Step Pipeline

### 1) Demultiplex by barcode

```bash
conda activate cutadapt
cutadapt -g file:forward_barcodes.fasta -e 0.1 --rc -j 64 \
  -o bamboo22.{name}.fastq bamboo22.fastq.gz
```

- `-g file:forward_barcodes.fasta` reads multiple adapters from a FASTA; `{name}` is replaced with the FASTA header names to produce one output per barcode.  
- `--rc` searches the reverse complement if present.  
- `-e 0.1` allows up to 10% error rate in matches.  
- **Outputs** like `bamboo22.S4.fastq`, `bamboo22.S5.fastq`, etc.

### 2) Multi‑step 5’ trimming (barcodes/primers)

For each sample, sequentially trim the three adapter/primer motifs (allowing RC and 10% error):

```bash
# Barcode (round 1)
cutadapt -g GGTAGTATATACAGAGAG -e 0.1 --rc -j 64 \
  -o bamboo22.S4.b.fastq bamboo22.S4.fastq

# Primer (round 2)
cutadapt -g AGRGTTYGATYMTGGCTCAG -e 0.1 --rc -j 64 \
  -o bamboo22.S4.c.fastq bamboo22.S4.b.fastq

# Primer (round 3)
cutadapt -g RGYTACCTTGTTACGACTT -e 0.1 --rc -j 64 \
  -o bamboo22.S4.d.fastq bamboo22.S4.c.fastq
```

### 3) Length and quality trimming

```bash
cutadapt --minimum-length 500 -l 1550 --quality-cutoff 20,20 -j 64 \
  -o bamboo22.S4.e.fastq bamboo22.S4.d.fastq
```

- `--minimum-length 500`: retain reads ≥ 500 nt  
- `-l 1550`: crop reads to at most 1550 nt (3' end clipping to a fixed length)  
- `--quality-cutoff 20,20`: trim low-quality bases from both ends (Phred 20)

### 4) Quality control

```bash
conda activate fastqc
mkdir -p fastqc_output
fastqc *.e.fastq -o fastqc_output/
```

### 5) Mapping and primary alignments

```bash
# Ensure minimap2 and samtools are available (PATH or absolute paths)
sample=bamboo22.S4.e

minimap2 -t 32 -ax map-ont silva_nr99_v138.2_toSpecies_trainset_uq.mmi ${sample}.fastq > ${sample}.sam

samtools view -b -o ${sample}.bam ${sample}.sam
samtools sort -o ${sample}.sorted.bam ${sample}.bam
samtools index ${sample}.sorted.bam

# Primary alignments only (drop secondary 0x100 and supplementary 0x800)
samtools view -h -F 0x900 ${sample}.sorted.bam > ${sample}.sorted.primary.sam
samtools view -b -o ${sample}.sorted.primary.bam ${sample}.sorted.primary.sam

# List of reference IDs (RNAME) from primary alignments
samtools view ${sample}.sorted.primary.bam | awk '{print $3}' > ${sample}.sorted.primary.alignedseqs.txt
```

> **Note:** If your tools are installed in specific locations (e.g., `/shares/...` or `../../../../tools/...`), replace `minimap2` and `samtools` above with those full paths.

---

## Run All Samples in a Loop

To avoid repeating commands for each sample (`S4 S5 S6 S10 S11 S12`), you can use:

```bash
# Demultiplex once (produces bamboo22.{name}.fastq)
conda activate cutadapt
cutadapt -g file:forward_barcodes.fasta -e 0.1 --rc -j 64 \
  -o bamboo22.{name}.fastq bamboo22.fastq.gz

# Define samples present after demultiplexing
SAMPLES=(S4 S5 S6 S10 S11 S12)

# Round 1–3 trimming + length/quality
for S in "${SAMPLES[@]}"; do
  cutadapt -g GGTAGTATATACAGAGAG    -e 0.1 --rc -j 64 -o bamboo22.${S}.b.fastq bamboo22.${S}.fastq
  cutadapt -g AGRGTTYGATYMTGGCTCAG  -e 0.1 --rc -j 64 -o bamboo22.${S}.c.fastq bamboo22.${S}.b.fastq
  cutadapt -g RGYTACCTTGTTACGACTT   -e 0.1 --rc -j 64 -o bamboo22.${S}.d.fastq bamboo22.${S}.c.fastq
  cutadapt --minimum-length 500 -l 1550 --quality-cutoff 20,20 -j 64 \
           -o bamboo22.${S}.e.fastq bamboo22.${S}.d.fastq

done

# QC
conda activate fastqc
mkdir -p fastqc_output
fastqc *.e.fastq -o fastqc_output/

# Mapping & primary extraction (ensure minimap2 and samtools are available)
for S in "${SAMPLES[@]}"; do
  sample="bamboo22.${S}.e"

  minimap2 -t 32 -ax map-ont silva_nr99_v138.2_toSpecies_trainset_uq.mmi ${sample}.fastq > ${sample}.sam

  samtools view -b -o ${sample}.bam ${sample}.sam
  samtools sort -o ${sample}.sorted.bam ${sample}.bam
  samtools index ${sample}.sorted.bam

  samtools view -h -F 0x900 ${sample}.sorted.bam > ${sample}.sorted.primary.sam
  samtools view -b -o ${sample}.sorted.primary.bam ${sample}.sorted.primary.sam
  samtools view ${sample}.sorted.primary.bam | awk '{print $3}' > ${sample}.sorted.primary.alignedseqs.txt

done
```

---

## Parameter Notes

- **Cutadapt demultiplexing with `file:` and `{name}`**  
  Using `-g file:forward_barcodes.fasta` with `-o bamboo22.{name}.fastq` creates one output per FASTA header (the header becomes `{name}`), which is ideal for barcode demultiplexing.  
- **`--rc`** checks both the given adapter and its reverse complement (useful for ONT reads where orientation can vary).  
- **Error rate `-e 0.1` (10%)** is commonly used for ONT adapters/primers, balancing sensitivity and specificity.  
- **`--minimum-length 500` and `-l 1550`** enforce a read length window suitable for your amplicon design; adjust to your target if different.  
- **Mapping preset `-ax map-ont`** is optimized for ONT data.  
- **Primary alignments** are extracted with `samtools view -F 0x900` (removes secondary `0x100` and supplementary `0x800`).

---

## Troubleshooting

- **Tool not found:** Ensure `minimap2`, `samtools`, `cutadapt`, and `fastqc` are available on your `PATH` or use absolute paths.  
- **Demultiplexing yields unexpected sample names:** Check FASTA headers in `forward_barcodes.fasta`; those become `{name}` in outputs.  
- **Very low retention after trimming:** Relax `-e` (e.g., `0.15`) or review adapter sequences and `--rc`. Confirm expected amplicon size vs. `--minimum-length` and `-l`.  
- **Reference index missing:** If you only have a FASTA (e.g., `ref.fasta`), build an index:  
  ```bash
  minimap2 -d silva_nr99_v138.2_toSpecies_trainset_uq.mmi ref.fasta
  ```

---

## Citations

Please cite the tools used in this pipeline:

- **Cutadapt** – Adapter trimming and demultiplexing  
  Martin, M. (2011). *Cutadapt removes adapter sequences from high-throughput sequencing reads.* EMBnet.journal, 17(1), 10–12.

- **Minimap2** – Long-read alignment  
  Li, H. (2018). *Minimap2: pairwise alignment for nucleotide sequences.* Bioinformatics, 34(18), 3094–3100.

- **Samtools** – Alignment processing  
  Li, H. et al. (2009). *The Sequence Alignment/Map format and SAMtools.* Bioinformatics, 25(16), 2078–2079.

- **FastQC** – Read quality control  
  Andrews, S. (2010). *FastQC: A quality control tool for high throughput sequence data.* (Babraham Bioinformatics)

---

## License

Add your chosen license here (e.g., MIT, GPL‑3.0). If this repository accompanies a paper, you may also include a citation to the manuscript.
