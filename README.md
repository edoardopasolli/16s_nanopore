# Full-length 16S rRNA preprocessing and mapping pipeline  

This repository provides a reproducible workflow to preprocess full-length 16s rRNA nanopore sequencing (demultiplexing and adapter/primer trimming), run QC, and map reads to a reference database to extract primary alignments.

**Highlights**
- Demultiplex with barcode FASTA  
- Multi‑step 5' trimming (barcodes/primers) with reverse‑complement matching  
- Length & quality trimming  
- QC with FastQC  
- Mapping with minimap2 (`map-ont`) and primary-alignment extraction with samtools

---

## Contents
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Step-by-Step Pipeline](#step-by-step-pipeline)
- [Parameter Notes](#parameter-notes)
- [Citation](#citations)
- [License](#license)

---

## Requirements

Tested on Linux (x86_64). Recommended sufficient RAM for indexing/mapping.

**Software:**
- cutadapt
- FastQC
- minimap2
- samtools
- conda (for environment management)

**Conda install:**

```bash
conda create -y -n cutadapt cutadapt
conda create -y -n fastqc fastqc
conda create -y -n mapping minimap2 samtools
```

If you already have binaries (e.g. on an HPC), ensure they are on your `PATH` or use absolute paths.

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
  > Minimap2 prebuilt index

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

## Step-by-Step Pipeline

### 1) Demultiplex by barcode

```bash
conda activate cutadapt
cutadapt -g file:forward_barcodes.fasta -e 0.1 --rc -j 64 -o bamboo22.{name}.fastq bamboo22.fastq.gz
```

- `-g file:forward_barcodes.fasta` reads multiple adapters from a FASTA; `{name}` is replaced with the FASTA header names to produce one output per barcode.  
- `--rc` searches the reverse complement if present.  
- `-e 0.1` allows up to 10% error rate in matches.  
- **Outputs** like `bamboo22.S4.fastq`, `bamboo22.S5.fastq`, etc.

### 2) Multi‑step 5’ trimming (barcodes/primers)

For each sample, sequentially trim the three adapter/primer motifs (allowing RC and 10% error):

```bash
# Barcode (round 1)
cutadapt -g GGTAGTATATACAGAGAG -e 0.1 --rc -j 64 -o bamboo22.S4.b.fastq bamboo22.S4.fastq

# Primer (round 2)
cutadapt -g AGRGTTYGATYMTGGCTCAG -e 0.1 --rc -j 64 -o bamboo22.S4.c.fastq bamboo22.S4.b.fastq

# Primer (round 3)
cutadapt -g RGYTACCTTGTTACGACTT -e 0.1 --rc -j 64 -o bamboo22.S4.d.fastq bamboo22.S4.c.fastq
```

### 3) Length and quality trimming

```bash
cutadapt --minimum-length 500 -l 1550 --quality-cutoff 20,20 -j 64 -o bamboo22.S4.e.fastq bamboo22.S4.d.fastq
```

- `--minimum-length 500`: retain reads ≥ 500 nt  
- `-l 1550`: crop reads to at most 1550 nt (3` end clipping to a fixed length)  
- `--quality-cutoff 20,20`: trim low-quality bases from both ends (Phred 20)

### 4) Quality control

```bash
conda activate fastqc
mkdir -p fastqc_output
fastqc *.e.fastq -o fastqc_output/
```

### 5) Mapping and primary alignments

```bash
conda activate mapping
sample=bamboo22.S4.e

minimap2 -t 32 -ax map-ont silva_nr99_v138.2_toSpecies_trainset_uq.mmi ${sample}.fastq > ${sample}.sam

samtools view -b -o ${sample}.bam ${sample}.sam
samtools sort -o ${sample}.sorted.bam ${sample}.bam
samtools index ${sample}.sorted.bam

# Primary alignments only (drop secondary 0x100 and supplementary 0x800)
samtools view -h -F 0x900 ${sample}.sorted.bam > ${sample}.sorted.primary.sam
samtools view -b -o ${sample}.sorted.primary.bam ${sample}.sorted.primary.sam

# List of reference IDs (RNAME) from primary alignments
samtools view ${sample}.sorted.primary.bam | awk `{print $3}` > ${sample}.sorted.primary.alignedseqs.txt
```

## 6) Filter high‑quality primary alignments

```bash
conda activate mapping
sample=bamboo22.S4.e

samtools view -h -q 20 -F 0x904 "${sample}.sorted.bam" | \
awk '
  BEGIN{OFS="\t"}
  /^@/ { print; next }  # pass header lines
  {
    cigar = $6
    # compute aligned length from CIGAR (sum of M,=,X,I,D; ignore S/H/P/N)
    aln = 0
    tmp = cigar
    while (match(tmp, /([0-9]+)([MIDNSHP=X])/, a)) {
      n = a[1]
      op = a[2]
      if (op=="M" || op=="=" || op=="X" || op=="I" || op=="D") aln += n
      tmp = substr(tmp, RSTART + RLENGTH)
    }
    # parse NM tag (edit distance)
    nm = 0
    if (match($0, /NM:i:([0-9]+)/, m)) nm = m[1]
    pid = (aln > 0) ? (1 - nm/aln) : 0
    # filters: aligned length >= 500 bp AND PID >= 0.85
    if (aln >= 500 && pid >= 0.85) print
  }' | \
samtools view -b -o "${sample}.filtered.q20.pid85.len500.bam" -

samtools index "${sample}.filtered.q20.pid85.len500.bam"

samtools view "${sample}.filtered.q20.pid85.len500.bam" | awk '{print $3}' > "${sample}.filtered.q20.pid85.txt"
```

---

## Parameter Notes

- **Cutadapt demultiplexing with `file:` and `{name}`**  
  Using `-g file:forward_barcodes.fasta` with `-o bamboo22.{name}.fastq` creates one output per FASTA header (the header becomes `{name}`), which is ideal for barcode demultiplexing.  
- **`--rc`** checks both the given adapter and its reverse complement (useful for ONT reads where orientation can vary).  
- **Error rate `-e 0.1` (10%)** can be used for ONT adapters/primers, balancing sensitivity and specificity.  
- **`--minimum-length 500` and `-l 1550`** enforce a read length window suitable for your amplicon design; adjust to your target if different.  
- **Mapping preset `-ax map-ont`** is optimized for ONT data.  
- **Primary alignments** are extracted with `samtools view -F 0x900` (removes secondary `0x100` and supplementary `0x800`).

---

## Citations

Please cite the pipeline as follows:

- Ida Romano, Edoardo Pasolli, Jean-Claude Walser, Valeria Ventorino, Sonja Reinhard, Giuseppina Magaraci, Olimpia Pepe, Natacha Bodenhausen. *A hybrid and cost-efficient barcoding strategy for full-length 16s rRNA.* Under review.
  
---

## License

This repository is licensed under **GNU General Public License v3.0 or later (GPL‑3.0‑or‑later)**. This is a strong *copyleft* license that allows use, modification, and redistribution, provided that any distributed derivative works are also licensed under GPL‑3.0 (or a later GPL version), the source code is provided, and the same freedoms are preserved. It also includes protections against *tivoization* (hardware lockdown that prevents running modified versions) and contains explicit patent provisions; the software is provided **without warranty**.
