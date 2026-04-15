# Practical 1 — From FASTQ to SNPs: a resequencing workflow

This practical is a **basic end-to-end short-read resequencing workflow**:

1. inspect FASTQ files
2. run QC with FastQC and MultiQC
3. trim reads with fastp
4. map reads to the regional reference
5. check alignment depth and coverage
6. call variants per sample with GATK HaplotypeCaller
7. joint-genotype all samples with GenomicsDBImport and GenotypeGVCFs
8. optionally inspect the provided regional VCF


This tutorial uses a **small subset of a real WGS dataset** from the species **Biscutella laevigata**, sampled from **lowland** and **high-elevation** populations.

## What you are given in this tutorial

The tutorial dataset provides:

- **sixteen paired-end samples** (i.e., 8 samples from high-elevation (A2-S) and 8 from lowland (A2-G) populations) x paired end R1/R2 FASTQ files
- **One reference genome** (FASTA)
- A **gene annotation** file (GFF/GFF3)
- A **BED file of 4-fold degenerate sites** (so you don’t need to generate them from scratch)

Suggested folder layout (adjust names to whatever you were given):

```text
dataset-tutorial1/reads/                 # FASTQ(.gz) for the sixteen samples
dataset-tutorial1/ref/genome.fa          # reference FASTA
dataset-tutorial1/ref/genes.gff3         # annotation
dataset-tutorial1/ref/4fold.bed          # provided 4-fold sites (BED)
```

## Sample details

At the start of the practical, identify the sample names directly from the FASTQ files.

Each sample should have:

- one `R1` file
- one `R2` file
- the same sample name before `_R1` or `_R2`

For example, if you have:

- `reads/sample1_R1.fastq.gz`
- `reads/sample1_R2.fastq.gz`

then the sample name is `sample1`.


List all samples:

```bash
ls dataset-tutorial1/reads/*_R1.fastq.gz
```

If you want just the sample names:

```bash
for r1 in dataset-tutorial1/reads/*_R1.fastq.gz; do
  basename "$r1" _R1.fastq.gz
done
```

The same sample-name pattern is used throughout the tutorial loops below.

## VS Code defaults in the Renku image

In this Renku image, `Ctrl+Enter` is set to run `Terminal: Run Selected Text in Active Terminal`.

If that shortcut stops working in an existing session, run the same command once from the command palette or restart the Renku session.

---

## 1) FASTQ: what it is

A FASTQ file stores sequencing reads:

- **Sequence** (A/C/G/T/N)
- **Quality scores** per base (Phred scores)

The Phred score $Q$ relates to the probability of an incorrect base call $p$:

$$Q = -10 \log_{10}(p)$$

So:

- $Q20 \Rightarrow p \approx 1\%$
- $Q30 \Rightarrow p \approx 0.1\%$

Paired-end reads usually come as:

- `sample_R1.fastq.gz`
- `sample_R2.fastq.gz`


### Inspecting a FASTQ file

Before running any analysis, it's good practice to open a FASTQ file and check its contents to understand the format and verify the data. You can use `less`, `head`, or any text editor to view the file. For example:


```bash
zcat dataset-tutorial1/reads/sample1_R1.fastq.gz | head
```

Look for the four-line FASTQ structure:  
1. Sequence identifier (starts with `@`)  
2. Sequence (A/C/G/T/N)  
3. Separator line (starts with `+`)  
4. Quality scores (ASCII characters)  

---

## 2) QC of raw FASTQ (FastQC + MultiQC)

### Why

Before trimming/mapping, you want to know whether problems are:

- **Global** (bad sequencing run) or **sample-specific** (library issues)
- Fixable by trimming or filtering

### Commands

For one file:

```bash
mkdir -p results/qc/fastqc
fastqc -o results/qc/fastqc dataset-tutorial1/reads/sample_R1.fastq.gz
```

For all files:

```bash
mkdir -p results/qc/fastqc
fastqc -o results/qc/fastqc -t 4 dataset-tutorial1/reads/*.fastq.gz
multiqc -o results/qc results/qc/fastqc
```

### Key FastQC plots to understand

- **Per base sequence quality**: falling quality at 3’ ends is common. Official FastQC help: [Per Base Sequence Quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)
- **Adapter content**: adapters indicate reads are longer than inserts. Official FastQC help: [Adapter Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)
- **Per sequence GC content**: odd shapes can indicate contamination. Official FastQC help: [Per Sequence GC Content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html)
- **Overrepresented sequences**: adapters, primers, rRNA, contamination. Official FastQC help: [Overrepresented Sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)
- **Sequence duplication levels**: can be normal (high depth) or PCR artifacts. Official FastQC help: [Sequence Duplication Levels](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)

### Red flags

- Very low median quality across read length (not just 3’ tail)
- Strong adapter content early in reads
- GC profile wildly different from expectation (contamination)
- Many reads very short *before* any trimming (library degradation)

---

## 3) Trimming and filtering (fastp)

### Why trimming is needed

Trimming/filters can:

- Remove **adapter sequences**
- Remove **low-quality tails** that cause mismatches in mapping
- Filter out **very short** reads that map ambiguously

This typically improves:

- Mapping rate
- Variant calling specificity
- Downstream QC consistency

### Example command

For one sample:

```bash
mkdir -p results/trimmed results/unpaired results/fastp_reports

fastp \
  -i dataset-tutorial1/reads/sample_R1.fastq.gz \
  -I dataset-tutorial1/reads/sample_R2.fastq.gz \
  -o results/trimmed/sample_R1.trim.fastq.gz \
  -O results/trimmed/sample_R2.trim.fastq.gz \
  --unpaired1 results/unpaired/sample_R1.unpaired.fastq.gz \
  --unpaired2 results/unpaired/sample_R2.unpaired.fastq.gz \
  --detect_adapter_for_pe \
  --cut_front --cut_tail \
  --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 50 \
  --thread 4 \
  --html results/fastp_reports/sample.fastp.html \
  --json results/fastp_reports/sample.fastp.json
```

For all samples:

```bash
mkdir -p results/trimmed results/unpaired results/fastp_reports

for r1 in dataset-tutorial1/reads/*_R1.fastq.gz; do
  sample=$(basename "$r1" _R1.fastq.gz)
  r2="dataset-tutorial1/reads/${sample}_R2.fastq.gz"

  echo "Trimming ${sample}"

  fastp \
    -i "$r1" \
    -I "$r2" \
    -o "results/trimmed/${sample}_R1.trim.fastq.gz" \
    -O "results/trimmed/${sample}_R2.trim.fastq.gz" \
    --unpaired1 "results/unpaired/${sample}_R1.unpaired.fastq.gz" \
    --unpaired2 "results/unpaired/${sample}_R2.unpaired.fastq.gz" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail \
    --cut_window_size 4 --cut_mean_quality 20 \
    --length_required 50 \
    --thread 4 \
    --html "results/fastp_reports/${sample}.fastp.html" \
    --json "results/fastp_reports/${sample}.fastp.json"
done

multiqc -o results/fastp_reports results/fastp_reports
```

### How to read fastp output

fastp produces an **HTML** (human) and **JSON** (machine) report.

Common metrics:

- **Reads passed / failed** filters
- **Adapter trimming** counts
- **Base quality before vs after**
- **Read length distribution** after trimming
- **Duplication rate** (rough indicator)

In this example, paired reads stay in `trimmed/`, while reads whose mate is removed during trimming are written to `unpaired/`.

### Red flags

- Large fraction of reads become very short (e.g., most < 50 bp)
- Enormous read loss due to quality filters
- Duplication rate extremely high (PCR duplicates / low-complexity library)

---

## 4) Mapping to a reference (BWA-MEM)

### Why

Mapping assigns reads to coordinates on a reference genome so you can:

- Compute coverage
- Call variants
- Compare individuals at homologous sites

### Basic workflow

1. Index the reference
2. Align
3. Convert to BAM, sort, index
4. Mark duplicates (recommended for variant calling)

Run for all samples:

```bash
REF=dataset-tutorial1/ref/ref.fa
mkdir -p results/ref results/bam
cp "$REF" results/ref/ref.fa
REF=results/ref/ref.fa

bwa index "$REF"

for r1 in results/trimmed/*_R1.trim.fastq.gz; do
  sample=$(basename "$r1" _R1.trim.fastq.gz)
  r2="results/trimmed/${sample}_R2.trim.fastq.gz"

  echo "Mapping ${sample}"

  bwa mem -t 8  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    "$REF" "$r1" "$r2" \
    | samtools sort -@ 4 -o "results/bam/${sample}.sorted.bam"

  samtools index "results/bam/${sample}.sorted.bam"
done
```

### Mark duplicates (Picard)

When this should not be done, or should be treated cautiously:

- do not remove duplicates blindly in ddRAD, GBS, or other reduced-representation datasets, where many reads naturally start at the same positions
- do not treat duplicate counts as a simple artifact signal in very high-depth targeted sequencing without checking the library design first
- do not remove duplicates blindly in RNA-seq, where highly expressed transcripts can naturally produce many reads from the same regions and duplicate counts actually  reflect biology not just technical artifacts.

For all samples:

```bash
for bam in results/bam/*.sorted.bam; do
  sample=$(basename "$bam" .sorted.bam)

  echo "Marking duplicates for ${sample}"

  picard MarkDuplicates \
    I="$bam" \
    O="results/bam/${sample}.sorted.markdup.bam" \
    M="results/bam/${sample}.dup_metrics.txt" \
    CREATE_INDEX=true
done

```

Check the format of the bam files (or uncompressed sam files):

```bash

samtools view results/bam/sample1.sorted.markdup.bam | less

```

### Red flags

- Low overall mapping rate
- Very high fraction of reads with MAPQ=0 (multi-mapping)
- High duplicate rate (PCR artifacts)

---

## 5) Coverage report

### Why coverage matters

Coverage affects:

- Power to call variants
- Confidence (depth, allele balance)
- Missing data rate

You usually care about:

- **Mean depth**
- **Breadth of coverage** (fraction of genome covered at ≥1× / ≥10×)
- **Uniformity** (how even coverage is across the genome)

### Quick depth summaries

For all samples:

```bash
for bam in results/bam/*.sorted.markdup.bam; do
  sample=$(basename "$bam" .sorted.markdup.bam)

  echo "Summarizing coverage for ${sample}"

  # Per-base depth (filtered)
  sambamba depth base -t 4 -c 0 \
    -F "mapping_quality > 0 and not (duplicate or failed_quality_control or unmapped or secondary_alignment)" \
    -o "results/bam/${sample}.depth.txt" "$bam"

  # Mean depth
  echo -e "mean_depth\t$(awk '{sum+=$3; n++} END{print (n?sum/n:0)}' results/bam/${sample}.depth.txt)" \
    > "results/bam/${sample}.depth.summary.txt"

  # Breadth >0x
  echo -e "%_covered_gt0x\t$(awk '{c++; if($3>0) total+=1} END{print (total/c)*100}' results/bam/${sample}.depth.txt)" \
    >> "results/bam/${sample}.depth.summary.txt"

  # Breadth >10x
  echo -e "%_covered_gt10x\t$(awk '{c++; if($3>10) total+=1} END{print (total/c)*100}' results/bam/${sample}.depth.txt)" \
    >> "results/bam/${sample}.depth.summary.txt"

  # BAM stats
  bamtools stats -in "$bam" | grep -v "*" > "results/bam/${sample}_bamtools_stats.txt"
  samtools flagstat "$bam" > "results/bam/${sample}_flagstat.txt"


done
```

### Red flags

- Very uneven coverage (peaks/holes) not explained by biology
- Large uncovered regions (could be contamination, poor mapping, reference issues)
- Samples with systematically lower depth than others (library prep / sequencing imbalance)

---

## 6) Variant calling and joint genotyping

### Note: different dataset types need different expectations

The steps after mapping (coverage summaries, variant calling, joint genotyping, filtering, annotation) are written mostly with **WGS (whole-genome short-read resequencing)** in mind.

If your data are *not* WGS, the same tools can often be used, but you should adjust expectations and sometimes change the workflow:

- **ddRAD / GBS (reduced-representation) mapped to a reference**
  - Coverage is concentrated at restriction-site loci; many genomic positions have **zero coverage** by design.
  - Expect **high missingness across loci** (between individuals) and more variable depth between loci.
  - Duplicate handling is trickier: some apparent duplicates can be genuine sampling of the same locus; PCR duplicates can also be common. 
  - Variant filtering often needs stronger missing-data and depth-based rules (e.g., per-site missingness thresholds; per-genotype depth thresholds) because “low depth” is more frequent.
  - Many ddRAD/GBS projects use purpose-built toolchains (e.g., Stacks/ipyrad) for locus construction and SNP calling; reference-based calling can work, but interpretation of coverage and missingness differs from WGS.

- **Pool-seq (sequencing pooled individuals)**
  - You typically estimate **allele frequencies** rather than individual genotypes; standard diploid genotype calling and joint genotyping are often not appropriate.
  - ANGSD-style likelihood-based summaries or pileup/frequency-based pipelines are often a better fit than per-individual genotype calling.

- **Target capture / exome sequencing**
  - Compute and watch **on-target rate** and coverage uniformity across targets; off-target reads can be substantial.
  - Use intervals/targets for variant calling and QC (you care about depth on targets, not genome-wide breadth).


- **Low-coverage resequencing**
  - Hard genotype calls become much less reliable when coverage per sample is low, because distinguishing `0/0`, `0/1`, and `1/1` is harder with only a few reads.
  - Expect more uncertainty, more missing data, and noisier allele-balance patterns.
  - In these cases, a genotype-likelihood framework is often better than forcing confident genotype calls too early.
  - A common recommendation is an **ANGSD-like pipeline**, where you work with genotype likelihoods, site-frequency spectra, allele frequencies, PCA, and population-genetic summaries without overcommitting to hard genotype calls.

- **RNA-seq**
  - Do not use BWA-MEM for spliced transcripts; use a splice-aware aligner (STAR/HISAT2). Variant calling from RNA-seq requires additional care (expression, splicing, allele-specific expression).

### Concepts

- **Variant calling**: identify positions where sample differs from reference.
- **Genotyping**: assign genotype (0/0, 0/1, 1/1) and confidence.
- **Joint genotyping**: genotype multiple samples together for consistency.

Joint genotyping helps:

- analyze all samples against the same set of candidate variant sites, instead of letting each sample decide on its own which sites exist
- make genotype calls more comparable across samples, which is important when you later compare populations or calculate allele frequencies
- reduce cases where one sample has a variant call but another sample is recorded as if nothing happened there, when the real issue is only lower confidence or lower coverage or borderline evidence

For this practical, after you create and check the BAM files, this is the main next step.

### GATK-style gVCF workflow (common)

Create the indexed reference FASTA:
```bash
REF=results/ref/ref.fa
samtools faidx "$REF"
gatk CreateSequenceDictionary -R "$REF" -O "${REF%.fa}.dict"
```

Then call variants per sample with HaplotypeCaller in GVCF mode:

```bash
mkdir -p results/vcf

REF=results/ref/ref.fa

for bam in results/bam/*.sorted.markdup.bam; do
  sample=$(basename "$bam" .sorted.markdup.bam)

  echo "Calling variants for ${sample}"

  gatk HaplotypeCaller \
    -R "$REF" \
    -I "$bam" \
    -O "results/vcf/${sample}.g.vcf.gz" \
    --all-site-pls \
    -ERC GVCF
done
```

Then jointly genotype with **GenomicsDB**:

```bash
mkdir -p results/gendb  results/tmp

# Import gVCFs into a GenomicsDB workspace
gatk GenomicsDBImport \
  -R "$REF" \
  $(printf -- '-V %s ' results/vcf/*.g.vcf.gz) \
  -L $(grep "^>" results/ref/ref.fa | sed 's/>//') \
  --genomicsdb-workspace-path results/gendb/cohort_db \
  --tmp-dir results/tmp \
  --reader-threads 4

# Jointly genotype from the GenomicsDB workspace
gatk GenotypeGVCFs \
  -R "$REF" \
  -all-sites \
  --standard-min-confidence-threshold-for-calling 20 \
  -V gendb://results/gendb/cohort_db \
  -O results/vcf/cohort.raw.vcf.gz

```

Notes:

- For real datasets, you typically add intervals (e.g. `-L chr1` or an intervals list) and/or split by contig to make GenomicsDBImport faster and more memory-friendly.
- For very small cohorts, `CombineGVCFs` can also work, but GenomicsDB scales better.


### Red flags

- Many variants with extreme strand bias metrics
- Huge number of variants in one sample only (sample contamination or mapping artifacts)
- Abnormally high heterozygosity across the genome (contamination, polyploidy, reference mismatch)



### Note: indel realignment (and why GATK HaplotypeCaller doesn’t need it)

Older variant-calling workflows often included an explicit **indel realignment** step after mapping. The reason is that read alignments around true indels can contain mismatches that look like clusters of false SNPs (and they can distort allele counts).

With **GATK HaplotypeCaller**, a separate indel realignment step is **not necessary**, because HaplotypeCaller performs **local re-assembly/realignment** of reads in “active regions” as part of calling (and the old `IndelRealigner` tool is deprecated/removed in GATK4). See the GATK docs describing HaplotypeCaller’s local assembly and the deprecation of IndelRealigner:

- https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
- https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847

For **Pool-seq** (or any workflow where you rely heavily on **read-count-based allele frequencies** from alignments/pileups), indel-associated misalignments can be especially painful because they:

- Create **false SNPs** near indels.
- Bias **allele frequency estimates** by shifting reads between reference/alternate support.

If you are not using an assembly-based caller (or if you’re doing Pool-seq allele-frequency analyses), it can be worth either (a) doing an indel-aware realignment step with an appropriate tool, or (b) using conservative filters such as masking SNPs within a window around called indels and applying strict mapping/base-quality thresholds.

---

## 7) Variant filtering (hard filters)

### Why filtering is needed

Raw VCFs include artifacts from:

- Mis-mapping in repeats
- Local alignment errors around indels
- Strand bias
- Low depth / low genotype quality

Filtering aims to remove variants that look unreliable.

### Hard filtering vs model-based filtering

- **Hard filtering**: apply thresholds (simple, transparent; good for small datasets)
- **VQSR** (model-based): requires large training sets; often not feasible for non-model organisms or small cohorts

### What to pay attention to

Per-variant and per-genotype fields:

- `DP` (depth)
- `GQ` (genotype quality)
- `AD` (allele depths; supports allele balance checks)
- `MQ` (mapping quality)
- Strand bias metrics (e.g., `FS`, `SOR`)

### Strategies for setting thresholds

Exact thresholds depend on coverage, and pipeline.

Use `bcftools` for exploring distributions:

```bash
mkdir -p results/plots
bcftools stats -s - results/vcf/cohort.raw.vcf.gz > results/vcf/cohort.stats.txt

# Explore depth distribution
bcftools query -f '[%DP\n]' results/vcf/cohort.raw.vcf.gz | \
python3 -c "
import sys, numpy as np, matplotlib.pyplot as plt
x = np.array([float(l.strip()) for l in sys.stdin if l.strip() != '.'])
plt.hist(x, bins=100, density=True)
plt.xlabel('DP'); plt.ylabel('Density')
plt.savefig('results/plots/depth_distribution.pdf')"

# Explore genotype quality distribution
bcftools query -f '[%GQ\n]' results/vcf/cohort.raw.vcf.gz | \
python3 -c "
import sys, numpy as np, matplotlib.pyplot as plt
x = np.array([float(l.strip()) for l in sys.stdin if l.strip() != '.'])
plt.hist(x, bins=100, density=True)
plt.xlabel('GQ'); plt.ylabel('Density')
plt.savefig('results/plots/genotype_quality_distribution.pdf')"

# Explore allele balance distribution for heterozygotes
python3 -c "import allel; import matplotlib.pyplot as plt; import numpy as np

callset = allel.read_vcf('results/vcf/cohort.raw.vcf.gz', fields=['AD', 'DP','GT', 'numalt'])


gt = callset['calldata/GT']
ad = callset['calldata/AD']

numalt = callset['variants/numalt']

# biallelic sites only
biallelic = numalt == 1

gt = gt[biallelic]
ad = ad[biallelic]

# heterozygotes only
het = (gt[:, :, 0] != gt[:, :, 1])

# calculate allele balance (alt / (ref + alt) for heterozygotes)
denom = ad.sum(axis=2).astype(float)

mask = het & (denom > 0)

ab = ad[:, :, 1].astype(float)[mask] / denom[mask]

plt.hist(ab, bins=100)
plt.xlabel('Allele Balance (heterozygotes)')
plt.ylabel('Frequency')
plt.title('AB distribution (biallelic het sites)')
plt.savefig('results/plots/ab_distribution.pdf')"
```


A common strategy is:

1. Filter obvious low-quality variants (variant-level)
2. Then apply genotype-level masks (e.g., set low-depth genotypes to missing)

[Basic hard filter recommended by GATK best practices:](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

```bash
gatk VariantFiltration \
    -V results/vcf/cohort.raw.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O results/vcf/cohort.filtered.vcf.gz
```
Select only passing variants:

```bash
gatk SelectVariants \
  -V results/vcf/cohort.filtered.vcf.gz \
  --exclude-filtered \
  -O results/vcf/cohort.filtered.pass.vcf.gz

bcftools stats -s - results/vcf/cohort.filtered.pass.vcf.gz > results/vcf/cohort.filtered.pass.stats.txt
``` 

Then set low-depth genotypes to missing:

```bash
bcftools filter \
  -e 'FMT/DP<10 || FMT/DP>100' \
  -S . \
  -Oz -o results/vcf/cohort.filtered.depth_masked.vcf.gz \
  results/vcf/cohort.filtered.pass.vcf.gz
tabix -p vcf results/vcf/cohort.filtered.depth_masked.vcf.gz

bcftools stats -s - results/vcf/cohort.filtered.depth_masked.vcf.gz > results/vcf/cohort.filtered.depth_masked.stats.txt
```


Filter on missingness (e.g., remove sites with >20% missing genotypes):

```bash
bcftools view \
  -i 'F_MISSING<0.2' \
  -Oz -o results/vcf/cohort.filtered.miss.vcf.gz \
  results/vcf/cohort.filtered.depth_masked.vcf.gz
tabix -p vcf results/vcf/cohort.filtered.miss.vcf.gz

bcftools stats -s - results/vcf/cohort.filtered.miss.vcf.gz > results/vcf/cohort.filtered.miss.stats.txt
```

Select only biallelic SNPs + invariant sites (if needed):

```bash
bcftools view \
  -m2 -M2 \            # keep only biallelic sites
  -v snps,ref \       # include SNPs and invariant (reference) sites
  -Oz -o results/vcf/cohort.biallelic.with_invariant.vcf.gz \
  results/vcf/cohort.filtered.miss.vcf.gz

tabix -p vcf results/vcf/cohort.biallelic.with_invariant.vcf.gz

bcftools stats -s - results/vcf/cohort.biallelic.with_invariant.vcf.gz > results/vcf/cohort.biallelic.with_invariant.stats.txt
```

---

## 8) Variant annotation and extraction of 4-fold degenerate sites

### Conceptual overview

Variant annotation assigns functional consequences to each polymorphism relative to gene models. Tools such as snpEff classify variants into categories including:

* synonymous vs nonsynonymous
* stop gained / loss
* splice site disruption
* UTR, intronic, intergenic

This step is typically essential when the objective is to partition variation by functional class and infer selection.

A specific subset of annotated sites are **4-fold degenerate positions**, defined as codon positions where all possible nucleotide substitutions are synonymous. These sites are often used as a proxy for near-neutral evolution.

### Strategy in this workflow

In this tutorial, a precomputed BED file of 4-fold degenerate sites in the reference genomeis already provided. 
Here we will directly extracts variants overlapping these known 4-fold sites.

### Practical steps

1. Start from a filtered cohort VCF
2. Subset variants using the provided 4-fold BED

```bash
FOURFOLD_BED=ref/4fold.bed

bcftools view \
  -T "$FOURFOLD_BED" \
  -Oz -o results/vcf/cohort.filtered.4fold.vcf.gz \
  results/vcf/cohort.filtered.vcf.gz

tabix -p vcf results/vcf/cohort.filtered.4fold.vcf.gz

bcftools stats -s - results/vcf/cohort.filtered.4fold.vcf.gz > results/vcf/cohort.filtered.4fold.stats.txt
```

---
