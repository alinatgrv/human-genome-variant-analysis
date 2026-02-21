# H+ / Human Genome Variant Analysis

## 0. Project Overview

This project focuses on the analysis of human SNP-based genotyping data with the following objectives:

1. Prepare and validate genotype data (VCF format, genome build compatibility).
2. Identify maternal (mtDNA) and paternal (Y chromosome) haplogroups.
3. Infer selected phenotypic traits from genotype data.
4. Annotate variants using public databases (ClinVar, GWAS Catalog, dbSNP).
5. Identify clinically relevant variants.
6. Propose a small set of hypothetical CRISPR-based genome edits supported by literature evidence.

The analysis is performed as a reproducible bioinformatics workflow with clear separation of raw data, intermediate files, and documented commands.

Genome build compatibility and data privacy considerations are strictly maintained throughout the analysis.

---

## 1. Dataset Description

Two independent SNP-array datasets were used:

1. **Genotek VCF dataset**
2. **23andMe raw genotype dataset (v4)**

Both datasets are based on reference human genome build 37 (GRCh37 / hg19).

### 1.1 Genotek Dataset

Dataset source: Genotek  
File format: VCF  
Reference genome: hg19 (GRCh37)

Genome build was determined from the VCF header:

```bash
head -n 10 data/intermediate/genotek.vcf
```

Header line: 

##reference=hg19

Total number of variants:
```bash
grep -v "^#" data/intermediate/genotek.vcf | wc -l
```
Result: 618255 variants


---

### 1.2 23andMe Dataset

Raw genotype file:  
`data/raw/SNP_raw_v4_Full_20170514175358.txt`

Header inspection confirmed:
- We are using reference human assembly build 37


The file was converted to VCF using PLINK v1.9:

```bash
plink --23file data/raw/SNP_raw_v4_Full_20170514175358.txt \
  --recode vcf --out data/intermediate/23andme_snps_clean \
  --output-chr MT --snps-only just-acgt
```

Summary statistics:

- Total SNPs: 595401

- Mitochondrial SNPs (MT): 3286

- Y chromosome SNPs: 2084

The presence of Y chromosome variants confirms a male genotype.

## 2. Sex Inference

Sex was inferred based on the presence of Y chromosome variants.

Command used:
```bash
grep -w "^chrY" data/intermediate/genotek.vcf | wc -l
```

Result: 3556 Y chromosome variants detected.

Conclusion: The individual is genetically male (XY karyotype).

## 3. Origins
### 3.1 Mitochondrial Haplogroup (Genotek dataset)

To determine haplogroups, mitochondrial (chrM) and Y chromosome variants were extracted from the full VCF.

Initially, variants were extracted using:
```bash
(
grep "^##" data/intermediate/genotek.vcf
grep "^#CHROM" data/intermediate/genotek.vcf
grep "^chrM" data/intermediate/genotek.vcf
) > data/intermediate/mtDNA_fixed.vcf

(
grep "^##" data/intermediate/genotek.vcf
grep "^#CHROM" data/intermediate/genotek.vcf
grep "^chrY" data/intermediate/genotek.vcf
) > data/intermediate/Y_fixed.vcf
```
Number of mitochondrial variants:

```bash
wc -l data/intermediate/mtDNA_fixed.vcf
```
Result: 986 variants

Mitochondrial variants were classified using HaploGrep 3  
(PhyloTree 17 – Forensic Update 1.2).

**Input file:**  
`data/intermediate/mtDNA_fixed.vcf`

**Results:**

- Haplogroup: **U5b2a4**
- Number of mtDNA variants: 986
- Tree overlap: 6.85%
- Reference mismatches: 757
- Strand flips: 274

**Interpretation:**

Haplogroup U5 belongs to one of the oldest mitochondrial lineages in Europe (estimated age ~30,000–35,000 years). It is associated with Upper Paleolithic hunter-gatherer populations and remains common in Northern and Eastern Europe today.

The detected subclade U5b2a4 suggests deep maternal ancestry linked to ancient European populations.

**Note:**  
The relatively low tree overlap is expected for SNP-array data, since genotyping chips do not provide full mitochondrial genome coverage.
**HaploGrep classification result:**

![HaploGrep result](assets/haplogrep.jpg)

### 3.2 Y-Chromosome Haplogroup (23andMe dataset)


Y haplogroup inference was performed using the MorleyDNA Y-SNP Subclade Predictor (Experimental tree).

The 23andMe raw genotype file was pre-processed using the Morley "extractFromAutosomal" tool, which extracts recognized Y-SNP markers from autosomal SNP-array data.

Preprocessing summary:

- 166 recognized positive Y-SNP mutations
- 733 recognized negative Y-SNP mutations
- 1086 no-calls

Result:

**Suggested haplogroup:** R1a1a  
**Subclade:** R1a-L168 (R1a-M17, R1a-M198)

Haplogroup R1a is a major Eurasian Y-chromosome lineage commonly observed in Eastern Europe, Central Asia, and parts of South Asia.

The resolution is limited by SNP-array coverage; deeper subclade determination would require high-coverage Y sequencing (e.g., BigY or whole-genome sequencing).

![Morley result](assets/morley.jpg)

## 4. ClinVar Annotation (23andMe dataset)

To assess potential clinical relevance of detected variants, the 23andMe-derived VCF (GRCh37) was intersected with the ClinVar VCF (GRCh37 build) downloaded from NCBI.

### 4.1 ClinVar Intersection

ClinVar VCF (GRCh37) was downloaded and indexed:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi
```

The 23andMe VCF was compressed and indexed:
```bash
bgzip -c data/intermediate/23andme_snps_clean.vcf > data/intermediate/23andme_snps_clean.vcf.gz
tabix -p vcf data/intermediate/23andme_snps_clean.vcf.gz
```
Variant positions were extracted and used to query ClinVar:
```bash
/usr/local/bin/bcftools query -f '%CHROM\t%POS\n' data/intermediate/23andme_snps_clean.vcf.gz \
  | sort -u > data/intermediate/23andme_sites.tsv

/usr/local/bin/bcftools view -R data/intermediate/23andme_sites.tsv \
  -O v -o data/intermediate/clinvar_hits.vcf \
  data/reference/clinvar.vcf.gz
```
Total ClinVar-matched records:
```bash
wc -l results/tables/clinvar_hits.tsv
```
Result: 69,522 ClinVar entries overlapping by genomic position.

### 4.2 Distribution of Clinical Significance (CLNSIG)

Top categories observed:

- Pathogenic: 17,411

- Benign: 13,302

- Uncertain significance: 11,090

- Likely pathogenic: 6,221

- Conflicting classifications: 4,811

- Pathogenic/Likely pathogenic: 4,740

### 4.3 Pathogenic / Likely Pathogenic Filtering

Non-conflicting pathogenic variants were extracted:
```bash
awk -F'\t' '
BEGIN{IGNORECASE=1}
{
  clnsig=$7
  if (clnsig ~ /pathogenic/ && clnsig !~ /conflict/) print
}' results/tables/clinvar_hits.tsv > results/tables/clinvar_pathogenic_or_likely.tsv
```

Result: 28,384 records.

Further prioritization was performed based on ClinVar review status (CLNREVSTAT).

High-confidence subset (reviewed by expert panel or practice guideline):

```bash
awk -F'\t' '
BEGIN{IGNORECASE=1}
{
  clnsig=$7; rev=$8
  if (clnsig ~ /pathogenic/ && clnsig !~ /conflict/ &&
      (rev ~ /reviewed_by_expert_panel/ || rev ~ /practice_guideline/))
    print
}' results/tables/clinvar_hits.tsv > results/tables/clinvar_pathogenic_expert_or_guideline.tsv
```
Result: 3,180 high-confidence records.

**Interpretation**

The large number of ClinVar matches reflects the fact that most SNP-array variants are represented in ClinVar, including benign and uncertain variants.

Further prioritization is required to:

identify clinically actionable variants

remove benign polymorphisms

focus on variants with strong review status and clear disease association

The next step involves structured prioritization and selection of a limited set of clinically relevant variants for detailed interpretation.

### 4.4 Allele-Level Annotation and Genotype Filtering

The initial ClinVar intersection was performed at the genomic position level, resulting in 69,522 overlapping records. However, positional overlap alone does not indicate that the individual carries the corresponding pathogenic allele.

To obtain clinically meaningful results, ClinVar annotations were incorporated directly into the 23andMe VCF at the allele level using bcftools annotate.
```bash
bcftools annotate \
  -a data/reference/clinvar.vcf.gz \
  -c CHROM,POS,REF,ALT,INFO \
  -O z -o data/intermediate/23andme_plus_clinvar.vcf.gz \
  data/intermediate/23andme_snps_clean.vcf.gz

tabix -p vcf data/intermediate/23andme_plus_clinvar.vcf.gz
```

Only variants with non-reference genotypes (GT ≠ 0/0) were retained for clinical prioritization.
```bash
bcftools view -i 'GT!="0/0" && INFO/CLNSIG!=""' \
  -O z -o results/vcf/23andme_clinvar_nonref.vcf.gz \
  data/intermediate/23andme_plus_clinvar.vcf.gz
  ```

#### Summary of Filtering Steps

| Step | Description | Number of Variants |
|------|------------|-------------------|
| 1 | Position-level ClinVar overlap | 69,522 |
| 2 | Allele-matched ClinVar annotations | 40,495 |
| 3 | Non-reference genotypes (GT ≠ 0/0) | 3,939 |
| 4 | P/LP (non-conflicting) | 12 |
| 5 | Expert panel / guideline | 0 |

This stepwise filtering reduced tens of thousands of positional matches to a small, interpretable set of clinically relevant candidate variants.

---

### 4.5 Pathogenic / Likely Pathogenic Variants (Non-reference Genotypes)

A total of 12 non-conflicting Pathogenic or Likely Pathogenic variants were identified in heterozygous state (0/1).

#### Distribution by ClinVar Review Status

| Review Status | Count |
|---------------|-------|
| Multiple submitters, no conflicts | 4 |
| Single submitter | 6 |
| No assertion criteria provided | 2 |

No variants were classified as reviewed by expert panel or practice guideline.

### 4.6 Prioritized Candidate Variants

The most reliable variants (multiple submitters, no conflicts) are listed below:

| Gene | Variant ID | Genotype | Clinical Association | Review Status |
|------|------------|----------|---------------------|---------------|
| CFH | rs460897 | 0/1 | Atypical HUS susceptibility / AMD | Multiple submitters, no conflicts |
| MYO7A | rs1052030 | 0/1 | Usher syndrome type 1 (AR) | Multiple submitters, no conflicts |
| VWF | i5049301 | 0/1 | Von Willebrand disease type 2 | Multiple submitters, no conflicts |
| CYP21A2 | i5005436 | 0/1 | Congenital adrenal hyperplasia (AR) | Multiple submitters, no conflicts |

All variants were detected in heterozygous state (0/1).

For autosomal recessive conditions (e.g., *MYO7A*, *CYP21A2*), the findings most likely represent carrier status rather than disease manifestation.

The *VWF* variant may be clinically relevant depending on phenotype and inheritance pattern and would require clinical-grade validation.

The *CFH* variant is associated with susceptibility phenotypes and does not imply deterministic disease development.

### 4.7 Variants with Limited or Conflicting Evidence

The remaining Pathogenic/Likely Pathogenic variants were supported by either single submitters or lacked assertion criteria in ClinVar.

These variants include genes associated with:

- SERPING1 (hereditary angioedema)
- MAP1A (neurodevelopmental disorder)
- SLC12A1 (renal tubular disorder, autosomal recessive)
- NAGLU (Mucopolysaccharidosis type IIIB, autosomal recessive)
- RGS9 (retinal dystrophy, autosomal recessive)
- ESR1 (estrogen resistance susceptibility)

Most of these findings were present in heterozygous state (0/1) and are consistent with potential carrier status for autosomal recessive disorders.

Two variants were annotated with limited or ambiguous evidence:

- UGT1A cluster variant (no assertion criteria provided)
- CDKN2B (Likely_pathogenic | protective), indicating mixed or conflicting interpretation

Due to limited evidence strength and lack of expert-panel review, these findings should be interpreted cautiously and require independent validation before any clinical consideration.

## 5. Hypothetical Genome Editing Strategy

To illustrate potential translational applications of genome analysis, a set of hypothetical CRISPR-based corrections was designed for selected Pathogenic/Likely Pathogenic variants identified in this dataset.

These examples are theoretical and intended solely to demonstrate conceptual feasibility of allele correction strategies [1–3].

---

### 5.1 VWF (chr12:6128067 G>A)

- Gene: VWF  
- Genotype: 0/1  
- Clinical association: Von Willebrand disease type 2  
- Pathogenic allele: A  
- Proposed correction: A → G (reference allele)

**Editing strategy:**  
A cytosine or adenine base editor could theoretically be used to revert the pathogenic nucleotide to the reference sequence without introducing double-strand breaks. Alternatively, homology-directed repair (HDR) could be employed in dividing cells [1].

**Limitations:**  
Clinical significance depends on variant-specific functional impact and inheritance pattern. Validation in a clinical sequencing setting and functional assays would be required.

---

### 5.2 CYP21A2 (chr6:32008312 C>T)

- Gene: CYP21A2  
- Genotype: 0/1  
- Clinical association: Congenital adrenal hyperplasia (autosomal recessive)  
- Pathogenic allele: T  
- Proposed correction: T → C

**Editing strategy:**  
Single-base correction via base editing [2,3] or HDR-based CRISPR system.

**Interpretation note:**  
As the condition is autosomal recessive and detected in heterozygous state, this represents likely carrier status rather than active disease.

---

### 5.3 MYO7A (chr11:76853783 T>A)

- Gene: MYO7A  
- Genotype: 0/1  
- Clinical association: Usher syndrome type 1 (autosomal recessive)  
- Pathogenic allele: A  
- Proposed correction: A → T

**Editing strategy:**  
Precision base editing targeting the mutant allele [2].

**Limitations:**  
Carrier state does not require therapeutic intervention; included here for conceptual demonstration only.

---

### 5.4 CFH (chr1:196716319 C>T)

- Gene: CFH  
- Genotype: 0/1  
- Clinical association: Atypical hemolytic uremic syndrome susceptibility  
- Pathogenic allele: T  
- Proposed correction: T → C

**Editing strategy:**  
Base editing may theoretically revert the susceptibility-associated allele [2,3].

**Important note:**  
CFH variants often act as susceptibility modifiers rather than deterministic mutations. Therapeutic editing would require strong functional evidence.

---

### 5.5 SERPING1 (chr11:57365745 T>C)

- Gene: SERPING1  
- Genotype: 0/1  
- Clinical association: Hereditary angioedema  
- Pathogenic allele: C  
- Proposed correction: C → T

**Editing strategy:**  
Allele-specific CRISPR correction with HDR template or base editing [1–3].

**Limitations:**  
Variant supported by single submitter; clinical relevance requires independent validation.

---

## References

1. Jinek M, et al. (2012). A programmable dual-RNA–guided DNA endonuclease in adaptive bacterial immunity. *Science*, 337(6096), 816–821. https://doi.org/10.1126/science.1225829

2. Komor AC, et al. (2016). Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. *Nature*, 533, 420–424. https://doi.org/10.1038/nature17946

3. Anzalone AV, et al. (2019). Search-and-replace genome editing without double-strand breaks or donor DNA. *Nature*, 576, 149–157. https://doi.org/10.1038/s41586-019-1711-4

4. Sadler JE. (2005). von Willebrand factor: two sides of a coin. *Journal of Thrombosis and Haemostasis*, 3(8), 1702–1709.

5. Speiser PW, et al. (2018). Congenital adrenal hyperplasia due to steroid 21-hydroxylase deficiency. *Endocrine Reviews*, 39(4), 389–434.

6. Weil D, et al. (1995). Defective myosin VIIA gene responsible for Usher syndrome type 1B. *Nature*, 374, 60–61.
