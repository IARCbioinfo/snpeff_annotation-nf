# snpeff_annotation-nf
## Nextflow pipeline to annotate VCF files with SnpEff and dbSnp

This repository contains a Nextflow pipeline for annotating genetic variants in VCF (Variant Call Format) files. The pipeline processes input VCF files, performs various annotations, and generates a comprehensive annotation file.

## Prerequisites

Make sure you have the following dependencies installed before running the pipeline:

- [Nextflow](https://www.nextflow.io/)
- [conda](https://conda.io/projects/conda/en/latest/index.html)
- [dbSNP database](https://ftp.ncbi.nlm.nih.gov/snp/organisms)
- [dbNSFP database](https://pcingola.github.io/SnpEff/ss_dbnsfp/)

## Pipeline Overview

1. **FilterInputFiles:** Filters input VCF files using PLINK 2 to retain PASS variants with a maximum of 2 alleles.

2. **AnnotateWithRSID:** Annotates variants with RSID using SnpSift and the dbSNP database.

3. **AnnotateWithImpact:** Annotates variants with functional impact using snpEff and a specified reference genome.

4. **FullyAnnotateWithDbSNP:** Performs comprehensive annotation using SnpSift and dbNSFP database, including information on gene impact, gnomAD data, REVEL scores, ClinVar information, and more.

5. **ExtractFields:** Extracts relevant fields from the annotated VCF files and creates a tab-separated text file with a header for downstream analysis.

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/IARCbioinfo/snpeff_annotation-nf
   cd snpeff_annotation-nf
   ```

2. Adjust the `nextflow.config` file if necessary.

3. Run the pipeline with:

   ```bash
   nextflow run main.nf -profile conda
   ```

## Input

| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| `--output_foinput_folder_with_VCF_fileslder`    |  `${baseDir}/VCFs/`  | Folder containing `*vcf.gz` files |
 

## Parameters

  * #### Optional

| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| `--reference_genome`    |  `GRCh37.75`  | Folder containing `*vcf.gz` files |
| `--dbNSF_path`     |  `${baseDir}/dbNSFP4.1a.txt.gz` | [dbNSFP database](https://pcingola.github.io/SnpEff/ss_dbnsfp/) |
| `--dbSNP_path`    |  `${baseDir}/dbsnp150.vcf.gz`  |    [dbSNP database](https://ftp.ncbi.nlm.nih.gov/snp/organisms) |
| `--output_path`    |  `${baseDir}/output` |  Output folder |

## Output

The final annotated and extracted information will be available in the output directory as `full_annotation.txt`.

## Customization

- Adjust the memory requirements etc in the `nextflow.config` file.
- Customize the annotation processes in the `main.nf` script based on your specific requirements.

## Acknowledgments

- This pipeline utilizes various bioinformatics tools and databases, including PLINK, bcftools, SnpSift, snpEff, dbNSFP, and dbSNP.
