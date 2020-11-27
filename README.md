# Lightweight, Low-overhead RNAseq Data Processing
### A [NextFlow](https://nextflow.io) pipeline to pseudoalign and determine abundances from RNAseq data using [Kallisto](https://pachterlab.github.io/kallisto/)

### About the pipeline
The [nf-core RNAseq pipeline](https://github.com/nf-core/rnaseq) does not offer Kallisto as an option. This pipeline only offers Kallisto.

### How to Setup
Install [NextFlow](https://www.nextflow.io/index.html#GetStarted) and [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)

### Running the pipeline
```
nextflow run brucemoran/RNAseq_kallisto \
          -profile <genome or sonic, or DIY> \
          --sampleCsv "path/to/sample.csv" \
          --cdna  "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
```

### Arguments
```
nextflow run brucemoran/RNAseq_kallisto --help

Mandatory arguments:
    --sampleCsv       [file]  CSV format, header to include "sampleID,read1,read2" in case of paired data; no read2 for single-end

    --runID           [str]   string naming run and output files

  One of:
    --kallistoindex   [str]   suitable kallisto index
    OR
    --cdna            [str]   suitable cDNA fasta file URL

Optional arguments:
    --email           [str]   email address to send RNAseqR, multiQC reports

    --stranded        [str]   if data is fr- (first read forward) or rf-stranded (first read reverse), or not stranded (default: "")

    --metadataCsv     [file]  CSV format, header to include "sample" and any other covariates to use in DE analysis

    --metadataDesign  [str]   string indicating design formula for DE analysis; do not include "~ 0 +", only terms separated by "+" with last term as the covariate to use in DE (all contrasts are made)

    --controlRef      [str]   "reference" level of last covariate in metadataDesign (default: NULL)

    --genomePrefix    [str]   prefix for "_gene_ensembl" datasets found in biomaRt::listDatasets(biomaRt::useMart("ensembl"))\$dataset (default: "hsapiens")

    --msigdbrSpecies  [str]   species name from msigdbr::msigdbr_show_species() (default:"Homo sapiens")

    --msigdbCategory  [str]   category for MsigDB (one of c("H", paste0("C", 1:7))), see gsea-msigdb.org/gsea/msigdb/collections.jsp for details (default: "H")
```

### Nota Benes
```
`rf-stranded` is [typical for Illumina stranded protocols](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/); Novogene RNAseq (standard) is unstranded, and so use default: ""

Sample format of --sampleCsv input file:

  sampleID,read1,read2
  germ1,/full/path/to/germ1.R1.fastq.gz,/full/path/to/germ1.R2.fastq.gz
  soma1,/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
  soma2,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz

headers of `sampleCsv` argument must match above exactly. For each `sampleID`, an output directory is written to `./RNAseq_output/samples` with a `kallisto` directory (tagged with `stranded` argument if used) and the output of Kallisto therein.

When running RNAseqR R package for DE, include both --metadataCsv with a 'sample' column
and --metadataDesign which must specify at least one column header from metadataCsv file.
For designs with more than one covariate, the covariate of interest is last, and the others
are joined by a '+', e.g. --metadataDesign "cov_1 + weight + cov_of_interest".
```
