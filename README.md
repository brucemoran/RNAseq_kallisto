# Lightweight, Low-overhead RNAseq Data Processing
### A [NextFlow](https://nextflow.io) pipeline to pseudoalign and determine abundances from RNAseq data using [Kallisto](https://pachterlab.github.io/kallisto/)

### About the pipeline
The [nf-core RNAseq pipeline](https://github.com/nf-core/rnaseq) does not offer Kallisto as an option. This pipeline only offers Kallisto.

### How to Setup
Install [NextFlow](https://www.nextflow.io/index.html#GetStarted) and [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)

### Running the pipeline
```
nextflow run brucemoran/RNAseq_kallisto \
  -profile genome or sonic or DIY \
  --sampleCsv sample.rnaseq.csv \
  --cdna 'ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' \
  --stranded rf-stranded
```

### Arguments
```
nextflow run brucemoran/RNAseq_kallisto --help

Mandatory arguments:
  --sampleCsv     FILE      CSV format, header to include "sampleID,read1,read2" in case of paired data; no read2 for single-end
  --kallistoindex     STRING      suitable kallisto index
  OR
  --cdna     STRING      suitable cDNA fasta file URL

  Optional arguments (must include one!):
  --stranded     STRING     'fr-stranded' (first read forward) or 'rf-stranded' (first read reverse);
```
N.B. `rf-stranded` is [typical for Illumina stranded protocols](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/)

### Format of --sampleCsv:
```
sampleID,read1,read2
germ1,/full/path/to/germ1.R1.fastq.gz,/full/path/to/germ1.R2.fastq.gz
soma1,/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
soma2,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
```
N.B. headers of `sampleCsv` argument must match above exactly. For each `sampleID`, an output directory is written to `./RNAseq_output` with a `kallisto` directory (tagged with `stranded` argument if used) and the output of Kallisto therein.
