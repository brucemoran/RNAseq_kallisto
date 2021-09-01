#!/usr/bin/env nextflow

/* QC, trimming of paired- or single-end Illumina data
* fastqc, fastp; bbduk trim by quality into single fastq
* pseudoalign with kallisto
* metrics from picard suite, visualised by multiQC
*/

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                  RNASEQ with KALLISTO and RNAseqR PIPELINE
  -----------------------------------------------------------------------

    Usage:
    nextflow run brucemoran/RNAseq_kallisto \
              -profile Configuration profile (required: genome or sonic, or DIY) \
              --sampleCsv "path/to/sample.csv" \
              --cdna  "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
              --gtf "http://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz" \
              --stranded  "rf-stranded"

    Mandatory arguments:
        --sampleCsv       [file]  CSV format, header to include "sampleID,read1,read2" in case of paired data; no read2 for single-end

        --runID           [str]   string naming run and output files

        --gtf             [str]   gtf file (absolute path or URL)

      One of:
        --kallistoindex   [str]   suitable kallisto index
        OR
        --cdna            [str]   suitable cDNA fasta file URL

    Optional arguments:
        --email           [str]   email address to send RNAseqR, multiQC reports

        --stranded        [str]   if data is fr- (first read forward) or rf-stranded (first read reverse), or not stranded (default: "")

        --dupradar        [bool]  run dupRadar to determine PCR contaminants (default:TRUE)

        --metadataCsv     [file]  CSV format, header to include "sample" and any other covariates to use in DE analysis

        --metadataDesign  [str]   string indicating design formula for DE analysis; do not include "~ 0 +", only terms separated by "+" with last term as the covariate to use in DE (all contrasts are made)

        --controlRef      [str]   "reference" level of last covariate in metadataDesign (default: NULL)

        --genomePrefix    [str]   prefix for "_gene_ensembl" datasets found in biomaRt::listDatasets(biomaRt::useMart("ensembl"))\$dataset (default: "hsapiens")

        --msigdbrSpecies  [str]   species name from msigdbr::msigdbr_show_species() (default:"Homo sapiens")

        --msigdbCategory  [str]   category for MsigDB (one of c("H", paste0("C", 1:7))), see gsea-msigdb.org/gsea/msigdb/collections.jsp for details (default: "H")
    """.stripIndent()
}

if (params.help) exit 0, helpMessage()

//set outDir
params.outDir = "RNAseq_output"

//test kallistoindex
if(!params.kallistoindex){
  if(!params.cdna){
    exit 1, "Require --cdna or --kallistoindex input!"
  }
}

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

/* 0.00: Input using sample.csv
*/
sampleInputs = Channel.fromPath("$params.sampleCsv")

/* 0.00: find endedness
*/
process endedness {

  label 'low_mem'

  input:
  file(txt) from sampleInputs

  output:
  file('inputs.txt') into sampleCsvInput

  script:
  """
  COUNT=\$(head -n1 $txt | perl -ane '@s=split(/\\,/); print scalar(@s);')
  if [[ \$COUNT == 2 ]];then
    echo "sampleID,read1,read2" > "inputs.txt"
    tail -n+2 $txt | while read LINE; do
      echo \$LINE",inputs.txt" >> "inputs.txt"
    done
  else
    cat $txt > "inputs.txt"
  fi
  """
}

sampleCsvInput.splitCsv( header: true )
              .map { row -> [row.sampleID, file(row.read1), file(row.read2)] }
              .set { bbduking }

/* 1.0: Input trimming
 */
process bbduk {

  label 'high_mem'

  input:
  tuple val(sampleID), file(read1), file(read2) from bbduking

  output:
  tuple val(sampleID), val(paired), file('*bbduk.fastq.gz') into (kallistoing, fastping)

  script:
  paired = "${read2}" == null ? "single" : "paired"
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  if [[ ${paired} == "single" ]]; then
   reformat.sh ${taskmem} \
      in=${read1} \
      out="stdin.fastq" \
      tossjunk=T | \
   bbduk.sh ${taskmem} \
      int=t \
      in="stdin.fastq" \
      out=${sampleID}".bbduk.fastq.gz" \
      k=${params.bbdkmerx} \
      mink=${params.bbdmink} \
      hdist=1 \
      ktrim=r \
      trimq=${params.bbdqtrim} \
      qtrim=rl \
      maq=20 \
      ref=${params.bbmapAdapters} \
      tpe \
      tbo \
      stats=${sampleID}".bbduk.adapterstats.txt" \
      overwrite=T
  else

    reformat.sh ${taskmem} \
      in1=${read1} \
      in2=${read2} \
      out="stdout.fastq" \
      tossjunk=T | \
    bbduk.sh ${taskmem} \
      int=t \
      in="stdin.fastq" \
      out1=${sampleID}".R1.bbduk.fastq.gz" \
      out2=${sampleID}".R2.bbduk.fastq.gz" \
      k=${params.bbdkmerx} \
      mink=${params.bbdmink} \
      hdist=1 \
      ktrim=r \
      trimq=${params.bbdqtrim} \
      qtrim=rl \
      maq=20 \
      ref=${params.bbmapAdapters} \
      tpe \
      tbo \
      stats=${sampleID}".bbduk.adapterstats.txt" \
      overwrite=T
  fi
  } 2>&1 | tee ${sampleID}".bbduk.runstats.txt"
  """
}

/* 1.1: fastp QC
*/

process fastp {

  label 'med_mem'

  input:
  tuple val(sampleID), val(paired), file(reads) from fastping

  output:
  file('*.json') into fastp_multiqc

  script:
  def paired = "${reads}"[1] == null ? "single" : "paired"
  """
  if [[ ${paired} == "paired" ]]; then
    fastp -w ${task.cpus} \
      -j ${sampleID}.bbduk.fastp.json \
      --in1 ${reads}[0] \
      --in2 ${reads}[1]
  else
    fastp -w ${task.cpus} \
      -j ${sampleID}.bbduk.fastp.json \
      --in1 ${reads}
  fi
  """
}

/* 3.0: kallisto index
 */
process kallistondx {

  output:
  file('*.kallisto') into kallisto_index

  script:
  """
  if [[ ${params.kallistoindex} == "null" ]];then
    wget ${params.cdna}
    CDNA=\$(ls)
    kallisto index -i \$CDNA".kallisto" \$CDNA
  else
    ln -s ${params.kallistoindex} ${params.kallistoindex}".kallisto"
  fi
  """
}

/* 3.1: Kallisto
*/
process kallisto {

  label 'med_mem'

  publishDir "${params.outDir}/samples", mode: "copy"

  input:
  tuple sampleID, val(paired), file(reads) from kallistoing
  file(kallistoindex) from kallisto_index

  output:
  file("${sampleID}") into de_kallisto
  file("${sampleID}/kallisto/${sampleID}.kallisto.log.txt") into kallisto_multiqc
  tuple val(paired), file("${sampleID}.bam") into dupradar

  script:
  def stranding = params.stranded == "" ? "" : "--${params.stranded}"
  """
  mkdir -p ${sampleID}/kallisto
  {
  COUNT=\$(ls | grep fastq | wc -l)
  if [[ \$COUNT != 1 ]]; then
    kallisto quant \
      -t 10 \
      -b 100 \
      -i ${kallistoindex} \
      --pseudobam \
      -o ${sampleID}/kallisto ${stranding} ${reads} | samtools view -Shb - > ${sampleID}.bam

  else
    kallisto quant \
      -t 10 \
      -b 100 \
      --single \
      -l 200 \
      -s 30
      -i ${kallistoindex} \
      --pseudobam \
      -o ${sampleID}/kallisto ${stranding} ${reads} | samtools view -Shb - > ${sampleID}.bam
  fi
  } 2>&1 | tee > ${sampleID}/kallisto/${sampleID}".kallisto.log.txt"
  """
}

/* 3.11: dupRadar
*/
process dupRadar {

  label 'high_mem'

  publishDir "${params.outDir}", mode: "copy"

  input:
  tuple val(paired), file(bam) from dupradar.collect()

  output:
  file("dupradar") into dupradared

  when:
  !params.dupradar

  script:
  """
  {
  if [[ ${params.stranded} == "" ]]; then
    strand=0
  else
    if [[ ${params.stranded} == "fr-stranded" ]]; then
      strand=1
    else
      strand=2
    fi
  fi

  ##test GTF is URL?
  GTFURL=\$(curl --head --silent ${params.gtf} | head -n 1 | grep "OK" | wc -l)
  if [[ \$GTFURL =~ 1 ]]; then
    wget -o use.gtf ${params.gtf}
    gtf="use.gtf"
  else
    gtf=${params.gtf}
  fi

  Rscript --vanilla ${params.dupradarRscript} \
    \$gtf \
    ${paired} \
    \$strand \
    ${task.cpus}
  } 2>&1 | tee > dupradar/dupradar.log.txt"
  """
}

/* 3.0: MultiQC
 * NB bbduk_multiqc output not available for multiqc yet 130718
*/
fastp_multiqc.mix( kallisto_multiqc ).set { multiqc_multiqc }

process mltiqc {

  label 'low_mem'

  publishDir "${params.outDir}/multiqc", mode: "copy", pattern: "*"

  input:
  file ('*') from multiqc_multiqc.collect()

  output:
  file('*') into completedmultiqc
  file("${params.runID}.multiqc_report.html") into sendmail_multiqc

  script:
  """
  multiqc -n ${params.runID}.multiqc_report.html -c ${params.multiqcConfig} -f .
  """
}

/* 4.0: DE analysis
*/

process RNAseqon {

  label 'low_mem'

  publishDir "${params.outDir}", mode: "copy", pattern: "*[!.zip]"
  publishDir "${params.outDir}/${params.runID}_RNAseqon", mode: "copy", pattern: "*[.zip]"

  input:
  file (kdirs) from de_kallisto.collect()

  output:
  file("${params.runID}_RNAseqon") into completedRNAseqR
  file("${params.runID}_RNAseqon.zip") into sendmail_RNAseqR

  when:
  params.metadataCsv && params.metadataDesign

  script:
  def control_ref = params.controlRef == null ? "NULL" : params.controlRef
  def genome_pref = params.genomePrefix == null ? "hsapiens" : params.genomePrefix
  """
  Rscript -e "RNAseqon::run_prep_modules_bm(metadata_csv = \\"${params.metadataCsv}\\", metadata_design = \\"${params.metadataDesign}\\", tag = \\"${params.runID}\\", output_dir = \\"${params.runID}_RNAseqon\\", data_dir = NULL, control_reference = \\"${control_ref}\\", genome_prefix = \\"${genome_pref}\\", msigdb_species = \\"${params.msigdbSpecies}\\")"
  rm Rplots.pdf
  zip -r ${params.runID}_RNAseqon.zip ${params.runID}_RNAseqon
  """
}

// 4.1: ZIP for sending on sendmail
// sendmail_RNAseqR.mix(sendmail_multiqc).set { sendmail_all }
//
// process zipup {
//
//   label 'low_mem'
//
//   input:
//   file(send_all) from sendmail_all.collect()
//
//   output:
//   file("${params.runID}.RNAseq_kallisto.zip") into send_zip
//
//   script:
//   """
//   zip ${params.runID}.RNAseq_kallisto.zip *
//   """
// }

// 5.0: Completion e-mail notification
if(params.email){
  workflow.onComplete {
    sleep(1000)
    def subject = """\
      [brucemoran/RNAseq_kallisto] SUCCESS: ${params.runID} [${workflow.runName}]
      """
      .stripIndent()
    if (!workflow.success) {
        subject = """\
          [brucemoran/RNAseq_kallisto] FAILURE: ${params.runID} [${workflow.runName}]
          """
          .stripIndent()
    }

    def msg = """\
      Pipeline execution summary
      ---------------------------
      RunID       : ${params.runID}
      RunName     : ${workflow.runName}
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """
      .stripIndent()

    def attachments = sendmail_multiqc.toList().getVal()

    sendMail(to: "${params.email}",
             subject: subject,
             body: msg,
             attach: attachments)
  }
}
