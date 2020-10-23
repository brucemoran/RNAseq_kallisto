#!/usr/bin/env nextflow

/* QC, trimming of paired- or single-end Illumina data
* fastqc, fastp; bbduk trim by quality into single fastq
* pseudoalign with kallisto
* metrics from picard suite, visualised by multiQC
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '-------------------------------------'
  log.info 'NEXTFLOW RNASEQ QC, TRIM, PSEUDOALIGN'
  log.info '-------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run brucemoran/RNAseq_kallisto \
            -profile Configuration profile (required: genome or sonic, or DIY) \
            --sampleCsv "path/to/sample.csv" \
            --cdna  "ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
            --stranded  "rf-stranded"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --sampleCsv     FILE      CSV format, header to include "sampleID,read1,read2" in case of paired data; no read2 for single-end'
  log.info ''
  log.info 'Optional arguments (must include one!):'
  log.info '    --stranded     STRING     if data is fr- (first read forward) or rf-stranded (first read reverse); '
  log.info '    --kallistoindex     STRING      suitable kallisto index'
  log.info '    OR'
  log.info '    --cdna     STRING      suitable cDNA fasta file URL'
  log.info ''
  exit 1
}

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

  input:
  file(txt) from sampleInputs

  output:
  file('inputs.txt') into sampleCsvInput

  script:
  """
  COUNT=\$(head -n1 $txt | perl -ane '@s=split(/\\,/);print scalar(@s);')
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

  input:
  set val(sampleID), file(read1), file(read2) from bbduking

  output:
  set val(sampleID), file('*pre.fastq.gz') into fastppreing
  set val(sampleID), file('*bbduk.fastq.gz') into (kallistoing, fastpposting)

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  TESTR2=\$(echo $read2 | perl -ane 'if(\$_=~m/q.gz/){print "FQ";}')
  if [[ \$TESTR2 != "FQ" ]]; then
   ln -s $read1 $sampleID".pre.fastq.gz"
   reformat.sh ${taskmem} \
      in=$sampleID".pre.fastq.gz" \
      out="stdin.fastq" \
      tossjunk=T | \
   bbduk.sh ${taskmem} \
      int=t \
      in="stdin.fastq" \
      out=$sampleID".bbduk.fastq.gz" \
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
      stats=$sampleID".bbduk.adapterstats.txt" \
      overwrite=T
  else
    ln -s $read1 $sampleID".R1.pre.fastq.gz"
    ln -s $read2 $sampleID".R2.pre.fastq.gz"

    reformat.sh ${taskmem} \
      in1=$sampleID".R1.pre.fastq.gz" \
      in2=$sampleID".R2.pre.fastq.gz" \
      out="stdout.fastq" \
      tossjunk=T | \
    bbduk.sh ${taskmem} \
      int=t \
      in="stdin.fastq" \
      out1=$sampleID".R1.bbduk.fastq.gz" \
      out2=$sampleID".R2.bbduk.fastq.gz" \
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
      stats=$sampleID".bbduk.adapterstats.txt" \
      overwrite=T
  fi
  } 2>&1 | tee $sampleID".bbduk.runstats.txt"
  """
}

/* 1.1: fastp QC
*/
fastppreing.concat(fastpposting).set { fastping }

process fastp {

  input:
  set val(sampleID), file(reads) from fastping

  output:
  file('*.json') into fastp_multiqc

  script:
  """
  PREPOST=\$(ls | grep R1 | grep fastq | perl -ane 'if(\$_=~m/bbduk/){print "post";} else {print "pre";}')
  COUNT=\$(ls *fastq.gz | wc -l)

  if [[ \$COUNT == 2 ]]; then
    fastp -w ${task.cpus} \
      -j $sampleID"."\$PREPOST".fastp.json" \
      --in1 \$(ls | grep R1 | grep fastq.gz) \
      --in2 \$(ls | grep R2 | grep fastq.gz)
  else
    fastp -w ${task.cpus} \
      -j $sampleID"."\$PREPOST".fastp.json" \
      --in1 $reads
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

  publishDir "$params.outDir/$sampleID/kallisto_fr", mode: "copy", pattern: "*"

  input:
  set sampleID, file(reads) from kallistoing
  file(kallistoindex) from kallisto_index

  output:
  file('*') into completedkallisto
  file('*.kallisto.log.txt') into kallisto_multiqc

  script:
  """
  {
  COUNT=\$(ls | grep fastq | wc -l)
  if [[ \$COUNT != 1 ]]; then
    kallisto quant \
      -t 10 \
      -b 100 \
      -i $kallistoindex \
      -o ./ --${params.stranded} $reads

  else
    kallisto quant \
      -t 10 \
      -b 100 \
      --single \
      -l 200 \
      -s 30
      -i $kallistoindex \
      -o ./ --${params.stranded} $reads
  fi
  } 2>&1 | tee > $sampleID".kallisto.log.txt"
  """
}

/* 3.0: MultiQC
 * NB bbduk_multiqc output not available for multic yet 130718
*/
fastp_multiqc.mix( kallisto_multiqc ).set { multiqc_multiqc }

process mltiqc {

  publishDir "multiqc", mode: "copy", pattern: "*"

  input:
  file ('*') from multiqc_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  multiqc -n multiqc_report.html -c ${params.multiqcConfig} -f .
  """
}
