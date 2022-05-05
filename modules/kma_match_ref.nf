process fastp {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_fastp.{json,csv}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: csv

    script:
    """
    fastp \
      -t ${task.cpus} \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}_fastp.json

    parse_fastp_json.py \
      --fastp_json ${sample_id}_fastp.json \
      --sample_id "${sample_id}" \
      > ${sample_id}_fastp.csv
    """
}

process kma_align {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_kma.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(read_1), path(read_2), path(ref_db)

    output:
    tuple val(sample_id), path("${sample_id}_kma.csv")

    script:
    """
    kma index \
      -i ${ref_db} \
      -k ${params.kmer_size} \
      -m ${params.minimizer_size} \
      -ht ${params.homology_template}  \
      -hq ${params.homology_query}

    kma \
      -t ${task.cpus} \
      -1t1 \
      -and \
      -o ${sample_id} \
      -t_db ${ref_db} \
      -ipe ${read_1} ${read_2} 

    cat ${sample_id}.res | awk '\$1 ~ /^#/ {print substr(tolower(\$0), 2)}; \$1 ~ !/^#/ {print \$0}' \
      | tr -d ' ' | tr \$'\\t' ',' > ${sample_id}_kma_unsorted.csv
    head -n 1 ${sample_id}_kma_unsorted.csv > ${sample_id}_kma.csv
    sort -t ',' -k10,10nr <(tail -qn+2 ${sample_id}_kma_unsorted.csv) >> ${sample_id}_kma.csv
    """
}

process choose_best_ref {

  tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_ref.fa", mode: 'copy'

  input:
  tuple val(sample_id), path(kma_output), path(ref_db)

  output:
  tuple val(sample_id), path("${sample_id}_ref.fa"), optional: true

  script:
  """
  choose_best_ref.py ${kma_output} \
    --db ${ref_db} \
    --min-template-identity ${params.ref_matching_min_template_identity} \
    --min-template-coverage ${params.ref_matching_min_template_coverage} \
    > best_ref_id.tsv
  num_ref_ids=\$(wc -l < best_ref_id.tsv)
  if [ \${num_ref_ids} -gt 0 ]; then
    seqtk subseq ${ref_db} best_ref_id.tsv > ${sample_id}_ref.fa
  fi
  """
}

process bwa_align {

  tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}.bam*", mode: 'copy'

  input:
  tuple val(sample_id), path(reads_1), path(reads_2), path(ref)

  output:
  tuple val(sample_id), path("${sample_id}.bam*")

  script:
  """
  bwa index ${ref}

  bwa mem -t ${task.cpus} ${ref} ${reads_1} ${reads_2} | \
    samtools sort -o ${sample_id}.bam -

  samtools index ${sample_id}.bam
  """
}

process samtools_depth {

    tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_depth_by_position.tsv", mode: 'copy'


  input:
  tuple val(sample_id), file(alignment)

  output:
  tuple val(sample_id), path("${sample_id}_depth_by_position.tsv")
    
  script:
  """
  echo "ref,position,depth" | tr ',' \$'\t' > ${sample_id}_depth_by_position.tsv
  samtools depth \
    -aa \
    -J \
    --min-MQ ${params.depth_min_mapping_quality} \
    ${alignment[0]} >> ${sample_id}_depth_by_position.tsv
  """
}

process qualimap_bamqc {

    tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_bamqc", mode: 'copy'


  input:
  tuple val(sample_id), file(alignment)

  output:
  tuple val(sample_id), path("${sample_id}_bamqc/genome_results.txt"), emit: genome_results
  tuple val(sample_id), path("${sample_id}_bamqc"), emit: bamqc_dir
    
  script:
  """
  qualimap bamqc -bam ${alignment} --outdir ${sample_id}_bamqc
  """
}

process qualimap_bamqc_genome_results_to_csv {

  tag { sample_id }

  executor 'local'

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_qualimap_bamqc_genome_results.csv", mode: 'copy'

  input:
  tuple val(sample_id), path(qualimap_bamqc_genome_results)

  output:
  tuple val(sample_id), path("${sample_id}_qualimap_bamqc_genome_results.csv")

  script:
  """
  qualimap_bamqc_genome_results_to_csv.py -s ${sample_id} ${qualimap_bamqc_genome_results} > ${sample_id}_qualimap_bamqc_genome_results.csv
  """
}

