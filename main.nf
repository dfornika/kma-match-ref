#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp } from './modules/kma_match_ref.nf'
include { kma_align } from './modules/kma_match_ref.nf'
include { choose_best_ref } from './modules/kma_match_ref.nf'
include { bwa_align } from './modules/kma_match_ref.nf'
include { samtools_depth } from './modules/kma_match_ref.nf'
include { qualimap_bamqc } from './modules/kma_match_ref.nf'
include { qualimap_bamqc_genome_results_to_csv } from './modules/kma_match_ref.nf'

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

workflow {
    if (params.samplesheet_input != 'NO_FILE') {
      ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
    } else {
      ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_ref_db = Channel.fromPath( "${params.ref_db}")
    
    main:
      fastp(ch_fastq)

      kma_align(fastp.out.trimmed_reads.combine(ch_ref_db))

      choose_best_ref(kma_align.out.combine(ch_ref_db))

      bwa_align(ch_fastq.join(choose_best_ref.out))

      samtools_depth(bwa_align.out)

      qualimap_bamqc(bwa_align.out)

      qualimap_bamqc_genome_results_to_csv(qualimap_bamqc.out.genome_results)
}