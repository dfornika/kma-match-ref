# KMA Match Reference
This pipeline uses the [kma](https://bitbucket.org/genomicepidemiology/kma) aligner to align raw sequence data against
a database of similar reference genomes, and choose the one that best matches the input sequence data. Choice of best
reference sequence is made based on the Q score provided by kma for each reference.

For more details on kma parameters and outputs, see the [kma specification](https://bitbucket.org/genomicepidemiology/kma/raw/8dd45bf6e8e92eb143865433b09ef7b572f8762c/KMAspecification.pdf)

## Usage

```
nextflow run dfornika/kma-match-ref \
  --fastq_input </path/to/fastqs> \
  --ref_db </path/to/ref_db.fa> \
  --outdir </path/to/output_dir>
```

## Outputs

```
.
├── sample-01.bam
├── sample-01.bam.bai
├── sample-01_bamqc
│   ├── css
│   ├── genome_results.txt
│   ├── images_qualimapReport
│   ├── qualimapReport.html
│   └── raw_data_qualimapReport
├── sample-01_depth_by_position.tsv
├── sample-01_fastp.csv
├── sample-01_fastp.json
├── sample-01_kma.csv
├── sample-01_qualimap_bamqc_genome_results.csv
└── sample-01_ref.fa
```
