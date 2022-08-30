version 1.0

import "https://raw.githubusercontent.com/DNAstack/monkeypox-processing-pipeline/main/workflows/common/common.wdl" as common

workflow monkeypox_illumina_PE {
  input {
    String accession
    String NCBI_API_KEY
    String container_registry

    File ref
    File ref_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
  }

  call common.download_fastqs {
    input:
      accession = accession,
      NCBI_API_KEY = NCBI_API_KEY,
      container_registry = container_registry
  }

  call align {
    input:
      fastq_R1 = download_fastqs.fastq_R1,
      fastq_R2 = download_fastqs.fastq_R2,
      accession = accession,
      container_registry = container_registry,
      ref = ref,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_sa = ref_sa
  }

  call mark_duplicates {
    input:
      accession = accession,
      container_registry = container_registry,
      aligned_sorted_bam = align.aligned_sorted_bam
  }

  call call_variants {
    input:
      accession = accession,
      markdup_bam = mark_duplicates.markdup_bam,
      markdup_bam_index = mark_duplicates.markdup_bam_index,
      ref = ref,
      ref_index = ref_index,
      ref_dict = ref_dict
  }

  call assemble_genome {
    input:
      accession = accession,
      container_registry = container_registry,
      aligned_sorted_bam = align.aligned_sorted_bam
  }

  call common.assign_lineage {
    input:
      accession = accession,
      container_registry = container_registry,
      assembly = assemble_genome.assembly
  }

  output { 
    Array [File] bam_files = [mark_duplicates.markdup_bam, mark_duplicates.markdup_bam_index]
    Array [File] vcf_files = [call_variants.vcf, call_variants.vcf_index] 
    Array [File] assembly_files = [assemble_genome.assembly, assemble_genome.assembly_quality]
    File sample_metadata = download_fastqs.sample_metadata
    File lineage_metadata = assign_lineage.lineage_metadata
  }

  meta {
    author: "Karen Fang"
    email: "karen@dnastack.com"
  }
}

task align {
  input {
    String accession
    String container_registry
    File fastq_R1
    File fastq_R2

    File ref
    File ref_amb # Used by bwa !UnusedDeclaration
    File ref_ann # Used by bwa !UnusedDeclaration
    File ref_bwt # Used by bwa !UnusedDeclaration
    File ref_pac # Used by bwa !UnusedDeclaration
    File ref_sa # Used by bwa !UnusedDeclaration
  }

  String platform = "ILLUMINA"
  Int disk_size_fastq = ceil((size(fastq_R1, "GB") + size(fastq_R2, "GB") + size(ref, "GB")) * 2 + 20)

  command <<<
    bwa mem \
      -R "@RG\\tID:~{accession}\\tSM:~{accession}\\tPL:~{platform}" \
      "~{ref}" \
      "~{fastq_R1}" \
      "~{fastq_R2}" \
    | samtools view \
      -b \
    > "~{accession}.aligned.bam"

    samtools sort \
      -o "~{accession}.aligned.sorted.bam" \
      "~{accession}.aligned.bam"
  >>>

  output {
    File aligned_sorted_bam = "~{accession}.aligned.sorted.bam"
  }

  runtime {
    docker: "~{container_registry}/bwa_samtools:0.0.1"
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk " + disk_size_fastq + " HDD"
    preemptible: 2
  }
}

task mark_duplicates {
  input {
    String accession
    String container_registry
    File aligned_sorted_bam
  }

  Int disk_size = ceil(size(aligned_sorted_bam, "GB") + 20)

  command <<<
    # Sort and add MS and MC tags for markdup to use
    # Only mark duplicates because GATK is duplicate-aware and will remove
    samtools collate \
      -O \
      "~{aligned_sorted_bam}" \
    | samtools fixmate \
      -m \
      -O BAM \
      - \
      - \
    | samtools sort \
    | samtools markdup \
      -O BAM \
      - \
      "~{accession}.markdup.bam"

    samtools index "~{accession}.markdup.bam"
  >>>

  output {
    File markdup_bam = "~{accession}.markdup.bam"
    File markdup_bam_index = "~{accession}.markdup.bam.bai"
  }

  runtime {
    docker: "~{container_registry}/bwa_samtools:0.0.1"
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 2
  }
}

task call_variants {
  input {
    String accession
    File markdup_bam
    File markdup_bam_index

    File ref
    File ref_index
    File ref_dict
  }

  Int disk_size = ceil(size(markdup_bam, "GB") + 20)

  command <<<
    gatk HaplotypeCaller \
      -R "~{ref}" \
      -I "~{markdup_bam}" \
      -O "~{accession}.vcf.gz" \
  >>>

  output {
    File vcf = "~{accession}.vcf.gz"
    File vcf_index = "~{accession}.vcf.gz.tbi"
  }

  runtime {
    docker: "broadinstitute/gatk:4.2.6.1" 
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 2 
    bootDiskSizeGb: 20
  }
}

task assemble_genome {
  input {
    String accession
    String container_registry
    File aligned_sorted_bam
  }

  Int mpileup_depth = 100000
  Float freq_threshold = 0.25
  Float min_coverage_depth = 10

  Int disk_size = ceil(size(aligned_sorted_bam, "GB") * 2 + 20)

  command <<<
    samtools mpileup \
      -aa \
      -A \
      -d ~{mpileup_depth} \
      -Q0 \
      ~{aligned_sorted_bam} \
    | ivar consensus \
      -t ~{freq_threshold} \
      -m ~{min_coverage_depth} \
      -n \
      -N \
      -p ~{accession}

    # Change header name to accession
    sed -i "1s;.*;>~{accession};" "~{accession}.fa"
  >>>

  output {
    File assembly = "~{accession}.fa"
    File assembly_quality = "~{accession}.qual.txt"
  }

  runtime {
    docker: "~{container_registry}/ivar:1.3.1"
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 2
  }
}
