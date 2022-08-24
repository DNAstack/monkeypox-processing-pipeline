# Monkeypox assembly, variant calling, and lineage analysis pipeline

This repository contains workflows for processing monkeypox data from Illumina.

The [FASTQ-based workflows](#fastq-based-workflows) produce variant calls, assembled genomes, and lineage assignments from raw 
sequencing reads. The FASTQ-based workflows are able to process reads originating from Illumina (paired-end) sequencing data.


## Workflows

### FASTQ-based workflows

These workflows can be used to process FASTQ files into variant calls, assembled genomes, and lineage metadata.

The inputs and outputs for the FASTQ-based workflow is outlined below for processing Illumina paired-end monkeypox sequencing 
data.


#### Workflow inputs

An input template file with some defaults pre-defined can be found 
[here](https://github.com/DNAstack/monkeypox-processing-pipeline/blob/add_monkeypox_workflow/workflows/illumina_PE/inputs.json).

| Input | Description |
|:-|:-|
| `accession` | Sample ID |
| `ref` | [The monkeypox reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_063383.1/) |
| `ref_index` | The fai index file for the reference genome. Required input for variant calling with GATK |
| `ref_dict` | The dictionary file for the reference genome. Required input for variant calling with GATK |
| `ref_amb`, `ref_ann`, `ref_bwt`, `ref_pac`, `ref_sa` | Output files after indexing with bwa  |


#### Workflow outputs

| Output | Description |
|:-|:-|
| `vcf`, `vcf_index` | Variant calls and index in VCF format |
| `sample_metadata` | Associated sample metdata (technical aspects of sequencing experiments) from [NCBI 
SRA](https://www.ncbi.nlm.nih.gov/sra) |
| `lineage_metadata` | Lineage assignment and associated metadata (tool versions, etc.) output by `Nextclade` |


## Running workflows

### Required software

- [Docker](https://docs.docker.com/get-docker/)
- [Cromwell](https://github.com/broadinstitute/cromwell/releases) & Java (8+) OR [miniwdl](https://github.com/chanzuckerberg/miniwdl/releases) & python3

### Running using Cromwell

From the root of the repository, run:

```bash
java -jar /path/to/cromwell.jar run /path/to/workflow.wdl -i /path/to/inputs.json
```

Output and execution files will be located in the `cromwell-executions` directory. When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.


### Running using miniwdl

This command assumes you have `miniwdl` available on your command line. If `miniwdl` is not available, try installing using `pip install miniwdl`.

```bash
miniwdl run /path/to/workflow.wdl -i /path/to/inputs.json
```

Output and execution files will be located in a dated directory (e.g. named `20200704_073415_main`). When the workflow finishes successfully, it will output JSON (to stdout) specifying the full path to each output file.


## Containers

Docker image definitions can be found in 
[docker](https://github.com/DNAstack/monkeypox-processing-pipeline/tree/add_monkeypox_workflow/docker). 

All containers are publicly hosted in [DNAstack's container registry](https://hub.docker.com/u/dnastack), with the exception of 
[GATK](https://hub.docker.com/r/broadinstitute/gatk/).
