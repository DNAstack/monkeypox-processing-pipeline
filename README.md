# Monkeypox assembly, variant calling, and lineage analysis pipeline

This repository contains workflows for processing [NCBI SRA monkeypox 
data](https://www.ncbi.nlm.nih.gov/sra/?term=(monkeypox)+AND+%22Monkeypox+virus%22%5Borgn%3A__txid10244%5D) from different sequencing platforms including Illumina. The 
2022 monkeypox outbreak was declared as a Public Health Emergency of International Concern by the World Health Organization (WHO). The Monkeypox virus can cause a rash 
that can present as pimples or blisters. Other symptoms can include fever, chills, headache, and exhaustion. It stems from the same family of viruses that cause smallpox; 
however, monkeypox is milder and rarely fatal, lasting for 2-4 weeks. There are two distinct genetic clades of the monkeypox virus: the Central African (Congo Basin) 
clade and the West African clade, with the latter being implicated in the current outbreak. The most updated information states that monkeypox can be spread through 
direct contact, touching objects, and respiratory secretions. 

## Workflows

These workflows can be used to process FASTQ files into variant calls, assembled genomes, and lineage assignments from raw sequencing reads.

The inputs and outputs for the FASTQ-based workflow is outlined in detail in the repository for that workflow, linked below [Illumina](#illumina).


### Illumina

This workflow is used for processing Illumina (paired-end) monkeypox sequencing data.

- [Workflow inputs](https://github.com/DNAstack/monkeypox-processing-pipeline/tree/add_monkeypox_workflow/workflows/illumina_PE#workflow-inputs)

- [Workflow outputs](https://github.com/DNAstack/monkeypox-processing-pipeline/tree/add_monkeypox_workflow/workflows/illumina_PE#workflow-outputs)


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

