# Monkeypox virus analysis using Illumina PE reads

This workflow pulls paired-end Illumina Monkeypox reads from NCBI's SRA and runs assembly, variant calling, and lineage analysis.


## Tools

- [sra toolkit](https://github.com/ncbi/sra-tools) and [sra-toolkit docker](https://github.com/DNAstack/covid-processing-pipeline/tree/master/dockerfiles) 2.10.7
- [bwa mem](https://github.com/lh3/bwa) 0.7.17
- [samtools](https://github.com/samtools/samtools) 1.15
- [GATK](https://github.com/broadinstitute/gatk) 4
- [iVar](https://github.com/andersen-lab/ivar) 1.3.1
- [nextclade](https://github.com/nextstrain/nextclade) 2.4.0


## Workflow

### Workflow inputs

An input template file with some defaults pre-defined can be found [here](https://github.com/DNAstack/monkeypox-processing-pipeline/blob/main/workflows/illumina_PE/inputs.json).
   
| Input | Description |
|:-|:-|
| `NCBI_API_KEY` | An API key [from the NCBI](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) |
| `accession` | [NCBI](https://www.ncbi.nlm.nih.gov/sra/?term=%22Monkeypox+virus%22%5Borgn%3A__txid10244%5D) run accession ([SED]RRxxxx); Sample ID |
| `container_registry` | Location of container image registry. Default is set to `dnastack` |
| `ref` | [The monkeypox reference genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_063383.1/) |
| `ref_index` | The fai index file for the reference genome. Required input for variant calling with GATK |
| `ref_dict` | The dictionary file for the reference genome. Required input for variant calling with GATK |
| `ref_amb`, `ref_ann`, `ref_bwt`, `ref_pac`, `ref_sa` | Output files after [indexing with bwa](#preparing-bwa-index-files) |

#### Test inputs

Test inputs for use GCP can be found [here](gcp_test_inputs.json). The NCBI_API_KEY will need to be provided. Alternatively, reference files can be downloaded by visiting [this page](https://console.cloud.google.com/storage/browser/public_workflow_resources/monkeypox) (Google login required).


### Workflow outputs
   
| Output | Description |
|:-|:-|
| `fastq_R1`, `fastq_R2` | Paired-end raw reads in fastq files |
| `aligned_sorted_bam` | Sorted alignments in BAM format |
| `markdup_bam`, `markdup_bam_index` | Marked duplicate and sorted alignments and index in BAM format |
| `assembly`, `assembly_quality` | Assembled monkeypox genome and corresponding quality metrics |
| `vcf`, `vcf_index` | Variant calls and index in VCF format |
| `sample_metadata` | Associated sample metdata (technical aspects of sequencing experiments) from [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) |
| `lineage_metadata` | Lineage assignment and associated metadata (tool versions, etc.) output by `Nextclade` |


## Methods

### Preparing bwa index files

This command will prepare the amb, ann, bwt, pac, and sa files required by bwa.

```bash
bwa index NC_063383.1.fa
```

### Preparing reference index and dict files

These files are required by GATK HaplotypeCaller.

```bash
# Make ref index file
samtools faidx NC_063383.1.fa

# Make ref dict file
samtools dict NC_063383.1.fa -o NC_063383.1.dict
```


## Containers

Docker image definitions can be found in our [bioinformatics-public-docker-images](https://github.com/DNAstack/bioinformatics-public-docker-images) repo.

All containers are publicly hosted in [DNAstack's container registry](https://hub.docker.com/u/dnastack), with the exception of
[GATK](https://hub.docker.com/r/broadinstitute/gatk/).
