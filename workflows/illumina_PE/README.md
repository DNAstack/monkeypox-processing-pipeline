# Monkeypox virus analysis using Illumina PE reads

This workflow pulls paired-end Illumina Monkeypox reads from NCBI's SRA and runs assembly, variant calling, and lineage analysis.


## Tools

- [sra toolkit](https://github.com/ncbi/sra-tools)
- [bwa mem](https://github.com/lh3/bwa) 0.7.17
- [samtools](https://github.com/samtools/samtools) 1.15
- [GATK](https://github.com/broadinstitute/gatk) 4
- [iVar](https://github.com/andersen-lab/ivar) 1.3.1
- [nextclade](https://github.com/nextstrain/nextclade) 2.4.0


## Workflow inputs

- `accession`: [NCBI](https://www.ncbi.nlm.nih.gov/sra/?term=%22Monkeypox+virus%22%5Borgn%3A__txid10244%5D) run accession ([SED]RRxxxx)
- `ref`, `ref_amb`, `ref_ann`, `ref_bwt`, `ref_pac`, `ref_sa`: Monkeypox reference genome ([NC_063383](https://www.ncbi.nlm.nih.gov/nuccore/NC_063383)) and [associated bwa index files](#preparing-bwa-index-files)
- `ref_index`, `ref_dict`: Monkeypox reference genome [index and dict files](#preparing-reference-index-and-dict-files); required by GATK HaplotypeCaller


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


## Workflow outputs

_For each workflow output:_

_- `output_name`: description_

_Alternatively, a markdown table can be used._
