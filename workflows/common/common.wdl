version 1.0

task download_fastqs {
    input {
        String accession
    }

    String NCBI_API_KEY = "bd01816fd91fb4c867a002f5e3e613f0fd08"

    command <<<
        NCBI_API_KEY=~{NCBI_API_KEY}
        export NCBI_API_KEY

        use_fastq_dump=false

        retries=5
        status=1

        fastq_dump=""
        reads_invalid=""

        while [ $retries -gt 0 ] && [ $status -ne 0 ] || [ -n "$reads_invalid" ]; do
            echo "Trying download with [$retries] retries (fasterq-dump)"

            fasterq-dump \
                --mem 3GB \
                --split-3 \
                ~{accession} \
                2> dump.stderr.txt
            status=$?
            reads_invalid=$(grep "reads invalid" dump.stderr.txt)
            fastq_dump=$(grep "fastq-dump" dump.stderr.txt)
            if [ -n "$fastq_dump" ]; then
                cat dump.stderr.txt
                use_fastq_dump=true
                break
            fi
            ((retries--))
        done

        # Try using fastq-dump rather than fasterq-dump
        if [ "$use_fastq_dump" = "true" ] || [ $status -ne 0 ] || [ -n "$reads_invalid" ]; then
            retries=5
            status=1
            reads_invalid=""
            while [ $retries -gt 0 ] && [ $status -ne 0 ] || [ -n "$reads_invalid" ]; do
                echo "Trying download with [$retries] retries (fastq-dump)"

                fastq-dump \
                    --split-e \
                    ~{accession} \
                    2> dump.stderr.txt \
                status=$?
                reads_invalid=$(grep "reads invalid" dump.stderr.txt)
                ((retries--))
            done
        fi

        gzip "~{accession}_1.fastq"
        gzip "~{accession}_2.fastq"

        get_sample_metadata.sh \
            -a ~{accession}
    >>>

    output {
        File fastq_R1 = "~{accession}_1.fastq.gz"
        File fastq_R2 = "~{accession}_2.fastq.gz"
        File sample_metadata = "~{accession}.meta.csv"
    }

    runtime {
        docker: "dnastack/sra-toolkit:2.10.7"
        cpu: 2
        memory: "7.5 GB"
        disks: "local-disk 50 HDD"
    }
}

task assign_lineage {
    input {
        String accession
        File assembly
    }

    Int disk_size = ceil(size(assembly, "GB")) + 20

    command <<<
        nextclade run \
            --input-dataset "$MONKEYPOX_DATASET" \
            --output-tsv ~{accession}.nextclade_assignment.tsv \
            ~{assembly}

        nextclade_dataset=$(< "$MONKEYPOX_DATASET/tag.json"  jq -r '.name')
        nextclade_dataset_tag=$(< "$MONKEYPOX_DATASET/tag.json"  jq -r '.tag')
        nextclade_schema_version=$(< "$MONKEYPOX_DATASET/virus_properties.json" jq -r '.schemaVersion')

        paste \
            ~{accession}.nextclade_assignment.tsv \
            <(echo -e "nextclade_dataset\\tnextclade_dataset_tag\\tnextclade_schema_version\\n$nextclade_dataset\\t$nextclade_dataset_tag\\t$nextclade_schema_version") \
        > ~{accession}.nextclade.tsv

        lineage_column=$(head -1 ~{accession}.nextclade.tsv | tr '\t' '\n' | grep -n "lineage" | cut -d ":" -f 1)
        tail -1 ~{accession}.nextclade.tsv | cut -f "$lineage_column" > ~{accession}.lineage.txt
    >>>

    output {
        File lineage_metadata = "~{accession}.nextclade.tsv"
        String lineage = read_string("~{accession}.lineage.txt")
    }

    runtime {
        docker: "dnastack/nextclade:2.4.0"
        cpu: 2
        memory: "7.5 GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
