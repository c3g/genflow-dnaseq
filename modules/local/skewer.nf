process SKEWER_TRIMMING {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("trim-*.fastq.gz"), emit: reads
    path  "version.yml",                      emit: versions

    """
    cat > adapters.tsv << END
    >Adapter1
    ${meta["adapter1"]}
    >Adapter2
    ${meta["adapter2"]}
    END
    skewer -x adapters.tsv --threads ${task.cpus} --min 25 --end-quality 25 --format sanger --compress --output trim ${reads.join(' ')}
    skewer --version | awk '/version/ {print("    skewer:", \$NF)}' > version.yml
    """
}