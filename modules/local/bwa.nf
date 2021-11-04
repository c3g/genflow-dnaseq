import sun.nio.fs.UnixPath
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNTAR_BWA_INDEX {
    tag "$archive"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path archive

    output:
    path "bwa_index",    emit: index
    path "versions.yml", emit: versions

    """
    mkdir bwa_index
    tar $options.args \\
        --extract \\
        --gzip \\
        --directory bwa_index \\
        --file $archive

    cat <<-END_VERSIONS > versions.yml
    UNTAR_BWA_INDEX:
        tar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml",            emit: versions

    script:
    def rg_id = "ID:" + meta.id
    def rg_sm = "SM:" + meta.sample
    def rg_lb = "LB:" + (meta.run ?: meta.id)
    def rg_pu = (meta.sample && meta.run && meta.lane) ? "PU:" + [meta.sample, meta.run, meta.lane].join(".") : null
    def rg_cn = params.seq_center ? "CN:" + params.seq_center : null
    def rg_pl = params.seq_platform ? "PL:" + params.seq_platform : null
    def rg = ["@RG", rg_id, rg_sm, rg_lb, rg_pu, rg_cn, rg_pl].findAll().join("\\t")
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bwa_threads = Math.max(1, task.cpus - 2)
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem $options.args \\
        -t ${bwa_threads} \\
        -R '$rg' \\
        \$INDEX \\
        $reads \\
    | samtools view --bam \\
    | sambamba sort \\
        /dev/stdin \\
        -m ${task.memory.toGiga()}G \\
        --out ${meta.id}.sorted.bam
    sambamba index ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    BWA_MEM:
        bwa: \$(bwa 2>&1 | grep Version | sed 's/Version: //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        sambamba: \$(sambamba --version 2>&1 | grep sambamba | head -n1 | sed 's/^sambamba //g')
    END_VERSIONS
    """
}

process SAMBAMBA_MERGE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.bam")

    script:
    if(bams instanceof Collection)
    "sambamba merge --nthreads ${task.cpus} ${meta.id}.raw.bam $bams"
    else
    "ln -s $bams ${meta.id}.raw.bam"
}

// process BWA_STATS
//     tag "$meta.id"
//     label 'process_low'
//     publishDir "${params.outdir}",
//         mode: params.publish_dir_mode,
//         saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//     input:
//     tuple val(meta), path(bam)

//     output:
//     tuple val(meta), path("*.stats.txt"), emit: stats
//     path  "versions.yml"                , emit: versions

//     """
//     samtools stats $bam > ${meta.id}.mapping.stats.txt

//     cat <<-END_VERSIONS > versions.yml
//     samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
//     END_VERSIONS
//     """
// }