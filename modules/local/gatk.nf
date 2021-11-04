include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MARKDUPLICATES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bai"), emit: bams
    path "versions.yml",                                             emit: versions

    """
    gatk MarkDuplicates $options.args \\
    --REMOVE_DUPLICATES false \\
    --VALIDATION_STRINGENCY SILENT \\
    --CREATE_INDEX true \\
    --METRICS_FILE ${meta.id}.markdup_metrics.txt \\
    --INPUT $bam \\
    --OUTPUT ${meta.id}.bam
    gatk --version | awk '/Picard/ {print("    picard:", \$NF)}' > versions.yml
    """
}

process BASERECALIBRATOR {
    tag "${meta.id}"
    time '12h'

    input:
    tuple val(meta), path(bam), path(bai), val(varmeta), path(vcfs), path(tbis), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.recal_data.table"), path(fasta), path(fai), path(dict), path(bam), path(bai), emit: all
    path "versions.yml", emit: versions

    script:
    def known_sites_args = vcfs.collect{ "--known-sites $it" }.join(' ')
    """
    gatk BaseRecalibrator $options.args \\
    --reference $fasta \\
    $known_sites_args \\
    --input $bam \\
    --output ${meta.id}.recal_data.table
    gatk --version | awk '/The Genome Analysis Toolkit/ {print("    GATK:", \$NF)}' > versions.yml
    """
}

process APPLYBQSR {
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(recal_table), path(fasta), path(fai), path(dict), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bai"), emit: bams
    path "versions.yml", emit: versions

    """
    gatk ApplyBQSR $options.args \\
    --create-output-bam-index true \\
    --reference $fasta \\
    --input $bam \\
    --bqsr-recal-file $recal_table \\
    --output ${meta.id}.recal.bam
    gatk --version | awk '/The Genome Analysis Toolkit/ {print("    GATK:", \$NF)}' > versions.yml
    """
}

process INDEX_VCF {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path(vcf), path("*.tbi"), emit: indexed_vcfs
    tuple path(vcf), path("*.tbi"),            emit: indexed_vcfs_nometa
    path("versions.yml"),                      emit: versions

    """
    tabix -p vcf $vcf
    tabix --version | awk '/tabix/ {print "    tabix:", \$NF}' > versions.yml
    """
}

process HAPLOTYPECALLER {
    tag "${meta.id}"
    cpus 4

    input:
    tuple val(meta), path(bam), path(bai), val(fasta_meta), path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.${fasta_meta.seqname}.g.vcf.gz"), emit: vcfs
    path "versions.yml", emit: versions

    when:
    fasta_meta.seqlength.toInteger() >= params.min_seqlength

    script:
    """
    gatk HaplotypeCaller $options.args \\
    --emit-ref-confidence GVCF --max-reads-per-alignment-start 0 -G StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads ${task.cpus} \\
    --input $bam \\
    --output ${meta.id}.${fasta_meta.seqname}.g.vcf.gz \\
    --intervals ${fasta_meta.seqname} \\
    --reference $fasta
    gatk --version | awk '/The Genome Analysis Toolkit/ {print("    GATK:", \$NF)}' > versions.yml
    """
}

process MERGE_VCFS_AND_CALL {
    tag "$meta.id"

    input:
    tuple val(meta), path(vcfs, stageAs: "vcfs/*"), path(tbis, stageAs: "vcfs/*")
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi")

    script:
    if(vcfs instanceof Collection)
    """
    ls vcfs/*.vcf.gz > variants.list
    gatk MergeVcfs \\
      --REFERNCE $fasta \\
      --INPUT variants.list \\
      --OUTPUT ${meta.id}.g.vcf.gz
    gatk GenotypeGVCFs \\
      --reference $fasta \\
      --variant ${meta.id}.g.vcf.gz \\
      --output ${meta.id}.vcf.gz
    """
    else
    """
    gatk GenotypeGVCFs \\
      --reference $fasta \\
      --variant vcfs/*.vcf.gz \\
      --output ${meta.id}.vcf.gz
    """
}