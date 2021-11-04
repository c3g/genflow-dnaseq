include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ANNOTATE {
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(vcf), path(tbi), path(dbsnp), path(dbsnp_index)

    output:
    tuple val(meta), path("${meta.id}.snpId.vcf.gz"), path("${meta.id}.snpId.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    """
    java -jar \$SNPEFF_HOME/SnpSift.jar annotate $dbsnp \\
    | bgzip -cf > ${meta.id}.snpId.vcf.gz
    tabix -pvcf ${meta.id}.snpId.vcf.gz
    java -jar /opt/snpEff/SnpSift.jar 2>&1 || true | awk '/^SnpSift version/ {print \$3}' > versions.yml
    """
}