// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DecomposeAndNormalise {
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf)
}