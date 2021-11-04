include { GUNZIP as GUNZIP_FASTA }   from '../../modules/nf-core/modules/gunzip/main'
include { UNTAR_BWA_INDEX } from '../../modules/local/bwa'
include { BWA_INDEX } from '../../modules/nf-core/modules/bwa/index/main'

process BUNDLE_BWA_INDEX {
    input:
    path(indices)

    output:
    path("bwa_index"), emit: 'index'

    """
    mkdir bwa_index
    cp $indices bwa_index/
    """
}

process INDEX_FASTA {
    input:
    path(fasta)

    output:
    tuple path(fasta), path("*.fai"), path("*.dict")

    """
    samtools faidx $fasta
    gatk CreateSequenceDictionary --REFERENCE $fasta
    """
}

workflow PREPARE_GENOME {
    ch_fasta   = Channel.empty()
    ch_versions = Channel.empty()

    ch_fasta = params.fasta.endsWith('.gz') ? GUNZIP_FASTA ( file(params.fasta) ).gunzip : file(params.fasta)

    INDEX_FASTA ( ch_fasta )

    if (!params.bwa) {
        ch_bwa_index = BWA_INDEX ( ch_fasta ).index
    } else if ( params.bwa.endsWith('.tar.gz') || params.bwa.endsWith('.tgz') ) {
        ch_bwa_index = UNTAR_BWA_INDEX ( file(params.bwa) ).index
    } else {
        ch_bwa_index = BUNDLE_BWA_INDEX ( Channel.fromPath(params.bwa).toList() ).index
    }

    ch_versions  = ch_versions.mix ( BWA_INDEX.out.versions.ifEmpty(null) )

    INDEX_FASTA.out \
    | map { fasta, fai, dict -> fai } \
    | splitCsv ( sep: '\t', header: ["seqname", "seqlength", "byteindex", "width", "bytewidth"] ) \
    | set { ch_fasta_records }

    emit:
    fasta         = ch_fasta
    indexed_fasta = INDEX_FASTA.out
    bwa_index     = ch_bwa_index
    versions      = ch_versions.ifEmpty(null)
    fasta_records = ch_fasta_records
}