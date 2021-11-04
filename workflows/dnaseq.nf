/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDnaseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome'] : [publish_files: false]
def bwa_mem_options = modules['bwa_mem']
def skewer_options = modules['skewer']
def sambamba_merge_options = modules['sambamba_merge']
def gatk_markduplicates_options = modules['gatk_markduplicates']
def gatk_baserecalibrator_options = modules['gatk_baserecalibrator']
def gatk_applybqsr_options = modules['gatk_applybqsr']
def gatk_haplotypecaller_options = modules['gatk_haplotypecaller']

// MODULE: Local to the pipeline
include { GET_SOFTWARE_VERSIONS                      } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { SKEWER_TRIMMING                            } from '../modules/local/skewer'                addParams( options: skewer_options )
include { BWA_MEM                                    } from '../modules/local/bwa'                   addParams( options: bwa_mem_options )
include { SAMBAMBA_MERGE                             } from '../modules/local/bwa'                   addParams( options: sambamba_merge_options )
include { MARKDUPLICATES as GATK_MARKDUPLICATES      } from '../modules/local/gatk'                  addParams( options: gatk_markduplicates_options )
include { BASERECALIBRATOR as GATK_BASERECALIBRATOR  } from '../modules/local/gatk'                  addParams( options: gatk_baserecalibrator_options )
include { APPLYBQSR as GATK_APPLYBQSR                } from '../modules/local/gatk'                   addParams( options: [publish_files : ['bam':'alignments']] )
include { HAPLOTYPECALLER  as GATK_HAPLOTYPECALLER   } from '../modules/local/gatk'                  addParams( options: gatk_haplotypecaller_options )
include { MERGE_VCFS_AND_CALL                        } from '../modules/local/gatk'
include { INDEX_VCF as INDEX_KNOWN_VCFS              } from '../modules/local/gatk'
include { INDEX_VCF as INDEX_CHROM_VCFS              } from '../modules/local/gatk'
include { INDEX_VCF as INDEX_DBSNP_VCFS              } from '../modules/local/gatk'
include { ANNOTATE as SNPEFF_ANNOTATE                } from '../modules/local/snpeff'

// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
include { INPUT_CHECK    } from '../subworkflows/local/input_check'    addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams( genome_options: publish_genome_options )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

// MODULE: Installed directly from nf-core/modules
def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DNASEQ {

    INPUT_CHECK (ch_input)
    INPUT_CHECK.out.reads | FASTQC
    INPUT_CHECK.out.reads | SKEWER_TRIMMING

    PREPARE_GENOME()

    // Map all reads to the genome
    BWA_MEM ( SKEWER_TRIMMING.out.reads, PREPARE_GENOME.out.bwa_index ).bam \
    | map { meta, bam -> [[id:meta.sample], bam] } \
    | groupTuple() \
    | SAMBAMBA_MERGE \
    | GATK_MARKDUPLICATES

    // If the user supplies vcf files of known variants,
    // use these variants to recalibrate the mapping.
    if(params.known_variants) {
        Channel.from ( params.known_variants.split(",").collect { path -> file(path) } ) \
        | map { [ [id: it.getBaseName() ], it] } \
        | INDEX_KNOWN_VCFS

        GATK_MARKDUPLICATES.out.bams \
        | combine ( INDEX_KNOWN_VCFS.out.indexed_vcfs ) \
        | combine ( PREPARE_GENOME.out.indexed_fasta ) \
        | GATK_BASERECALIBRATOR

        ChBamsForGenotyping = GATK_APPLYBQSR ( GATK_BASERECALIBRATOR.out.all ).bams
    } else {
        ChBamsForGenotyping = MARK_DUPLICATES.out.bams
    }

    ChBamsForGenotyping \
    | combine ( PREPARE_GENOME.out.fasta_records ) \
    | combine ( PREPARE_GENOME.out.indexed_fasta ) \
    | GATK_HAPLOTYPECALLER

    INDEX_CHROM_VCFS ( GATK_HAPLOTYPECALLER.out.vcfs )

    INDEX_CHROM_VCFS.out.indexed_vcfs \
    | map { meta, vcf, tbi -> [[id: meta.id], vcf, tbi]} \
    | groupTuple() \
    | set { ChSampleVCFs }

    MERGE_VCFS_AND_CALL ( ChSampleVCFs, PREPARE_GENOME.out.indexed_fasta )

    Channel.fromPath(params.dbsnp) \
    | map { path -> [[id:'dbsnp'], path]} \
    | INDEX_DBSNP_VCFS

    MERGE_VCFS_AND_CALL.out \
    | combine ( INDEX_DBSNP_VCFS.out.indexed_vcfs_nometa ) \
    | SNPEFF_ANNOTATE

    ch_software_versions = Channel.empty()
    ch_software_versions = ch_software_versions.mix(
        FASTQC.out.versions.first().ifEmpty(null),
        SKEWER_TRIMMING.out.versions.first().ifEmpty(null),
        BWA_MEM.out.versions.first().ifEmpty(null),
        GATK_BASERECALIBRATOR.out.versions.first().ifEmpty(null),
    )

    // MODULE: Pipeline reporting
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect() )

    // MODULE: MultiQC
    workflow_summary    = WorkflowDnaseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (ch_multiqc_files.collect())
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
