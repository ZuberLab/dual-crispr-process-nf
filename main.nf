#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage() {
    log.info"""
    ================================================================
     dual-crispr-process-nf
    ================================================================
     DESCRIPTION

     Process CRISPR and shRNA functional genetic screening data.

     Usage:
     nextflow run zuberlab/crispr-process-nf

     Options:
        --inputDir                          Input directory containing raw files. Either BAM or FASTQ files (R1 & R2).
                                            The FASTQ files must be named <lane>_R1.fastq.gz and <lane>_R2.fastq.gz.
                                            The BAM files must be named <lane>.bam.
                                            (default: '01_raw')

        --outputDir                         Output directory for processed files.
                                            (default: '02_processed')

        --library                           Path to sgRNA library file.
                                            (default: 'library.txt')
                                            The following columns are required:
                                                - id:       unique name of sgRNA
                                                - gene:     gene targeted by sgRNA
                                                - sequence: nucleotide sequence of sgRNA

         --barcodes                         Path to file containing barcodes for demultiplexing.
                                            (default: 'barcodes.fasta')
                                            The following columns are required:
                                                - lane:         name of BAM / FASTQ input file
                                                - sample_name:  name of demultiplexed sample 
                                                - barcode_R1:   nucleotide sequence of the sample barcode on R1
                                                - barcode_R2:   nucleotide sequence of the sample barcode on R2


        --barcode_demux_mismatches          Number of mismatches allowed during demultiplexing
                                            of barcode. (default: 0)

        --demux_unknown                     Perform demultiplexing of unknown reads with all possible 
                                            barcode combinations. (default: false)

        --barcode_demux_location            Read location of the sample barcode. Only the specified read is used for demultiplexing.
                                            Either 'both', 'forward' or 'reverse' (default: 'both')

        --barcode_length                    Number of nucleotides in sample barcode.
                                            (default: 4)

        --spacer_seq_R1                     Nucleotide sequence of spacer on R1 between
                                            barcodes and sgRNA sequence. 
                                            (default: ATATCCCTTGGAGAAAAGCCTTGTTT)

        --spacer_seq_R2                     Nucleotide sequence of spacer on R2 between
                                            barcodes and sgRNA sequence. 
                                            (default: CTTGCTATGCACTCTTGTGCTTAGCTCTGAAAC)

        --spacer_error_rate                 Error rate for spacer sequence. (default: 0.4)

        --guide_length                      Number of nucleotides in guide sequence. (default: 21)

        --padding_bases_first_guide         Nucleotides used for padding if first sgRNA are of
                                            unequal length. Must be one of G, C, T, and A.
                                            (default: GTT)

        --padding_bases_matching_guide      Nucleotides used for padding if matching sgRNA are of
                                            unequal length. Must be one of G, C, T, and A.
                                            (default: ACC)

        --read_length                       Length of reads. (default: 150)

        --random_barcode_length             Length of random barcode. (default: 4)

     Profiles:
        standard                    local execution
        apptainer                   local execution with apptainer
        cbe                         SLURM execution with apptainer on CBE cluster


     Docker:
     zuberlab/dual-crispr-nf:1.0

     Author:
     Florian Andersch (florian.andersch@imp.ac.at)
     
     
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

log.info ""
log.info " parameters "
log.info " =================================="
log.info " input directory                  : ${params.inputDir}"
log.info " output directory                 : ${params.outputDir}"
log.info " library file                     : ${params.library}"
log.info " barcode file                     : ${params.barcodes}"
log.info " barcode length                   : ${params.barcode_length}"
log.info " spacer seq R1 (nt)               : ${params.spacer_seq_R1}"
log.info " spacer seq R2 (nt)               : ${params.spacer_seq_R2}"
log.info " spacer error rate                : ${params.spacer_error_rate}"
log.info " demultiplex mismatches           : ${params.barcode_demux_mismatches}"
log.info " demultiplex unknown              : ${params.demux_unknown}"
log.info " sample barcode location          : ${params.barcode_demux_location}"
log.info " first guide padding base         : ${params.padding_bases_first_guide}"
log.info " matching guide padding base      : ${params.padding_bases_matching_guide}"
log.info " read length                      : ${params.read_length}"
log.info " random barcode length            : ${params.random_barcode_length}"
log.info " =================================="
log.info ""

// Import modules
include { BAM_TO_FASTQ } from './modules/bam_to_fastq'
include { TRIM_RANDOM_BARCODE } from './modules/trim_random_barcode'
include { PROCESS_BARCODES } from './modules/process_barcodes'
include { DEMULTIPLEX } from './modules/demultiplex'
include { DEMULTIPLEX_UNKNOWN } from './modules/demultiplex_unknown'
include { TRIM_BARCODE_AND_SPACER } from './modules/trim_barcode_and_spacer'
include { PROCESS_LIBRARY } from './modules/process_library'
include { BOWTIE_INDEX } from './modules/bowtie_index'
include { ALIGN } from './modules/align'
include { COUNT } from './modules/count'
include { COMBINE_COUNTS } from './modules/combine_counts'
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'

// Define input channels
ch_input_bam = Channel.fromPath("${params.inputDir}/*.bam")
    .map { file -> tuple(file.baseName, file) }

ch_input_fastq = Channel.fromPath("${params.inputDir}/*.fastq.gz")
    .map {
        file ->
            def lane = file.name.toString().replaceAll(/_R[12]\.fastq\.gz$/, '')
            tuple(lane, file)
        }
    .groupTuple()

ch_barcodes = Channel.fromPath(params.barcodes)
ch_library = Channel.fromPath(params.library)
all_3mer_fasta = Channel.fromPath("$projectDir/bin/all_3mer.fasta")

// Main workflow
workflow {
    // BAM to FASTQ conversion
    ch_fastq_from_bam = BAM_TO_FASTQ(ch_input_bam)

    // Combine BAM-derived and direct FASTQ inputs
    ch_all_fastq = ch_fastq_from_bam.mix(ch_input_fastq)

    // Trim random barcode
    ch_trimmed_random = TRIM_RANDOM_BARCODE(ch_all_fastq)

    // Process barcodes
    ch_processed_barcodes = PROCESS_BARCODES(ch_barcodes)

    // Combine processed barcodes with trimmed random barcodes
    ch_trimmed_random_barcodes = ch_processed_barcodes
        .flatten()
        .map { barcode -> 
            def lane = barcode.name.toString().replaceAll(/_R[12]\.fasta$/, '')
            [lane, barcode]
        }
        .groupTuple()
        .join(ch_trimmed_random)
    
    // Demultiplex
    ch_demuxed = DEMULTIPLEX(ch_trimmed_random_barcodes)

    // Flatten demultiplexed files
    ch_demuxed_flattened = ch_demuxed
        .flatMap { lane, files ->
                files.collect { file ->
                    def id = file.name.toString().replaceAll(/_R[12]\.fastq\.gz$/, '')
                    [lane, id, file]
                }
        }.groupTuple(by: [0,1])

    // Split the channel based on whether id ends with "unknown"
    ch_demuxed_flattened
        .branch {
            unknown: it[1].toString().endsWith("#unknown")
            known: true
        }
        .set { ch_demuxed_flattened_split }

    ch_unknown_with_metadata = ch_demuxed_flattened_split.unknown
        .combine(all_3mer_fasta)
        .map { lane, id, files, all_3mer_fasta -> [lane, id, files, all_3mer_fasta] }

    // Demultiplex unknown samples
    if(params.demux_unknown) {
        // Demultiplex unknown samples
        ch_demuxed_unknown = DEMULTIPLEX_UNKNOWN(ch_unknown_with_metadata)
    }

    // Trim barcode and spacer
    ch_trimmed_spacer = TRIM_BARCODE_AND_SPACER(ch_demuxed_flattened_split.known)

    // Process library
    ch_library_out = PROCESS_LIBRARY(ch_library)

    // Bowtie index
    ch_bt2_index = BOWTIE_INDEX(ch_library_out.fasta)

    // Combine trimmed spacer files and bowtie index
    ch_trimmed_spacer_combined = ch_trimmed_spacer
        .combine(ch_bt2_index)
    
    // Align
    ch_aligned = ALIGN(ch_trimmed_spacer_combined)

    // Group aligned files by lane
    ch_grouped_aligned =  ch_aligned.alignedFiles
        .map { lane, id, file -> tuple(lane, file) }
        .groupTuple()
        .combine(ch_library_out.saf)

    // Count
    ch_counted = COUNT(ch_grouped_aligned)

    // Combine counts
    COMBINE_COUNTS(ch_counted.countedFiles.collect(), ch_library)

    // collect all fastq files
    ch_fastq_files = ch_all_fastq
        .mix(ch_demuxed_flattened.map { lane, baseName, file -> tuple(baseName, file) })

    // FastQC
    ch_fastqc = FASTQC(ch_fastq_files)

    // Combine FASTQC files, alignment results, and featureCounts results
    ch_multiqqc_files = ch_fastqc.collect()
        .mix(ch_aligned.alignResults.collect())
        .mix(ch_aligned.alignStats.collect())
        .mix(ch_aligned.alignFlagstats.collect())
        .mix(ch_counted.featureCountsResults.collect())
        .collect()

    // MultiQC
    MULTIQC(ch_multiqqc_files)
}

// On completion
workflow.onComplete {
    println ( workflow.success ? "COMPLETED!" : "FAILED" )
}