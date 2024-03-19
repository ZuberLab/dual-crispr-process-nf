#!/usr/bin/env nextflow

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
        --inputDir                          Input directory containing raw files.
                                            (default: '01_raw')

        --outputDir                         Output directory for processed files.
                                            (default: '02_processed')

        --library                           Path to sgRNA / shRNA library file.
                                            (default: 'library.txt')
                                            The following columns are required:
                                                - id:       unique name of sgRNA / shRNA
                                                - gene:     gene targeted by sgRNA / shRNA
                                                - sequence: nucleotide sequence of sgRNA / shRNA

         --barcodes                         Path to file containing barcodes for demultiplexing.
                                            (default: 'barcodes.fasta')
                                            The following columns are required:
                                                - lane:         name of BAM / FASTQ input file
                                                - sample_name:  name of demultiplexed sample followed by an integer number (1-4)
                                                                corresponding to the length of the stagger sequence used for this sample,
                                                                seperated by the keyword STAGGERLENGTH (e.g. STAGGERLENGTH1, STAGGERLENGTH2, ...)
                                                - barcode:      8 nucleotide long sequence of the stagger, followed by the sample barcode,
                                                                filled up to a length of 8 nucleotides with the subsequent sequence of the spacer

        --reverse_complement                Should reads be converted to reverse_complement
                                            before trimming (default: false)

        --barcode_random_length             Number of nucleotides in random barcode
                                            (default for sgRNA: 4)

        --barcode_demux_mismatches          Number of mismatches allowed during demultiplexing
                                            of barcode. (default: 1)

        --barcode_demux_location            Read location of the sample barcode. Only the specified read is used for demultiplexing.
                                            Either 'forward' or 'reverse' (default: 'forward')

        --barcode_length                    Number of nucleotides in sample barcode.
                                            (default: 4)

        --spacer_length_R1                  Number of nucleotides in spacer sequence between
                                            barcodes and sgRNA / shRNA sequence. (default: 26)

        --spacer_length_R2                  Number of nucleotides in spacer sequence between
                                            barcodes and sgRNA / shRNA sequence. (default: 33)

        --guide_length                      Number of nucleotides in guide sequence. (default: 21)

        --padding_bases_first_guide         Nucleotides used for padding if first sgRNA / shRNA are of
                                            unequal length. Must be one of G, C, T, and A.
                                            (default: GTT)

        --padding_bases_matching_guide      Nucleotides used for padding if matching sgRNA / shRNA are of
                                            unequal length. Must be one of G, C, T, and A.
                                            (default: ACC)

        --forward_read_length               Read length of the forward read, neccessary to determine post guide sequence (default: 65)

        --post_guide_sequence_nonEmpty      Sequence compositions downstream of guide sequence for non-empty library (default: 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATA')

        --post_guide_sequence_Empty         Sequence compositions downstream of guide sequence for empty library (default: 'GTTTGAGTCTTCGGTTTAAACGCGGCCGCGGATCCGAAGAC')

        --post_guide_sequence_Empty_Empty   Sequence compositions downstream of guide sequence for empty-empty library (default:'GACGATCTCAAGTCAAGC')

     Profiles:
        standard                    local execution
        singularity                 local execution with singularity
        low_ressources_cbe          SLURM execution with singularity on CBE (low ressources)
        medium_ressources_cbe       SLURM execution with singularity on CBE (medium ressources)
        high_ressources_cbe         SLURM execution with singularity on CBE (high ressources)


     Docker:
     zuberlab/crispr-nf:latest

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
log.info " barcode random (nt)              : ${params.barcode_random_length}"
log.info " barcode demultiplex (nt)         : ${params.barcode_length}"
log.info " spacer R1 (nt)                   : ${params.spacer_length_R1}"
log.info " spacer R2 (nt)                   : ${params.spacer_length_R2}"
log.info " demultiplex mismatches           : ${params.barcode_demux_mismatches}"
log.info " sample barcode location          : ${params.barcode_demux_location}"
log.info " first guide padding base         : ${params.padding_bases_first_guide}"
log.info " matching guide padding base      : ${params.padding_bases_matching_guide}"
log.info " reverse complement               : ${params.reverse_complement}"
log.info " post guide sequnce non-empty     : ${params.post_guide_sequence_nonEmpty}"
log.info " post guide sequnce empty         : ${params.post_guide_sequence_Empty}"
log.info " post guide sequnce empty-empty   : ${params.post_guide_sequence_Empty_Empty}"
log.info " forward read length              : ${params.forward_read_length}"

log.info " =============================="
log.info ""

Channel
    .fromPath( "${params.inputDir}/*.bam" )
    .map { file -> tuple( file.baseName, file ) }
    .set { rawBamFiles }

Channel
    .fromPath( "${params.inputDir}/*.{fastq,fq}.gz" )
    .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
    .set { rawFastqFiles }

Channel
    .fromPath(params.barcodes)
    .set { processBarcodeFiles }

Channel
    .fromPath(params.library)
    .into { processLibraryFile ; combineLibraryFile }

process bam_to_fastq {

    tag { lane }

    input:
    set val(lane), file(bam) from rawBamFiles

    output:
    set val(lane), file("${lane}.tar") into fastqFilesFromBam

    script:
    """
    samtools fastq -@ ${task.cpus} ${bam} -1 ${lane}_R1.fastq.gz -2 ${lane}_R2.fastq.gz
    tar -c --use-compress-program=pigz -f ${lane}.tar ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    rm ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    """
}

senseFastqFiles = Channel.create()
antisenseFastqFiles = Channel.create()

fastqFilesFromBam
    .mix(rawFastqFiles)
    .choice( senseFastqFiles, antisenseFastqFiles ) { params.reverse_complement ? 1 : 0 }

process reverse_complement_fastq {

    tag { lane }

    input:
    set val(lane), file(fastq) from antisenseFastqFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into rcfastqFiles

    script:
    """
    mv ${fastq} input.fastq.gz
    zcat input.fastq.gz | fastx_reverse_complement \
        -z \
        -o ${lane}.fastq.gz
    """
}

senseFastqFiles
    .mix(rcfastqFiles)
    .set { fastqFiles }

process trim_random_barcode {

    tag { lane }

    input:
    set val(lane), file(files) from fastqFiles

    output:
    set val(lane), file("${lane}.tar") into randomBarcodeTrimmedFiles

    script:
    position = params.barcode_random_length
    """
    tar -x --use-compress-program=pigz -f ${files}

    if [[ ${position} -gt 0 ]]
    then
      mv ${lane}_R1.fastq.gz input.fastq.gz
      cutadapt -u ${position} -j ${task.cpus} input.fastq.gz -o ${lane}_R1.fastq.gz

      mv ${lane}_R2.fastq.gz input.fastq.gz
      cutadapt -u ${position} -j ${task.cpus} input.fastq.gz -o ${lane}_R2.fastq.gz
    else
      mv ${lane}_R1.fastq.gz input.fastq.gz
      path_to_file=\$(readlink -f ${lane}_R1.fastq.gz)
      cp \${path_to_file} ${lane}_R1.fastq.gz

      mv ${lane}_R2.fastq.gz input.fastq.gz
      path_to_file=\$(readlink -f ${lane}_R2.fastq.gz)
      cp \${path_to_file} ${lane}_R2.fastq.gz
    fi

    rm input.fastq.gz
    rm ${lane}.tar
    tar -c --use-compress-program=pigz -f ${lane}.tar ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    """
}

process process_barcodes {

    tag { barcodes.baseName }

    input:
    file(barcodes) from processBarcodeFiles

    output:
    file("*.fasta") into demuxBarcodeFiles

    script:
    """
    process_barcodes.R ${barcodes}
    """
}

demuxBarcodeFiles
    .flatten()
    .map { file -> [file.baseName, file] }
    .concat(randomBarcodeTrimmedFiles)
    .groupTuple()
    .set { demuxFiles }

process demultiplex {

    tag { lane }

    input:
    set val(lane), file(files) from demuxFiles

    output:
    set val(lane), file('*.demux.tar') into splitFiles, splitFilesFastQC

    script:
    """
    tar -x --use-compress-program=pigz -f ${files[1]}

    if [[ ${params.barcode_demux_location} == 'forward' ]]
    then
        cutadapt -j ${task.cpus} -e ${params.barcode_demux_mismatches} --no-indels -g file:${files[0]} --action=none -o "{name}_R1.demux.fastq.gz" -p "{name}_R2.demux.fastq.gz" ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
    elif [[ ${params.barcode_demux_location} == 'reverse' ]]
    then
        cutadapt -j ${task.cpus} -e ${params.barcode_demux_mismatches} --no-indels -g file:${files[0]} --action=none -o "{name}_R2.demux.fastq.gz" -p "{name}_R1.demux.fastq.gz" ${lane}_R2.fastq.gz ${lane}_R1.fastq.gz
    else
        echo ${params.barcode_demux_location} not a valid option for the parameter barcode_demux_location.
        exit
    fi

    mv unknown_R1.demux.fastq.gz ${lane}_unknown_R1.demux.fastq.gz
    mv unknown_R2.demux.fastq.gz ${lane}_unknown_R2.demux.fastq.gz

    for file in *_R1.demux.fastq.gz;
    do
      filename="\${file%_R1.demux.fastq.gz}"
      tar -c --use-compress-program=pigz -f \${filename}.demux.tar \${filename}_R1.demux.fastq.gz \${filename}_R2.demux.fastq.gz
    done
    """
}

def ungroupTuple = {
    def result = []
    def name = it[0]
    it[1].each { result << [name, it] }
    return result
}

splitFilesFastQC
    .flatMap { it -> ungroupTuple(it) }
    .map { lane, file -> tuple(lane, file.name.replaceAll(/\.demux\.tar/, ''), file) }
    .set { fastqcSplitFiles }

splitFiles
    .flatMap { it -> ungroupTuple(it) }
    .filter { it[1].baseName =~ /^(?!.*unknown).*$/ }
    .map { lane, file -> tuple(lane, file.name.replaceAll(/\.demux\.tar/, ''), file) }
    .set { flattenedSplitFiles }

process trim_barcode_and_spacer {

    tag { id }

    publishDir path: "${params.outputDir}/fastq/${lane}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(lane), val(id), file(files) from flattenedSplitFiles

    output:
    set val(lane), val(id), file("${id}.trimmed.fastq.gz") into spacerTrimmedFiles

    script:
    barcode_spacer_length_R1 = params.spacer_length_R1 + params.barcode_length
    barcode_spacer_length_R2 = params.spacer_length_R2 + params.barcode_length

    """
    tar -x --use-compress-program=pigz -f ${files}

    str=${id}
    stagger_length=\${str: -1}
    remove_beginning_R1=\$(expr \${stagger_length} + ${barcode_spacer_length_R1})
    remove_beginning_R2=\$(expr \${stagger_length} + ${barcode_spacer_length_R2})

    cutadapt ${id}_R1.demux.fastq.gz -j ${task.cpus} -u \${remove_beginning_R1} -o ${id}_R1_trimmed_beginning.fastq.gz
    cutadapt ${id}_R2.demux.fastq.gz -j ${task.cpus} -u \${remove_beginning_R2} -o ${id}_R2_trimmed_beginning.fastq.gz

    cutadapt ${id}_R1_trimmed_beginning.fastq.gz -j ${task.cpus} -l ${params.guide_length} -o ${id}_R1.trimmed.fastq.gz
    rm ${id}_R1_trimmed_beginning.fastq.gz
    
    cutadapt ${id}_R2_trimmed_beginning.fastq.gz -j ${task.cpus} -l ${params.guide_length} -o ${id}_R2.trimmed.fastq.gz
    rm ${id}_R2_trimmed_beginning.fastq.gz

    java -ea -Xmx30g -cp /opt/conda/envs/dual-crispr-nf/opt/bbmap-39.01-1/current/ jgi.FuseSequence in1=${id}_R1.trimmed.fastq.gz in2=${id}_R2.trimmed.fastq.gz pad=0 out=${id}.trimmed.fastq.gz fusepairs=t

    """
}

process process_library {

    tag { library.baseName }

    input:
    file(library) from processLibraryFile

    output:
    file("${library.baseName}.saf") into librarySafFile
    file("${library.baseName}.fasta") into libraryFastaFile

    script:
    """
    process_library.R ${library} ${params.padding_bases_first_guide} ${params.padding_bases_matching_guide}
    """
}

process bowtie_index {

    tag { library_fasta.baseName }

    input:
    file(library_fasta) from libraryFastaFile

    output:
    file("bt2") into bt2Index

    script:
    """
    mkdir -p bt2
    bowtie2-build ${library_fasta} bt2/index
    """
}

process align {

    tag { id }

    input:
    set val(lane), val(id), file(fastq) from spacerTrimmedFiles
    each file(index) from bt2Index

    output:
    set val(lane), val(id), file("${id}.sam") into alignedFiles
    file "${id}.log" into alignResults
    file "${id}.sam.stats" into alignStats
    file "${id}.sam.flagstat" into alignFlagstats

    script:
    """
    bowtie2 \
        --threads \$((${task.cpus})) \
        -x ${index}/index \
        -L 31 \
        -N 0 \
        --seed 42 \
        <(zcat ${fastq}) 2> ${id}.log > ${id}.sam

    samtools stats ${id}.sam > ${id}.sam.stats
    samtools flagstats ${id}.sam > ${id}.sam.flagstat

    """
}

alignedFiles
    .map { lane, id, file -> tuple(lane, file) }
    .groupTuple()
    .set { groupedAlignedFiles }

process count {

    tag { lane }

    publishDir path: "${params.outputDir}/counts/${saf.baseName}/${lane}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(lane), file(sams) from groupedAlignedFiles
    each file(saf) from librarySafFile

    output:
    file("${lane}.txt") into countedFiles
    file("${lane}.txt.summary") into featureCountsResults

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${saf} \
        -F SAF \
        -o ${lane}.txt \
        ${sams}
    """
}

process combine_counts {

    tag { library.baseName }

    publishDir path: "${params.outputDir}/counts/${library.baseName}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file(counts) from countedFiles.collect()
    file(library) from combineLibraryFile

    output:
    file("counts_mageck.txt") into combinedMageckFile

    script:
    """
    combine_counts.R ${library} ${counts} > counts_mageck.txt
    """
}


fastqcSplitFiles
    .map { lane, id, file -> tuple(lane, file) }
    .groupTuple()
    .set { fastqcSplitFiles }

process fastqc {

    tag { lane }

    input:
    set val(lane), file(files) from fastqcSplitFiles

    output:
    file "*_fastqc.{zip,html}" into fastqcResults

    script:
    """
    for file in ${files};
    do
      tar -x --use-compress-program=pigz -f \${file}
    done

    fastqc -t ${task.cpus} -q *.fastq.gz
    """
}

process multiqc {

    tag { 'all' }

    publishDir path: "${params.outputDir}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file (fastqc_demux: 'fastqc_demux/*') from fastqcResults.collect()
    file (align: 'align/*') from alignResults.collect()
    file (alignStats: 'align/*') from alignStats.collect()
    file (alignFlagstats: 'align/*') from alignFlagstats.collect()
    file (featurecounts: 'featureCounts/*') from featureCountsResults.collect()

    output:
    file "*multiqc_report.html" into multiqc_report

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc --interactive -f -x *.run .
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
