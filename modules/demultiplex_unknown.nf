process DEMULTIPLEX_UNKNOWN {
    tag { lane }

    publishDir path: "${params.outputDir}/multiqc_unknown/${lane}",
        mode: 'copy',
        overwrite: 'true'

    input:
    tuple val(lane), val(id), path(fastq_files), path(all_3mer_fasta)

    output:
    path("*multiqc_report.html"), emit: multiqc_report

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    if [[ ${params.barcode_demux_location} == 'forward' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e 0 \
            --no-indels \
            -g file:${all_3mer_fasta} \
            --action=none \
            -o "${lane}#unknown#{name}_R2.fastq.gz" \
            -p "${lane}#unknown#{name}_R1.fastq.gz" \
            ${R1} ${R2}
    fi

    if [[ ${params.barcode_demux_location} == 'reverse' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e 0 \
            --no-indels \
            -g file:${all_3mer_fasta} \
            --action=none \
            -o "${lane}#unknown#{name}_R2.fastq.gz" \
            -p "${lane}#unknown#{name}_R1.fastq.gz" \
            ${R2} ${R1}
    fi

    if [[ ${params.barcode_demux_location} == 'both' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e 0 \
            --no-indels \
            --pair-adapters \
            -g file:${all_3mer_fasta} \
            -G file:${all_3mer_fasta} \
            --action=none \
            -o "${lane}#unknown#{name}_R2.fastq.gz" \
            -p "${lane}#unknown#{name}_R1.fastq.gz" \
            ${R1} ${R2}
    fi

    fastqc -t ${task.cpus} -q ${lane}#unknown#*.fastq.gz
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc --interactive -f -x *.run .
    """
}