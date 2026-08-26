process DEMULTIPLEX {
    tag { lane }

    input:
    tuple val(lane), path(barcodes), path(fastq_files)

    output:
    tuple val(lane), path("*.fastq.gz"), emit: demuxed

    script:
    barcode_R1 = barcodes[0]
    barcode_R2 = barcodes[1]
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    if [[ ${params.barcode_demux_location} == 'forward' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e ${params.barcode_demux_mismatches} \
            --no-indels \
            -g file:${barcode_R1} \
            --action=none \
            -o "${lane}#{name}_R1.fastq.gz" \
            -p "${lane}#{name}_R2.fastq.gz" \
            ${R1} ${R2}
    fi

    if [[ ${params.barcode_demux_location} == 'reverse' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e ${params.barcode_demux_mismatches} \
            --no-indels -g file:${barcode_R2} \
            --action=none \
            -o "${lane}#{name}_R2.fastq.gz" \
            -p "${lane}#{name}_R1.fastq.gz" \
            ${R2} ${R1}
    fi

    if [[ ${params.barcode_demux_location} == 'both' ]]
    then
        cutadapt \
            -j ${task.cpus} \
            -e ${params.barcode_demux_mismatches} \
            --no-indels \
            --pair-adapters \
            -g file:${barcode_R1} \
            -G file:${barcode_R2} \
            --action=none \
            -o "${lane}#{name}_R1.fastq.gz" \
            -p "${lane}#{name}_R2.fastq.gz" \
            ${R1} ${R2}
    fi
    """
}