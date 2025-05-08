process TRIM_RANDOM_BARCODE {
    tag { lane }

    input:
    tuple val(lane), path(fastq_files)

    output:
    tuple val(lane), path("output/${lane}*.fastq.gz"), emit: trimmed

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer_R1="\${barcode}${params.spacer_seq_R1}"    
    barcode_spacer_R2="\${barcode}${params.spacer_seq_R2}"
    length_barcode_spacer_R1=\${#barcode_spacer_R1}
    length_barcode_spacer_R2=\${#barcode_spacer_R2}
    min_length_read=\$(( ${params.read_length} - 4 - ${params.random_barcode_length} ))

    mkdir -p output
    
    cutadapt \
        -O \${length_barcode_spacer_R1} \
        -e ${params.spacer_error_rate} \
        -m \$min_length_read \
        -g \${barcode_spacer_R1} \
        --action=retain \
        -j ${task.cpus} \
        -o output/${lane}_R1.tmp.fastq.gz -p output/${lane}_R2.tmp.fastq.gz \
        ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
   
    cutadapt \
        -O \${length_barcode_spacer_R2} \
        -e ${params.spacer_error_rate} \
        -m \$min_length_read \
        -G \${barcode_spacer_R2} \
        --action=retain \
        -j ${task.cpus} \
        -o output/${lane}_R1.fastq.gz -p output/${lane}_R2.fastq.gz \
        output/${lane}_R1.tmp.fastq.gz output/${lane}_R2.tmp.fastq.gz

    # Clean up temporary files
    rm output/${lane}_R1.tmp.fastq.gz
    rm output/${lane}_R2.tmp.fastq.gz

    """
}