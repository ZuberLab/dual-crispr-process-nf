process TRIM_RANDOM_BARCODE {
    tag { lane }

    input:
    tuple val(lane), path(fastq_files)

    output:
    tuple val(lane), path("output/${lane}*.fastq.gz"), emit: trimmed
    tuple val(lane), path("output/untrimmed/${lane}*.fastq.gz"), emit: untrimmed

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    barcode=\$(printf "%${params.barcode_length}s" | tr ' ' "N")
    barcode_spacer_R1="\${barcode}${params.spacer_seq_R1}"
    barcode_spacer_R2="\${barcode}${params.spacer_seq_R2}"
    length_barcode_spacer_R1=\${#barcode_spacer_R1}
    length_barcode_spacer_R2=\${#barcode_spacer_R2}

    mkdir -p output
    mkdir -p output/untrimmed

    #get the read length from fastq and save it in the read_length
    read_length_R1=\$(zcat ${lane}_R1.fastq.gz | head -n 2 | tail -n 1 | wc -c)
    read_length_R1=\$(( read_length_R1 - 1 )) # subtract 1 for the newline character
    read_length_R2=\$(zcat ${lane}_R2.fastq.gz | head -n 2 | tail -n 1 | wc -c)
    read_length_R2=\$(( read_length_R2 - 1 )) # subtract 1 for the newline character

    # calculate the minimum length of the read after trimming, subtract 4 (maximum stagger length) and the random barcode
    min_length_read_R1=\$(( read_length_R1 - 4 - ${params.random_barcode_length} ))
    min_length_read_R2=\$(( read_length_R2 - 4 - ${params.random_barcode_length} ))

    cutadapt \
        -O \${length_barcode_spacer_R1} \
        -e ${params.spacer_error_rate} \
        -g \${barcode_spacer_R1} \
        --action=retain \
        --untrimmed-output output/untrimmed/${lane}_R1.pass1.fastq.gz \
        --untrimmed-paired-output output/untrimmed/${lane}_R2.pass1.fastq.gz \
        -j ${task.cpus} \
        -o output/${lane}_R1.tmp.fastq.gz -p output/${lane}_R2.tmp.fastq.gz \
        ${lane}_R1.fastq.gz ${lane}_R2.fastq.gz
        #-m \$min_length_read_R1 \

    cutadapt \
        -O \${length_barcode_spacer_R2} \
        -e ${params.spacer_error_rate} \
        -g \${barcode_spacer_R2} \
        --action=retain \
        --untrimmed-output output/untrimmed/${lane}_R2.pass2.fastq.gz \
        --untrimmed-paired-output output/untrimmed/${lane}_R1.pass2.fastq.gz \
        -j ${task.cpus} \
        -o output/${lane}_R2.fastq.gz -p output/${lane}_R1.fastq.gz \
        output/${lane}_R2.tmp.fastq.gz output/${lane}_R1.tmp.fastq.gz
        #-m \$min_length_read_R2 \

    # Combine reads rejected in either pass into a single set of untrimmed reads
    cat output/untrimmed/${lane}_R1.pass1.fastq.gz output/untrimmed/${lane}_R1.pass2.fastq.gz > output/untrimmed/${lane}_R1.fastq.gz
    cat output/untrimmed/${lane}_R2.pass1.fastq.gz output/untrimmed/${lane}_R2.pass2.fastq.gz > output/untrimmed/${lane}_R2.fastq.gz
    rm output/untrimmed/${lane}_R1.pass1.fastq.gz output/untrimmed/${lane}_R1.pass2.fastq.gz
    rm output/untrimmed/${lane}_R2.pass1.fastq.gz output/untrimmed/${lane}_R2.pass2.fastq.gz

    # Clean up temporary files
    rm output/${lane}_R1.tmp.fastq.gz
    rm output/${lane}_R2.tmp.fastq.gz

    """
}
