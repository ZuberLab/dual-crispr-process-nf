process ALIGN {
    tag { "${id}" }

    input:
    tuple val(lane), val(id), path(fastq_files), path(index)

    output:
    tuple val(lane), val(id), path("${id}.sam"), emit: alignedFiles
    path("${id}.log"), emit: alignResults
    path("${id}.sam.stats"), emit: alignStats
    path("${id}.sam.flagstat"), emit: alignFlagstats

    script:
    R1 = fastq_files[0]
    R2 = fastq_files[1]
    """
    bowtie2 \
        --seed 42 \
        --threads \$((${task.cpus})) \
        -x ${index}/index \
        --ignore-quals \
        -L 12 \
        -N 1 \
        --score-min L,-42,0 \
        --no-overlap \
        --no-contain \
        -I 43 -X 43 \
        --fr \
        -1 ${R1} \
        -2 ${R2} 2> ${id}.log > ${id}.sam


    # Generate statistics for the aligned SAM file
    samtools stats ${id}.sam > ${id}.sam.stats
    samtools flagstats ${id}.sam > ${id}.sam.flagstat

    """
}