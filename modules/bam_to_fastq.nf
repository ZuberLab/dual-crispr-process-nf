process BAM_TO_FASTQ {
    tag { lane }

    input:
    tuple val(lane), path(bam)

    output:
    tuple val(lane), path("*.fastq.gz"), emit: fastq

    script:
    """
    samtools fastq -@ ${task.cpus} ${bam} -1 ${lane}_R1.fastq.gz -2 ${lane}_R2.fastq.gz
    """
}