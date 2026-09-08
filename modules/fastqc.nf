process FASTQC {
    tag { id }

    input:
    tuple val(id), path(fastq_files)

    output:
    path "*_fastqc.{zip,html}", emit: fastqcResults

    script:
    """
    for f in ${fastq_files}; do
        ln -s "\$f" "${id}__\$(basename \$f)"
    done
    fastqc -t ${task.cpus} -q ${id}__*.fastq.gz
    """
}