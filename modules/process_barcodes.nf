process PROCESS_BARCODES {
    tag { barcodes.baseName }
    
    input:
    path(barcodes)

    output:
    tuple path("*_R1.fasta"), path("*_R2.fasta"), emit: processed_barcodes

    script:
    """
    process_barcodes.R ${barcodes}
    """
}