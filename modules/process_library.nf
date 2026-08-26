process PROCESS_LIBRARY {
    tag { library.baseName }

    input:
    path(library)

    output:
    path("${library.baseName}.saf"), emit : saf
    path("${library.baseName}.fasta"), emit : fasta
    
    script:
    """
    process_library.R ${library} ${params.padding_bases_first_guide} ${params.padding_bases_matching_guide}
    """
}