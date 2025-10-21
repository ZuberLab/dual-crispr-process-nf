#!/usr/bin/env Rscript

################################################################################
# process library text file
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp, extended and curated by Florian Andersch
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2023/03/07
################################################################################

# library file should be tab-separated text file with three columns:
# 1) id      : unique id of sgRNA / shRNA
# 2) group   : group targeted by sgRNA / shRNA (e.g gene id or domain)
# 3) sequence: sequence of sgRNA / shRNA as it appears in the sequencing reads

### command line parameters
args         <- commandArgs(trailingOnly = TRUE)
input_file   <- args[1]
padding_bases_first_guide <- args[2]
padding_bases_matching_guide <- args[3]

library_name <- stringr::str_replace(basename(input_file), ".txt", "")

### functions
`%>%` <- dplyr::`%>%`

### import
raw <- readr::read_tsv(input_file)

###########
# library #
###########

sequence_first_guide <- paste0(raw$sequence, padding_bases_first_guide) %>%
  substr(start=0, stop=21)

sequence_matching_guide <- paste0(padding_bases_matching_guide, raw$sequence_matching) %>%
  stringr::str_sub(start= -21)

sequence <- paste0(sequence_first_guide, "N", sequence_matching_guide)
id <- raw$id %>% make.unique()

### check for id and sequence duplication
# stopifnot(!any(duplicated(raw$id)))
# stopifnot(!any(duplicated(sequence)))

### generate fasta file for bowtie2 index
seq_length <- max(nchar(sequence))

sequence %>%
  # toupper %>%
  # stringr::str_pad(pad = padding_base, width = seq_length, side = "left") %>%
  purrr::set_names(id) %>%
  Biostrings::DNAStringSet() %>%
  Biostrings::writeXStringSet(paste0(library_name, ".fasta"), format = "fasta")

### generate SAF annotation file for featureCount (subread package)
tibble::tibble(GeneID = id,
               Chr    = id,
               Start  = 1,
               End    = seq_length,
               Strand = "*") %>%
  readr::write_tsv(paste0(library_name, ".saf"))

