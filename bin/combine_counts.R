#!/usr/bin/env Rscript

################################################################################
# combine processed counts by summing counts by sample name across lanes
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp, extended and curated by Florian Andersch
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2023/03/07
################################################################################

### command line parameters
args         <- commandArgs(trailingOnly = TRUE)
library_file <- args[1]
count_files  <- args[2:length(args)]

### functions
`%>%` <- dplyr::`%>%`

read_featurecounts <- function(path) {
  readr::read_tsv(path, comment = "#") %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
    dplyr::rename(id = Geneid)
}

### combine counts
library <- readr::read_tsv(library_file) %>%
  dplyr::select(id, group)

names(count_files) <- stringr::str_replace(basename(count_files), ".txt", "")
pattern <- paste(c(paste0(names(count_files), "_"), ".sam"), collapse = "|")

counts <- lapply(count_files, read_featurecounts) %>%
  purrr::map(tidyr::gather, sample_name, count, -id) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sample_name = stringr::str_replace_all(sample_name, pattern, "")) %>%
  dplyr::group_by(id, sample_name) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(sample_name, count) %>%
  dplyr::inner_join(library, by = "id") %>%
  dplyr::select(id, group, dplyr::everything())

# counts[rowMeans(as.matrix(counts[,3:ncol(counts)]))>5,]
counts %>%
  readr::format_tsv() %>%
  cat
