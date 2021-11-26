devtools::load_all('../Emily/breaktools/')
library(ggplot2)
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readr)
library(rtracklayer)
library(tidyr)
library(ggthemr)
library(ggfortify)

color_scheme = c("Chr5"="#1F78B4", "DMSO (Chr5)"="#A6CEE3", "Chr6"="#E31A1C", "DMSO (Chr6)"="#FB9A99")
samples_path = "samples.tsv"
samples_df = readr::read_tsv(samples_path) %>%
  # dplyr::mutate(path=file.path(dirname(samples_path), gsub("_result", "", path))) %>%
  # dplyr::mutate(path=file.path("results", gsub("_result.tlx", "", path), gsub("_result.tlx", "_result.tlx", path))) %>%
  dplyr::filter(!grepl("99", path)) %>%
  dplyr::mutate(group_name=paste0(ifelse(control, "DMSO (", ""), group, ifelse(control, ")", ""))) %>%
  dplyr::mutate(rowname=sample) %>%
  tibble::column_to_rownames("rowname")
tlx_df = tlx_read_many(samples_df)
tlx_df = tlx_mark_rand_chromosomes(tlx_df)
tlx_df = tlx_remove_rand_chromosomes(tlx_df)
tlx_df = tlx_mark_bait_chromosome(tlx_df)
tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)
baits_df = tlx_identify_baits(tlx_df, breaksite_size=sgRNA_length, genome_fasta="genomes/mm10/mm10.fa")
