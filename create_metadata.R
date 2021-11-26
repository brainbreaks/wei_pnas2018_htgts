library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('breaktools/')

#
# Download FASTQ files
#
runs_df = readr::read_tsv("wei_pnas2018_runs.tsv")
lastchar = base::substring(runs_df$Run, nchar(runs_df$Run), nchar(runs_df$Run))
commands_download = with(runs_df, paste0("wget -O ", Run, "_1_", Description, ".fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR628/00", lastchar, "/", Run, "/", Run, "_1.fastq.gz\nwget -O ", Run, "_2_", Description, ".fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR628/00", lastchar, "/", Run, "/", Run, "_2.fastq.gz"))
writeLines(commands_download, con="download_ena.sh")

#
# Test which number of lines in each sample (0 if TLX is missing)
#
fasta_df = data.frame(fasta_path=Sys.glob("fastq/*_1_*.fastq.gz")) %>%
  dplyr::mutate(Run=gsub("_.*", "", basename(fasta_path))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(fasta_size=file.info(fasta_path)$size) %>%
  dplyr::ungroup()
results_df = data.frame(tlx_path=Sys.glob("*/results/*/*_result.tlx")) %>%
  dplyr::mutate(Run=gsub("_.*", "", basename(tlx_path))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(tlx_junctions=file_count_lines(tlx_path)) %>%
  dplyr::ungroup()
runs2results_df = runs_df %>%
  dplyr::left_join(results_df, by="Run") %>%
  dplyr::left_join(fasta_df, by="Run") %>%
  dplyr::mutate(tlx_junctions=tidyr::replace_na(tlx_junctions, 0)) %>%
  dplyr::arrange(tlx_junctions) %>%
  dplyr::mutate(i=1:n())
ggplot(runs2results_df, aes(x=i, y=tlx_junctions, color=Treatment)) +
  geom_label(aes(label=gsub("Chr", "", Chrom)), size=2) +
  ggrepel::geom_text_repel(aes(label=Run), data=runs2results_df %>% dplyr::filter(tlx_junctions==0), size=2) +
  scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
  coord_cartesian(ylim=c(0, 30)*1e3)


baits_df = readr::read_tsv("wei_pnas2018_baits.tsv")
baits_positions_df = htgts_calculate_positions(sequences_df=baits_df, database_path = "genomes/mm10/mm10.fa")

#
# Create Metadata files
#
htgts_metadata = runs_df %>%
  dplyr::group_by(Chrom, Treatment) %>%
  dplyr::mutate(Replicate=LETTERS[1:n()]) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(baits_positions_df, by="Chrom") %>%
  dplyr::mutate(
    Library=Run, Sequencing="WEI2018", Researcher="WP", Assembly="mm10", Chr=gsub("Chr", "chr", Chrom),
    Breakseq=NA_character_, Breaksite=NA_character_,
    MID=MID,
    Primer=RED_primer,
    Adapter="CCACGCGTGCCCTATAGTC",
    Cutter=NA_character_,
    Description=paste0(Chrom, "_", Treatment, "_", "Exp", Replicate)) %>%
  dplyr::select(Library, Sequencing, Researcher, Assembly, Chr, Start, End, Strand, Breakseq, Breaksite, MID, Primer, Adapter, Cutter, Description) %>%
  data.frame()
htgts_metadata %>%
  dplyr::rowwise() %>%
  dplyr::do(readr::write_tsv(as.data.frame(.), file=paste0("metadata/metadata_", .$Library, "_", .$Description, ".txt"), na=""))

#
# Create run script an separate directory for each chromosome
#
htgts_metadata %>%
  dplyr::filter(Chr=="chr14") %>%
  dplyr::group_by(Chr) %>%
  dplyr::do((function(htgts_group) {
    chr = htgts_group$Chr[1]
    cmd_prepare = paste0("mkdir ", chr, "; cd ", chr, "; rm htgts_latest.sif genomes metadata fastq; ln -s ../genomes genomes; ln -s ../metadata metadata; ln -s ../fastq fastq; ln -s ../htgts_latest.sif htgts_latest.sif")
    system(cmd_prepare)

    commands = with(htgts_group, paste0(
      "singularity exec -B /mnt/sda1 htgts_latest.sif TranslocPreprocess.pl metadata/metadata_", Library, "_", Description, ".txt preprocess --read1 fastq/", Library, "_1_", Description, ".fastq.gz --read2 fastq/", Library, "_2_", Description, ".fastq.gz\n",
      "singularity exec -B /mnt/sda1 --env BOWTIE2_INDEXES='genomes/mm10' -B `pwd` htgts_latest.sif TranslocWrapper.pl metadata/metadata_", Library, "_", Description,  ".txt preprocess results --threads 1"
    ))
    writeLines(commands, con=paste0(chr, "/run_chrom.sh"))
    data.frame()
  })(.))



#
# Extract MID for each sample
#
fastq_paths = list.files(path="~/Workspace/wei_pnas2018_htgts/fastq", pattern=".*_1_.*.fastq.gz", full.names=T)
barcodes_df = htgts_barcodes_detect(fastq_paths, primers=data.frame(Chrom=gsub(".*(Chr[^_]+).*", "\\1", fastq_paths)) %>% dplyr::left_join(baits_df) %>% .$RED_primer, max_sequences=10000)
ggplot(barcodes_df %>% dplyr::arrange(barcode_percent) %>% dplyr::mutate(i=1:n())) +
  ggrepel::geom_text_repel(aes(x=i, y=barcode_percent, label=gsub("_.*", "", basename(barcode_fasta)))) +
  geom_point(aes(x=i, y=barcode_percent))


sample_names = paste0(gsub("_.*", "", basename(fastq_paths)), "\n", gsub("_", "\n", gsub(".*_1_|\\.gz|\\.fastq|_Exp\\w", "", basename(fastq_paths))))
pp = plot_logos_coordinates(fastq_paths, sample_names, widths=list("Beginning"=c(1,35)))
pdf("reports/htgts_logo_all.pdf", width=8.27, height=10*11.69)
cowplot::plot_grid(plotlist=pp, ncol=1)
dev.off()





#
# Rename files
#
# htgts_files = htgts_metadata %>%
#   tidyr::crossing(data.frame(pair=c(1, 2))) %>%
#   dplyr::mutate(path_original=paste0("fastq/", Library, "_", pair, ".fastq"), path_new=paste0("fastq/", Library, "_", pair, "_", Description, ".fastq")) %>%
#   dplyr::select(path_original, path_new)
# writeLines(paste("mv", htgts_files$path_original, htgts_files$path_new))
# writeLines(paste("mv", htgts_files$path_new, htgts_files$path_original))
# fastq_paths = list.files(path="~/Workspace/wei_pnas2018_htgts/fastq", pattern=".*_[12]_.*.fastq", full.names=F)
# writeLines(paste0("seqtk seq -a ", fastq_paths, " > fasta/", gsub("\\.fastq", ".fa", fastq_paths)))

# fasta_paths = list.files(path="~/Workspace/wei_pnas2018_htgts/fastq_ncbi", pattern=".*_[12]_.*.fastq", full.names=F)
# writeLines(paste0("gzip ", fasta_paths))

# runs_new = readr::read_tsv("wei_pnas2018_runs.tsv") %>%
#   dplyr::select(-dplyr::matches("MID|Description")) %>%
#   dplyr::left_join(mid_df %>% dplyr::mutate(Description=gsub("SRR[^_]+_1_", "", sample_name)) %>% dplyr::select(Run, Description, MID) %>% dplyr::mutate(MID=ifelse(MID=="", NA_character_, MID)), by="Run")
# readr::write_tsv(runs_new, file="wei_pnas2018_runs.tsv", na="")