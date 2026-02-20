#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)
library(Biostrings)
library(chromoMap)
library(htmlwidgets)
library(data.table)

# -------------------------------
# 1. Set path
# -------------------------------
genome_base <- "resources/genomes"
blast_base <- "68_results/03_locus_search"
meta_file <- "config/samples.tsv"
out_dir <- "68_results/test_chromomap"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 2. Set command-line options
# -------------------------------
option_list = list(
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help=" Path to sex specific genes list(One gene at a line)", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default="all", 
              help="'Choose among 'all', 'order', 'sample_id'", metavar="character"),
  make_option(c("-t", "--target"), type="character", default=NULL, 
              help=" If mode is set to 'order' or 'sample_id', target name (e.g. Dicriales or S001, S002, ... sep = , ", metavar="character"),
  make_option(c("-b", "--blast_type"), type="character", default="best_hits", 
              help=" Choose between best_hits and all_hits", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$genes)){
  print_help(opt_parser)
  stop("Gene list must be specified", call.=FALSE)
}

# -------------------------------
# 3. Sex linked gene list (for red color) and target sample list preparation
# -------------------------------
# Sex linked gene list (color = red)
sex_lniked_genes <- readLines(opt$genes)
sex_lniked_genes <- sex_lniked_genes[sex_lniked_genes != ""]

meta_df <- read.table(meta_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
target_samples <- c()   # This will be the species to visualize in chromoMap

if (opt$mode == "all") {
  target_samples <- meta_df$sample_id
} else if (opt$mode == "order") {
  if (is.null(opt$target)) stop("If mode is 'order' --target must be specified")
  target_samples <- meta_df$sample_id[meta_df$order == opt$target]
} else if (opt$mode == "sample_id") {
  if (is.null(opt$target)) stop("If mode is 'sample_id' --target must be specified")
  target_samples <- strsplit(opt$target, ",")[[1]]
} else {
  stop("Unknown mode. Please choose among 'all', 'order', 'sample_id'")
}

if (length(target_samples) == 0) {
  stop("No samples match the specified criteria. Please check metadata and target name.")
}

message(sprintf("Creating chromoMap for %d samples...", length(target_samples)))

# -------------------------------
# 4. Loop through target samples to create chromoMap visualizations
# -------------------------------
std_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
              "qstart", "qend", "sstart", "send", "evalue", "bitscore")

for (sample in target_samples) {
  message("Processing: ", sample)

  ref_fa <- meta_df$genome_filename[meta_df$sample_id == sample]
  blast_dir <- file.path(blast_base, sample, opt$blast_type, )
  
  if (is.na(ref_fa) || !file.exists(ref_fa)) {
    message("❌ [ERROR] Reference genome file not found for sample: ", sample)
    message("   -> Looked for path: ", ref_fa)
    message("   -> Skipping this sample and moving to the next...")
    next
  }
  
  # read reference genome and create chromosome length map
  ref <- readDNAStringSet(ref_fa)
  chrlen_map <- setNames(width(ref), names(ref))
  
  # create chromoMap input file for chromosome information
  chrom_file <- file.path(out_dir, paste0(sample, "_chrom.txt"))
  data.frame(
    chrom = names(chrlen_map),
    start = 1,
    end = as.integer(chrlen_map)
  ) %>% write_tsv(chrom_file, col_names = FALSE)
  
  # Read BLAST files
  blast_files <- list.files(blast_dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(blast_files) == 0) {
    warning("BLAST result files not found: ", blast_dir, " -> skipping sample.")
    next
  }
  
  blast_list <- lapply(blast_files, function(f){
    dt <- fread(f, header=FALSE, col.names = std_cols, sep="\t", fill=TRUE)    
    if (nrow(dt) > 0) {
      dt %>% mutate(
        start = pmin(sstart, send),
        end   = pmax(sstart, send),
        probe = qseqid
      )
    }
  })
  
  all_hits <- rbindlist(blast_list)
  
  if (nrow(all_hits) == 0) {
    warning("BLAST hits are 0 for sample: ", sample, " -> skipping sample.")
    next
  }
  
  # -------------------------------
  # 5. Prepare annotation file for chromoMap (probe locations)
  # -------------------------------
  annot_file <- file.path(out_dir, paste0(sample, "_annotation.txt"))
    all_hits %>% 
    arrange(sseqid, start, end) %>%
    transmute(
      ElementName    = probe,
      ChromosomeName = sseqid,
      Start          = as.integer(start),
      End            = as.integer(end),
      Category       = ifelse(probe %in% sex_lniked_genes, "1_Target", "2_Other")
    ) %>% 
    filter(ChromosomeName %in% names(chrlen_map)) %>%
    write_tsv(annot_file, col_names = FALSE)
  
  # -------------------------------
  # 6. Visualize with chromoMap and save svg
  # -------------------------------
  title <- paste0(sample, "_", opt$blast_type, "_chromomap")
  html_file <- file.path(out_dir, paste0(title, ".html"))
  
  # chromoMap
  w <- chromoMap(
    ch.files = chrom_file,
    data.files = annot_file,
    data_based_color_map = TRUE,
    data_type = "categorical",
    data_colors = list(c("red", "blue")),
    segment_ann = TRUE,
    chr_color   = "#bdbdbd",
    ploidy      = 1,
    canvas_width  = 1400,
    canvas_height = 700,
    title = title
  )
  
  saveWidget(w, file = html_file, selfcontained = TRUE)
  message("✔ HTML saved at: ", html_file)
}

message("All jobs completed.")