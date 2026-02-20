#!/usr/bin/env Rscript
# Users should create sex_linked_genes.txt file with one gene name per line (header: genes)(option -g or --genes) and specify mode (option -m or --mode) and target (option -t or --target) if needed.

library(optparse)
library(tidyverse)
library(Biostrings)
library(chromoMap)
library(htmlwidgets)
library(data.table)

# -------------------------------
# 1. Set path
# -------------------------------
genome_base <- "results/00_renamed"
blast_base <- "results/03_locus_search"
meta_file <- "config/samples.tsv"
out_dir <- "results/test_chromomap"
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
              help=" Choose between best_hit and all_hits", metavar="character")
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
              "qstart", "qend", "sstart", "send", "evalue", "bitscore", "rank", "sex_tag")

for (sample in target_samples) {
  message("Processing: ", sample)

  ref_fa <- file.path(genome_base, paste0(sample, "_renamed.fasta"))
  if (length(ref_fa) == 0 || is.na(ref_fa[1]) || !file.exists(ref_fa[1])) {
    message("❌ [ERROR] Reference genome file not found for sample: ", ref_fa, " -> skipping sample.")
    next
  }
  # create chrlen_map and chrom_file.
  ref <- readDNAStringSet(ref_fa)
  chrlen_map <- setNames(width(ref), names(ref))
  
  chrom_file <- file.path(out_dir, paste0(sample, "_chrom.txt"))
  data.frame(
    chrom = names(chrlen_map),
    start = 1,
    end = as.integer(chrlen_map)
  ) %>% write_tsv(chrom_file, col_names = FALSE)

  blast_dir <- file.path(blast_base, sample)

  target_file_name <- ifelse(opt$blast_type == "all_hits", 
                           "A_All_hits.tsv", 
                           "B_Best_hits.tsv")
  blast_file <- file.path(blast_dir, target_file_name)

  dt <- fread(blast_file, header=TRUE, col.names = std_cols, sep="\t", fill=TRUE)
  if (nrow(dt) == 0) {
    warning("BLAST result is empty for sample: ", sample, " -> skipping sample.")
    next
  } else {
    message("✔ BLAST results loaded: ", blast_file)
    all_hits <- dt %>% mutate(
      start = pmin(sstart, send),
      end   = pmax(sstart, send),
      probe = qseqid
    )
  }
  # -------------------------------
  # 5. Prepare annotation file for chromoMap (probe locations)
  # -------------------------------
  annot_df <- all_hits %>% 
    arrange(sseqid, start, end) %>%
    transmute(
      ElementName    = probe,
      ChromosomeName = sseqid,
      Start          = as.integer(start),
      End            = as.integer(end),
      Category       = ifelse(sub("\\|.*", "", probe) %in% sex_lniked_genes, "1_Target", "2_Other")
    ) %>% 
    filter(ChromosomeName %in% names(chrlen_map))

  if(nrow(annot_df) == 0) {
    warning("No annotations left after filtering for sample: ", sample, " -> skipping.")
    next
  }

  annot_file <- file.path(out_dir, paste0(sample, "_annotation.txt"))
  write_tsv(annot_df, annot_file, col_names = FALSE)
  
  present_categories <- sort(unique(annot_df$Category))
  ordered_categories <- unique(annot_df$Category)
  
  sample_colors <- c()
  for (cat in ordered_categories) {
    if (cat == "1_Target") {
      sample_colors <- c(sample_colors, "#FF0066")
    } else if (cat == "2_Other") {
      sample_colors <- c(sample_colors, "#0099CC")
    }
  }
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
    data_colors = list(sample_colors),
    segment_ann = TRUE,
    chr_color   = "#bdbdbd",
    ploidy      = 1,
    canvas_width  = 1400,
    canvas_height = 700,
    title = title, 
    export.options = list(
      format = "svg",
      filename = file.path(out_dir, paste0(title, ".svg"))
    )
  )
  
  saveWidget(w, file = html_file, selfcontained = TRUE)
  message("✔ HTML saved at: ", html_file)
}

message("All jobs completed.")