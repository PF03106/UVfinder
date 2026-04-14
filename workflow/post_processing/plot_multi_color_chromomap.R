#!/usr/bin/env Rscript
# Users should create sex_linked_genes.txt file as a tab-separated file 
# with Gene ID, Count, and Samples (option -g or --genes).

library(optparse)
library(tidyverse)
library(Biostrings)
library(chromoMap)
library(htmlwidgets)
library(data.table)

# -------------------------------
# 1. Set path
# -------------------------------
genome_base <- "/blue/mcdaniel/seyeonkim/UVfinder/results/00_renamed"
blast_base <- "/blue/mcdaniel/seyeonkim/UVfinder/results/03_locus_search"
meta_file <- "/blue/mcdaniel/seyeonkim/UVfinder/config/samples.tsv"
out_dir <- "/blue/mcdaniel/seyeonkim/UVfinder/results/Hypanales_chromomap" 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# 2. Set command-line options
# -------------------------------
option_list = list(
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="Path to gene list with counts (Tab separated: GeneID, Count, Samples)", metavar="character"),
  make_option(c("-m", "--mode"), type="character", default="all", 
              help="'Choose among 'all', 'order', 'sample_id'", metavar="character"),
  make_option(c("-t", "--target"), type="character", default=NULL, 
              help="If mode is set to 'order' or 'sample_id', target name (e.g. Dicriales or S001, S002, ... sep = , ", metavar="character"),
  make_option(c("-b", "--blast_type"), type="character", default="best_hits", 
              help="Choose between best_hits and all_hits", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$genes)){
  print_help(opt_parser)
  stop("Gene list must be specified", call.=FALSE)
}

# -------------------------------
# 3. Read multi-column gene list and prepare target samples
# -------------------------------
gene_info <- read.table(opt$genes, header=FALSE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
colnames(gene_info)[1:2] <- c("Gene", "Count") 
gene_info$Gene <- trimws(gene_info$Gene)

meta_df <- read.table(meta_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
target_samples <- c()   

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
  stop("No samples match the specified criteria.")
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
  if (!file.exists(ref_fa[1])) {
    message("❌ [ERROR] Reference genome file not found for sample -> skipping.")
    next
  }
  
  ref <- readDNAStringSet(ref_fa)
  chrlen_map <- setNames(width(ref), names(ref))
  
  chrom_file <- file.path(out_dir, paste0(sample, "_chrom.txt"))
  data.frame(
    chrom = names(chrlen_map),
    start = 1,
    end = as.integer(chrlen_map)
  ) %>% write_tsv(chrom_file, col_names = FALSE)

  blast_dir <- file.path(blast_base, sample)
  target_file_name <- ifelse(opt$blast_type == "all_hits", "A_All_hits.tsv", "B_Best_hits.tsv")
  blast_file <- file.path(blast_dir, target_file_name)

  dt <- fread(blast_file, header=TRUE, col.names = std_cols, sep="\t", fill=TRUE)
  if (nrow(dt) == 0) {
    warning("BLAST result is empty -> skipping.")
    next
  } else {
    all_hits <- dt %>% mutate(
      start = pmin(sstart, send),
      end   = pmax(sstart, send),
      probe = qseqid
    )
  }
  
  # -------------------------------
  # 5. Prepare annotation file (Modified with Dummy Lock Hack)
  # -------------------------------
  annot_df <- all_hits %>% 
    arrange(sseqid, start, end) %>%
    mutate(CleanProbe = trimws(sub("\\|.*", "", probe))) %>% 
    left_join(gene_info, by = c("CleanProbe" = "Gene")) %>%
    transmute(
      ElementName    = probe,
      ChromosomeName = sseqid,
      Start          = as.integer(start),
      End            = as.integer(end),
      Category       = ifelse(!is.na(Count), 
                              sprintf("Shared_%02d", as.integer(Count)), 
                              "Z_Other")
    ) %>% 
    filter(ChromosomeName %in% names(chrlen_map))

  if(nrow(annot_df) == 0) {
    warning("No annotations left -> skipping.")
    next
  }

  # -------------------------------
  # [Core Fix] Insert hidden dummy data to prevent chromoMap color shifting
  # -------------------------------
  # Create a fully defined color map
  color_map <- c(
    "Shared_01" = "#FFD166",
    "Shared_02" = "#B4CE84",
    "Shared_03" = "#48CAE4",
    "Shared_04" = "#0096C7",
    "Shared_11" = "#023E8A",
    "Z_Other"   = "#BDBDBD"
  )
  
  # Failsafe 1: Ensure color map is strictly alphabetically sorted
  color_map <- color_map[order(names(color_map))]

  # Failsafe 2: Create invisible 1bp dummy categories at the start of the 1st chromosome
  dummy_df <- data.frame(
    ElementName    = paste0("dummy_", names(color_map)),
    ChromosomeName = names(chrlen_map)[1],  
    Start          = 1L,
    End            = 1L,
    Category       = names(color_map),
    stringsAsFactors = FALSE
  )

  # Place dummy data at the very top of the actual dataset
  annot_df <- bind_rows(dummy_df, annot_df)

  annot_file <- file.path(out_dir, paste0(sample, "_annotation.txt"))
  write_tsv(annot_df, annot_file, col_names = FALSE)
  
  # -------------------------------
  # 6. Set Colors based on Count Category
  # -------------------------------
  # Pass the fixed, ordered 6 colors directly.
  sample_colors <- unname(color_map)

  # -------------------------------
  # 7. Visualize with chromoMap
  # -------------------------------
  title <- paste0(sample, "_", opt$blast_type, "_chromomap")
  html_file <- file.path(out_dir, paste0(title, ".html"))
  
  w <- chromoMap(
    ch.files = chrom_file,
    data.files = annot_file,
    data_based_color_map = TRUE,
    data_type = "categorical",
    data_colors = list(sample_colors),
    segment_ann = TRUE,
    chr_color   = "#f0f0f0",  
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