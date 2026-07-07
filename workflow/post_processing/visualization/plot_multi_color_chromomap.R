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
out_dir <- "/blue/mcdaniel/seyeonkim/UVfinder/results/hypnales_Dups_chromomap" 
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
  
  # Load the entire genome first
  ref <- readDNAStringSet(ref_fa)
  all_seq_names <- names(ref)

  # Load BLAST results
  blast_dir <- file.path(blast_base, sample)
  target_file_name <- ifelse(opt$blast_type == "all_hits", "A_All_hits.tsv", "B_Best_hits.tsv")
  blast_file <- file.path(blast_dir, target_file_name)

  if (file.exists(blast_file)) {
    dt <- fread(blast_file, header=TRUE, col.names = std_cols, sep="\t", fill=TRUE)
  } else {
    dt <- data.table()
  }

  # --- Filtering Logic ---
  # 1. Identify sequence names ending in 'Chr[0-9UV]+'
  chr_names <- all_seq_names[grepl("Chr[0-9UV]+$", all_seq_names)]
  
  # 2. Identify sequence names that have at least one BLAST hit
  hit_seq_names <- unique(dt$sseqid)
  
  # 3. Final selection: All Chr + Scaffolds with hits (that exist in ref)
  target_seq_names <- unique(c(chr_names, intersect(all_seq_names, hit_seq_names)))
  
  if (length(target_seq_names) == 0) {
    message("⚠️ [SKIP] No Chr or Scaffolds with hits found for: ", sample)
    next
  }

  # Filter reference to include only target sequences
  ref_filtered <- ref[target_seq_names]
  chrlen_map <- setNames(width(ref_filtered), names(ref_filtered))
  
  # Create chrom_file
  chrom_file <- file.path(out_dir, paste0(sample, "_chrom.txt"))
  data.frame(
    chrom = names(chrlen_map),
    start = 1,
    end = as.integer(chrlen_map)
  ) %>% write_tsv(chrom_file, col_names = FALSE)

  # Prepare annotation data
  if (nrow(dt) > 0) {
    all_hits <- dt %>% mutate(
      start = pmin(sstart, send),
      end   = pmax(sstart, send),
      probe = qseqid
    )
    
    annot_df <- all_hits %>% 
      filter(sseqid %in% target_seq_names) %>%
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
      )
  } else {
    annot_df <- data.frame(ElementName=character(), ChromosomeName=character(), 
                           Start=integer(), End=integer(), Category=character())
  }

  # -------------------------------
  # Insert hidden dummy data to prevent chromoMap color shifting
  # -------------------------------
  color_map <- c(
    "Shared_01" = "#FFD166",
    "Shared_02" = "#B4CE84",
    "Shared_03" = "#48CAE4",
    "Shared_04" = "#0096C7",
    "Shared_11" = "#023E8A",
    "Z_Other"   = "#BDBDBD" 
  )
  
  color_map <- color_map[order(names(color_map))]

  # Dummy data on the first selected sequence (usually a Chr)
  dummy_df <- data.frame(
    ElementName    = paste0("dummy_", names(color_map)),
    ChromosomeName = names(chrlen_map)[1],  
    Start          = 1L,
    End            = 1L,
    Category       = names(color_map),
    stringsAsFactors = FALSE
  )

  final_annot_df <- bind_rows(dummy_df, annot_df)
  annot_file <- file.path(out_dir, paste0(sample, "_annotation.txt"))
  write_tsv(final_annot_df, annot_file, col_names = FALSE)
  
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
    canvas_height = 800, # Increased height slightly for more bars
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