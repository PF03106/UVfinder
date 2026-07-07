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
genome_base <- "/blue/mcdaniel/seyeonkim/UVfinder/results/00_renamed"
blast_base <- "/blue/mcdaniel/seyeonkim/UVfinder/results/03_locus_search"
meta_file <- "/blue/mcdaniel/seyeonkim/UVfinder/config/samples.tsv"
out_dir <- "/blue/mcdaniel/seyeonkim/UVfinder/results/hypnales_chromomap_S0021" 
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
# Read and perfectly clean the sex-linked gene list (remove hidden '\r' and whitespaces)
sex_linked_genes <- readLines(opt$genes)
sex_linked_genes <- unlist(strsplit(sex_linked_genes, ","))
sex_linked_genes <- trimws(sex_linked_genes)
sex_linked_genes <- sex_linked_genes[sex_linked_genes != ""]

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
  
  # ==============================================================
  # [수정] 필터링하기 전에 유전체 파일과 BLAST 결과를 먼저 불러옵니다.
  # ==============================================================
  # 1. FASTA 파일 로드 및 시퀀스 이름 추출
  ref <- readDNAStringSet(ref_fa)
  all_seq_names <- names(ref)

  # 2. BLAST 결과 파일 로드
  blast_dir <- file.path(blast_base, sample)
  target_file_name <- ifelse(opt$blast_type == "all_hits", 
                             "A_All_hits.tsv", 
                             "B_Best_hits.tsv")
  blast_file <- file.path(blast_dir, target_file_name)

  if (file.exists(blast_file)) {
    dt <- fread(blast_file, header=TRUE, col.names = std_cols, sep="\t", fill=TRUE)
  } else {
    dt <- data.table()
  }

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

  # --- Filtering Logic ---
  # 이제 all_seq_names와 dt가 존재하므로 정상 작동합니다.
  
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

  chrom_file <- file.path(out_dir, paste0(sample, "_chrom.txt"))
  data.frame(
    chrom = names(chrlen_map),
    start = 1,
    end = as.integer(chrlen_map)
  ) %>% write_tsv(chrom_file, col_names = FALSE)

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
      Category       = ifelse(trimws(sub("\\|.*", "", probe)) %in% sex_linked_genes, "1_Target", "2_Other")
    ) %>% 
    filter(ChromosomeName %in% names(chrlen_map))

  if(nrow(annot_df) == 0) {
    warning("No annotations left after filtering for sample: ", sample, " -> skipping.")
    next
  }

  # ==============================================================
  # [수정] 색상 꼬임 현상 방지 (더미 데이터 삽입)
  # ==============================================================
  # chromoMap은 데이터에 존재하는 카테고리의 '알파벳 순서대로' 색상을 부여합니다.
  # 특정 샘플에 '1_Target'이 하나도 없을 경우 '2_Other'가 타겟 색상(핑크)으로 변하는 것을 막아줍니다.
  color_map <- c("1_Target" = "#FF0066", "2_Other" = "#0099CC")
  
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
  
  # 색상은 항상 고정된 순서로 적용됩니다.
  sample_colors <- unname(color_map)
  
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