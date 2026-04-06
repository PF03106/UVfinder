# 1. Check and install missing packages
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("svglite")) install.packages("svglite") # Needed for SVG output

library(ggplot2)
library(dplyr)

# 2. Define paths
Results_dir <- "results/04_filtered"
file_path <- file.path(Results_dir, "sex_linked_summary_by_order.tsv")

# Create the directory if it doesn't exist
if (!dir.exists(Results_dir)) dir.create(Results_dir, recursive = TRUE)

# 3. Load and process data
if (file.exists(file_path)) {
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  stop(paste("INPUT FILE NOT FOUND at:", file_path))
}

# (Data processing part - identical to your previous code)
df_processed <- df %>%
  mutate(Category = case_when(
    Count == 1 ~ "1 Sample (Unique)",
    Count == 2 ~ "2 Samples (Shared)",
    Count >= 3 ~ "3+ Samples (Highly Shared)"
  )) %>%
  mutate(Category = factor(Category, 
                           levels = c("1 Sample (Unique)", 
                                      "2 Samples (Shared)", 
                                      "3+ Samples (Highly Shared)"))) %>%
  group_by(Order, Category) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Order) %>%
  mutate(Total = sum(n)) %>%
  ungroup() %>%
  mutate(Order = reorder(Order, -Total))

# 4. Create Plot
p <- ggplot(df_processed, aes(x = Order, y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("1 Sample (Unique)" = "#E69F00", 
                               "2 Samples (Shared)" = "#56B4E9", 
                               "3+ Samples (Highly Shared)" = "#009E73")) +
  labs(
    title = "Distribution of Sex-linked Loci across Moss Orders",
    x = "Taxonomic Order",
    y = "Number of Sex-linked Genes",
    fill = "Evidence Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

save_name <- "sex_linked_summary_by_order.svg"
save_path <- file.path(Results_dir, save_name)

ggsave(save_path, p, width = 10, height = 7, device = "svg")

cat("Success! Plot saved to:", save_path, "\n")