# 0. Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggrepel")) install.packages("ggrepel") # For non-overlapping labels

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# 1. Set directory and file path
Results_dir <- "results/04_filtered"
input_file_name <- "sex_linked_species_lists.tsv"
file_path <- file.path(Results_dir, input_file_name)

save_name <- "sex_linked_scatter_plot.svg"  # Modify if necessary
save_path <- file.path(Results_dir, save_name)

# 2. Load and process data
if (file.exists(file_path)) {
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = c("", "NA"))
} else {
  stop(paste("INPUT FILE NOT FOUND at:", file_path))
}

# 3. Aggregate data by Order
final_summary <- df %>%
  group_by(Order) %>%
  summarise(
    Sample_Size = n(),                                 # 1) Number or taxa(sample) in each Order
    Total_Genes_Found = sum(Gene_Count, na.rm = TRUE), # 2) Total number of sex-linked genes found in each Order
    .groups = 'drop'
  )

# 4. Generate the Scatter Plot
scatter_plot <- ggplot(final_summary, aes(x = Sample_Size, y = Total_Genes_Found)) +
  # Add a linear regression line to show the general trend (Sampling Effort vs Result)
  geom_smooth(method = "lm", color = "lightgrey", fill = "whitesmoke", linetype = "dashed", alpha = 0.5) +
  # Add points colored by Order
  geom_point(aes(color = Order), size = 4, alpha = 0.8) +
  # Add non-overlapping labels
  geom_text_repel(aes(label = Order), 
                  size = 3.5, 
                  fontface = "italic", 
                  box.padding = 0.5) +
  # Labs and styling
  labs(
    title = "Sample Size vs. Number of Sex-linked Gene",
    x = "Number of Samples(Orders)",
    y = "Total Sex-linked Genes Found",
    color = "Order",
  ) +
  theme_minimal() +
  theme(
    legend.position = "right", 
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

# 5. Save and Print
cat("Saving scatter plot to:", save_path, "\n")
ggsave(save_path, scatter_plot, width = 10, height = 7, device = "svg")