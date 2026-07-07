# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(svglite) 

# ---------------------------------------------------------
# 1. Load the data
# ---------------------------------------------------------
# Set the input and output file paths
input_file <- "/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/04_filtered/duplication_summary_by_order.tsv" 
output_file <- "/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/04_filtered/duplication_summary_stacked_plot.svg"

# Read the tab-separated file
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ---------------------------------------------------------
# 2. Filter and Select Necessary Columns
# ---------------------------------------------------------
# Select only the 'Order' and the specific duplication count columns
plot_data <- df %>%
  select(
    Order,
    Local_Dups,
    Local_Sex_Chrom_Dups,
    Intra_Chrom_Dups,
    Intra_Sex_Chrom_Dups,
    Inter_Chrom_Dups,
    Inter_Sex_Chrom_Dups,
    Sex_Linked_Dups,
    Potential_Inter_Chrom_Dups,
    Unclassified_Complex
  )

# ---------------------------------------------------------
# 3. Reshape Data (Wide to Long)
# ---------------------------------------------------------
# Pivot the data to long format for ggplot2 stacking
long_data <- plot_data %>%
  pivot_longer(
    cols = -Order, 
    names_to = "Duplication_Type",
    values_to = "Count"
  )

# Set factor levels to maintain a specific order in the legend and plot
long_data$Duplication_Type <- factor(long_data$Duplication_Type, levels = c(
  "Local_Dups", 
  "Local_Sex_Chrom_Dups",
  "Intra_Chrom_Dups", 
  "Intra_Sex_Chrom_Dups",
  "Inter_Chrom_Dups", 
  "Inter_Sex_Chrom_Dups",
  "Sex_Linked_Dups",
  "Potential_Inter_Chrom_Dups",
  "Unclassified_Complex"
))

# ---------------------------------------------------------
# 4. Generate the Stacked Bar Plot with Custom Colors
# ---------------------------------------------------------
# Define a custom color palette for each duplication type
my_colors <- c(
  "Local_Dups"                 = "#E97132", 
  "Local_Sex_Chrom_Dups"       = "#9aab4b", 
  "Intra_Chrom_Dups"           = "#febd2b", 
  "Intra_Sex_Chrom_Dups"       = "#d6828c", 
  "Inter_Chrom_Dups"           = "#0099CC", 
  "Inter_Sex_Chrom_Dups"       = "#3e9896", 
  "Sex_Linked_Dups"            = "#233341", 
  "Potential_Inter_Chrom_Dups" = "#CAEEFB", 
  "Unclassified_Complex"       = "#999999"
)

# Create the plot using scale_fill_manual for custom colors
p <- ggplot(long_data, aes(x = reorder(Order, -Count, sum), y = Count, fill = Duplication_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Distribution of Gene Duplication Types by Order",
    x = "Taxonomic Order",
    y = "Total Number of Duplications",
    fill = "Duplication Type"
  ) +
  scale_fill_manual(values = my_colors) + # Apply custom colors here
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# ---------------------------------------------------------
# 5. Save the Plot
# ---------------------------------------------------------
# Save the plot as a high-resolution SVG vector file
ggsave(output_file, plot = p, width = 12, height = 7)

print("Plot successfully saved as SVG!")