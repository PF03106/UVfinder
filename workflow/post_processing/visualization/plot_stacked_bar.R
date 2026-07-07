# Load necessary packages
library(tidyverse) # For data manipulation and ggplot2
library(svglite)   # For exporting high-quality SVG files

# 1. Load the data 
# Assuming the file is tab-separated and named "data.txt" in your working directory.
# If it is a CSV, replace read_tsv with read_csv.
data <- read_tsv("/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/04_filtered/sex_linked_summary_by_order.tsv", show_col_types = FALSE)

# 2. Calculate the total number of unique samples for each Order
# This will be used to create the X-axis labels (e.g., "Order\n(N)").
sample_counts <- data %>%
  separate_rows(Samples, sep = ",") %>% # Split comma-separated samples (e.g., "S0014,S0027") into multiple rows
  group_by(Order) %>%
  summarise(Total_Samples = n_distinct(Samples), .groups = 'drop') # Count unique samples per Order

# 3. Process data for plotting
plot_data <- data %>%
  # Categorize based on the 'Count' value
  mutate(Category = case_when(
    Count == 1 ~ "1 Sample (Unique)",
    Count == 2 ~ "2 Samples",
    Count >= 3 ~ "3+ Samples"
  )) %>%
  # Aggregate the number of genes for each Order and Category
  group_by(Order, Category) %>%
  summarise(Num_Genes = n(), .groups = 'drop') %>%
  # Merge with the sample counts calculated above
  left_join(sample_counts, by = "Order") %>%
  # Create the final X-axis label combining Order name and total samples
  mutate(Order_Label = paste0(Order, "\n(", Total_Samples, ")"))

# 4. Factor treatments for sorting
# Sort the X-axis so the Order with the highest total number of genes comes first
order_totals <- plot_data %>%
  group_by(Order_Label) %>%
  summarise(Total_Genes = sum(Num_Genes)) %>%
  arrange(desc(Total_Genes))

plot_data$Order_Label <- factor(plot_data$Order_Label, levels = order_totals$Order_Label)

# Sort the fill categories so the bars stack correctly 
# In ggplot, the first level goes at the bottom. We want 3+ at the bottom and 1 at the top.
plot_data$Category <- factor(plot_data$Category, levels = c("3+ Samples", "2 Samples", "1 Sample (Unique)"))

# 5. Build the plot with ggplot2
p <- ggplot(plot_data, aes(x = Order_Label, y = Num_Genes, fill = Category)) +
  # Use stat="identity" to plot the actual values, and adjust width for better visual spacing
  geom_bar(stat = "identity", width = 0.7) +
  
  # --- DESIGN AESTHETICS ---
  # Apply custom colors to match the specific palette from the reference image
  scale_fill_manual(
    values = c(
      "1 Sample (Unique)" = "#FFD166", 
      "2 Samples"         = "#A9C773",        
      "3+ Samples"        = "#48CAE4"        
    ),
    # Force the legend order to read logically from left to right: 1, 2, 3+
    breaks = c("1 Sample (Unique)", "2 Samples", "3+ Samples") 
  ) +
  labs(
    x = "Order\n(Number of Samples)",
    y = "Number of Sex-linked Genes",
    fill = NULL # Remove the legend title
  ) +
  
  # Start with a minimal theme to remove background gray boxes
  theme_minimal() +
  
  # Fine-tune the design elements to exactly match the target image
  theme(
    legend.position = "top", # Move legend to the top
    
    # Gridline design: keep horizontal major lines, remove vertical lines and minor gridlines for a clean look
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # Axis line design: Add a solid black line only for the X-axis to ground the bars
    axis.line.x = element_line(color = "black", size = 0.5), 
    
    # Axis text design: Rotate X-axis text 45 degrees so long Order names don't overlap
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    
    # Axis title design: Make titles bold and add margins so they don't crowd the tick labels
    axis.title.x = element_text(margin = margin(t = 10), size = 13, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10), size = 13, face = "bold"),
    
    # Legend text design
    legend.text = element_text(size = 10)
  )

# Print the plot to the Viewer
print(p)

# 6. Save as a high-resolution SVG file
# Specify dimensions that fit the rotated text and legend comfortably
ggsave("/blue/mcdaniel/seyeonkim/UVfinder/results_moss_new/04_filtered/sex_linked_summary_by_order.svg", plot = p, width = 10, height = 6.5, device = "svg")