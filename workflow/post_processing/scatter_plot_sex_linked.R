# 1. Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggrepel")) install.packages("ggrepel") # For non-overlapping labels

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# 2. Create the dataset from your provided table
# This converts your raw sample-level data into a data frame
raw_samples <- data.frame(
  Sample = c("S0001", "S0002", "S0003", "S0004", "S0005", "S0006", "S0007", "S0008", "S0009", "S0010", 
             "S0011", "S0012", "S0013", "S0014", "S0015", "S0016", "S0017", "S0018", "S0019", "S0020", 
             "S0021", "S0022", "S0023", "S0024", "S0025", "S0026", "S0027", "S0028", "S0029", "S0030", 
             "S0031", "S0032", "S0033", "S0034", "S0035", "S0036", "S0037", "S0038", "S0039", "S0040", 
             "S0041", "S0042", "S0043", "S0044", "S0045", "S0046", "S0047", "S0048", "S0049", "S0050", "S0051"),
  Order = c("Polytrichales", "Hypnales", "Andreaeales", "Hypnales", "Rhizogoniales", "Bryales", "Pottiales", "Bryoxiphiales", "Ditrichales", "Ditrichales", 
            "Hypnales", "Hypnales", "Hypnales", "Dicranales", "Grimmiales", "Hypnales", "Hypnales", "Hypnales", "Hypnales", "Hypnales", 
            "Hypnales", "Hypnales", "Orthotrichales", "Polytrichales", "Polytrichales", "Hypnales", "Dicranales", "Hypnales", "Polytrichales", "Hypnales", 
            "Hypnales", "Polytrichales", "Bartramiales", "Grimmiales", "Hypnodendrales", "Hypnales", "Hypnales", "Sphagnales", "Pottiales", "Sphagnales", 
            "Sphagnales", "Sphagnales", "Sphagnales", "Sphagnales", "Sphagnales", "Hypnales", "Sphagnales", "Hypnales", "Tetraphidales", "Hypnales", "Pottiales"),
  Sex_Linked_Genes = c("G6924", "G5892", "", "G5892", "G5892,G6036", "G4796,G5812,G5842,G5899,G5950,G6041,G6227,G7248", "G5018,G5138,G5271,G5427,G5644,G5816,G5841,G5853,G5913,G5922,G5945,G5958,G5960,G6000,G6038,G6065,G6366,G6384,G6406,G6883,G6947,G6956,G6961,G7029,G7128,G7135,G7572", "G4603,G4989,G5138,G6639,G6717,G6947", "G4992,G5116,G5280,G5339,G5464,G5596,G5753,G5815,G5849,G5892,G6016,G6164,G6384,G6393", "G5116,G5280,G5339,G5464,G5596,G5702,G5721,G5753,G5815,G5849,G6016,G6393,G6995,G7128", 
                       "", "", "", "G5634,G5949,G6405,G6570,G7279", "G4889,G5335,G5366,G6237,G6487", "G5892", "G4890,G4992,G5177,G5299,G5335,G5357,G5366,G5596,G5702,G5933,G5943,G5981,G6035,G6237,G6258,G6487,G6639,G6913,G6968,G7324", "", "", "G4890,G4942,G4992,G5177,G5299,G5335,G5357,G5366,G5596,G5702,G5892,G5933,G5943,G5981,G6035,G6036,G6237,G6258,G6487,G6639,G6913,G6947,G6968,G7324", 
                       "G4793,G4806,G5111,G5318,G5406,G5469,G5477,G5822,G5865,G5892,G6034,G6051,G6119,G6139,G6265,G6318,G6460,G6499,G6500,G6679,G7174", "", "G5892,G6924", "", "G4603,G4796,G4951,G4989,G5138,G5463,G5599,G5703,G5770,G5866,G5943,G5944,G6056,G6098,G6139,G6164,G6299,G6572,G6660,G6924,G6954,G6969", "G5892", "G5528,G5634,G5899,G5949,G6570,G7279", "G5892", "G4603,G4796,G4951,G4989,G5138,G5463,G5599,G5703,G5770,G5866,G5943,G5944,G6056,G6098,G6139,G6164,G6299,G6447,G6572,G6660,G6924,G6954,G6969", "", 
                       "G5892", "", "G4603,G4724,G4744,G4848,G4989,G5034,G5138,G5264,G5427,G5599,G5770,G5802,G5815,G5853,G5863,G5870,G5894,G5921,G5922,G5945,G5958,G5980,G6038,G6065,G6098,G6130,G6164,G6221,G6320,G6376,G6384,G6406,G6412,G6435,G6450,G6457,G6462,G6494,G6559,G6660,G6717,G6955,G6961,G6969,G6995,G7067,G7128", "", "G4992,G5163,G5177,G5357,G5842,G6318,G6420", "G4848,G4890,G4942,G5335,G5357,G5366,G5513,G5596,G5634,G5702,G5842,G5892,G5899,G5949,G5981,G6035,G6036,G6226,G6237,G6238,G6274,G6405,G6451,G6487,G6548,G6854,G6913,G6924,G6968,G6995,G7279", "G4471,G4527,G4893,G5032,G5090,G5116,G5430,G5528,G5639,G5892,G5974,G6036,G6299,G6363,G6404,G6420,G6528,G6538,G6933,G7577", "G4603", "G5116,G5857,G5866,G6016", "", 
                       "", "", "", "", "G4603", "G5892", "", "", "", "", "G5116,G5857,G5866,G6016")
)

# 3. Aggregate data by Order
# We calculate: 1) Number of samples per order, 2) Number of unique sex-linked genes per order
order_summary <- raw_samples %>%
  # Filter out empty genes and split comma-separated strings
  mutate(Gene_List = strsplit(as.character(Sex_Linked_Genes), ",")) %>%
  unnest(Gene_List) %>%
  filter(Gene_List != "") %>%
  group_by(Order) %>%
  summarise(
    Unique_Genes_Found = n_distinct(Gene_List),
    .groups = 'drop'
  )

# Add sample counts (including those with zero genes) to the summary
sample_counts <- raw_samples %>%
  group_by(Order) %>%
  summarise(Sample_Size = n(), .groups = 'drop')

final_summary <- left_join(sample_counts, order_summary, by = "Order") %>%
  # Replace NA with 0 for orders where no sex-linked genes were found
  mutate(Unique_Genes_Found = ifelse(is.na(Unique_Genes_Found), 0, Unique_Genes_Found))

# 4. Generate the Scatter Plot
scatter_plot <- ggplot(final_summary, aes(x = Sample_Size, y = Unique_Genes_Found)) +
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
    title = "Correlation: Sample Size vs. Sex-linked Gene Discovery",
    subtitle = "Aggregated by Moss Order (Total 51 Samples)",
    x = "Number of Samples Analyzed (Sampling Effort)",
    y = "Total Unique Sex-linked Genes Found",
    caption = "Source: UVfinder analysis on GoFlag408 dataset"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none", # Hide legend since labels are on points
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 11)
  )

# 5. Save and Print
print(scatter_plot)
# ggsave("Order_Level_ScatterPlot.svg", scatter_plot, width = 8, height = 6)

save_name <- "sex_linked_scatter_plot.svg"
Results_dir <- "results/04_filtered"
save_path <- file.path(Results_dir, save_name)

ggsave(save_path, scatter_plot, width = 10, height = 7, device = "svg")
 