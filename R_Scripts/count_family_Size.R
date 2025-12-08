# Install the tidyverse package if you haven't already
# install.packages("tidyverse")
rm(list = ls())
# Load the library
library(tidyverse)
# Define the directory where your files are located
input_dir <- "Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/umi_celegans_consensus"


fgbio_umi_familiy_size = read_tsv(file.path(input_dir, "fgbio-groupreadsbyumi-plot.tsv"))

fgbio_umi_familiy_size.perentage = fgbio_umi_familiy_size |> pivot_longer(2:16,names_to = "Family size",values_to = "Percentage of Sample") |> 
  mutate(`Family size` = as.numeric(`Family size`))



# Create the ggplot
family_size_plot <- ggplot(fgbio_umi_familiy_size.perentage, aes(x = `Family size`, y = `Percentage of Sample`, color = Sample, group = Sample)) + 
  geom_line(linewidth = 1) +  # Add the lines
  geom_point(size = 2) +      # Add points to show actual data
  
  # Format the y-axis to show percentages
  scale_y_continuous(breaks = seq(0,100,by =10)) +scale_x_continuous(breaks = seq(0,15,by =1))+
  
  # Add labels and a title
  labs(
    title = "Proportion of Reads by UMI Family Size",
    x = "Family Size",
    y = "Percentage of Sample",
    color = "Sample" # Legend title
  ) + coord_cartesian()+
  
  # Use a clean theme
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right"
  )

# Display the plot
print(family_size_plot)

# Optional: Save the plot to a file
ggsave("Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/umi_celegans_consensus/family_size_plot.svg", plot = family_size_plot, width = 5, height = 3)




# List the full paths to your six TSV files
file_paths <- c(
  file.path(input_dir, "N2.30min.HS.1.family_size_summary.tsv"),
  file.path(input_dir, "N2.30min.HS.2.family_size_summary.tsv"),
  file.path(input_dir, "N2.30min.HS.3.family_size_summary.tsv"),
  file.path(input_dir, "PRDE1.30min.HS.1.family_size_summary.tsv"),
  file.path(input_dir, "PRDE1.30min.HS.2.family_size_summary.tsv"),
  file.path(input_dir, "PRDE1.30min.HS.3.family_size_summary.tsv")
)

# Create a clear name for each sample (for the plot legend)
sample_names <- c(
  "N2 HS 1", "N2 HS 2", "N2 HS 3",
  "PRDE1 HS 1", "PRDE1 HS 2", "PRDE1 HS 3"
)

# Read all files into a single, combined data frame
# The `map_dfr` function automatically adds a 'sample' column based on the names
all_samples_df <- map_dfr(set_names(file_paths, sample_names), read_tsv, .id = "sample")

# Display the first few rows of the combined data
head(all_samples_df)

total_family_size1 = all_samples_df |> filter(family_size==1) |> mutate(num_reads = formatC(num_reads, format = "e", digits = 2))
write.csv(total_family_size1,file = "../Data/output/UMI_data/Tables_UMI/total_family_size1.csv")

# It's often hard to see the trend with very large family sizes,
# so we'll filter to a more reasonable range for visualization.
plot_df <- all_samples_df %>%
  filter(family_size <= 10) # Adjust this value as needed

# Create the ggplot
family_size_plot <- ggplot(plot_df, aes(x = family_size, y = proportion, color = sample, group = sample)) +
  geom_line(linewidth = 1) +  # Add the lines
  geom_point(size = 2) +      # Add points to show actual data
  
  # Format the y-axis to show percentages
  scale_y_continuous(labels = scales::percent_format(),breaks = seq(0,1,by =0.1)) +scale_x_continuous(breaks = seq(0,10,by =1))+
  
  # Add labels and a title
  labs(
    title = "Proportion of Reads by UMI Family Size",
    x = "Family Size",
    y = "Proportion of Total Reads",
    color = "Sample" # Legend title
  ) + coord_cartesian(ylim = c(0,1.0))+
  
  # Use a clean theme
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right"
  )

# Display the plot
print(family_size_plot)

# Optional: Save the plot to a file
 ggsave("Y:/Johnny/Roswell_Projects/Sehee_mRNA_piRNA/UMI_Celegans/Data/output/BROAD/family_size_proportions.2.svg", plot = family_size_plot, width = 5, height = 3)
 
 