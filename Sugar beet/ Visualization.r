# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(purrr)

# Load your data
data <- read.csv("Data.csv")
data$Chl
# Choose the trait to plot
selected_trait <- "Chl"  # Replace with any other column name if needed

# Average replicates for each treatment and variety
data_avg <- data %>%
  group_by(Treatments, Varieties) %>%
  summarise(AverageValue = mean(get(selected_trait), na.rm = TRUE)) %>%
  ungroup()

# Perform LSD test for each trait and add LSD grouping
lsd_results <- data %>%
  pivot_longer(cols = -c(Treatments, Varieties), names_to = "Trait", values_to = "Value") %>%
  filter(Trait == selected_trait) %>%  # Filter for the selected trait only
  group_by(Trait) %>%
  nest() %>%
  mutate(lsd_results = map(data, ~ {
    model <- aov(Value ~ Treatments, data = .x)
    lsd <- LSD.test(model, "Treatments", group = TRUE)$groups
    lsd$Treatment <- rownames(lsd)
    return(lsd)
  })) %>%
  unnest(cols = lsd_results) %>%
  rename(Treatments = Treatment, lsd_group = groups) # Rename to match the main dataset

# Merge LSD results with averaged data for the plot
data_plot <- data_avg %>%
  left_join(lsd_results, by = c("Treatments"))

# Ensure Treatments is a factor to space out x-axis points
data_plot$Treatments <- factor(data_plot$Treatments, levels = unique(data_plot$Treatments))

# Customize colors for the varieties
custom_colors <- c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7")

# Line plot for the selected trait values across treatments with LSD grouping
ggplot(data_plot, aes(x = Treatments, y = AverageValue, color = Varieties, group = Varieties)) +
  geom_line(size = 1.5) +            # Lines connecting average values
  geom_point(size = 4) +             # Points representing average values
  geom_text(aes(label = lsd_group),  # Add LSD group annotations
            position = position_dodge(0.8), 
            vjust = -1.4, 
            size = 5) +              # Customize text size for clarity
  labs(title = "(A)", x = "Treatments", y = "Purity.") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title centered and bold
    axis.title.x = element_text(size = 14, face = "bold"),            # X-axis label formatting
    axis.title.y = element_text(size = 14, face = "bold"),            # Y-axis label formatting
    axis.text = element_text(size = 12),                              # Axis text size
    legend.position = "right"                                         # Position the legend
  ) +
  scale_color_manual(values = custom_colors)
