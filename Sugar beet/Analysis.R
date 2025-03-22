library(openxlsx)
# Setting the Working Directory 
setwd("/home/raafat/Documents/Sugar_Beet")

# Load Agronomic Data 
data_beat <- read.csv("Data.csv", header=T)

factor_columns <- c("Treatmemts", "Varieties")
traits <- setdiff(names(data_beat), factor_columns)  # Dynamically detect traits

#Exploratory Data Analysis
hist(data_beat$Chl, col="blue",xlab="Plant Height (cm)", vlab="Frequency")

attach(data_beat)

boxplot(Chl~Treatmemts*Varieties, col="red", main="Yield by Population", xlab="Synthetic Population",ylab="Plant Height (cm)")

# Create a new Excel workbook
wb <- createWorkbook()

# Loop through each trait, perform ANOVA, and save the results
for (trait in traits) {
  
  # Formula for ANOVA
  formula <- as.formula(paste(trait, "~ Varieties * Treatmemts"))
  
  # Perform ANOVA
  aov_result <- summary(aov(formula, data = data_beat))
  
  # Convert the ANOVA table to a data frame
  anova_table <- as.data.frame(aov_result[[1]])
  
  # Add a worksheet to the Excel file for each trait
  addWorksheet(wb, trait)
  
  # Write the ANOVA table to the Excel sheet
  writeData(wb, trait, anova_table)
}
# Save the workbook to an Excel file
saveWorkbook(wb, "D2_ANOVA_Table.xlsx", overwrite = TRUE)
# Output message
cat("ANOVA results have been saved to anova_results.xlsx")



# Initialize an empty list to store formatted ANOVA results
anova_results <- list()

# Loop through each trait and perform ANOVA
for (trait in traits) {
  
  # Formula for ANOVA
  formula <- as.formula(paste(trait, "~ Varieties * Treatmemts"))
  
  # Perform ANOVA
  aov_result <- summary(aov(formula, data = data))
  
  # Extract ANOVA table
  anova_table <- as.data.frame(aov_result[[1]])
  
  # Rename columns to match your example
  colnames(anova_table) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  
  # Add a column for significance symbols based on p-values
  anova_table$Significance <- ifelse(anova_table$`Pr(>F)` < 0.001, "***",
                                     ifelse(anova_table$`Pr(>F)` < 0.01, "**",
                                            ifelse(anova_table$`Pr(>F)` < 0.05, "*", "ns")))
  
  # Add the Trait name at the top of the table to identify the trait in the final Excel sheet
  trait_header <- data.frame(Df = NA, Sum.Sq = NA, Mean.Sq = NA, `F value` = NA, `Pr(>F)` = NA, Significance = trait)
  
  # Ensure column names in the header match the ANOVA table
  colnames(trait_header) <- colnames(anova_table)
  
  # Combine the trait header with the ANOVA table for this trait
  trait_anova <- rbind(trait_header, anova_table)
  
  # Append to the results list
  anova_results[[trait]] <- trait_anova
}

# Combine all ANOVA tables into one data frame
combined_anova <- do.call(rbind, anova_results)

# Create a new Excel workbook
wb <- createWorkbook()

# Add a worksheet to store the combined ANOVA results
addWorksheet(wb, "ANOVA D2 Results")

# Write the combined ANOVA table to the worksheet
writeData(wb, "ANOVA D2 Results", combined_anova)

# Save the workbook to an Excel file
saveWorkbook(wb, "ANOVA D2 Results.xlsx", overwrite = TRUE)

# Output message
cat("ANOVA results have been saved to formatted_anova_results.xlsx")

#############################################################################################
# 2 Projection
#############################################################################################

# Varieties * Treatmemts
# Load necessary libraries
library(ggplot2)
library(factoextra)  # For PCA visualization
library(dplyr)       # For data manipulation

data_beat$Varieties <- as.factor(data_beat$Varieties)
data_beat$Treatments <- as.factor(data_beat$Treatmemts)

factor_columns <- c("Varieties", "Treatments")
trait_data <- data_beat %>% select(-one_of(factor_columns))

pca_result <- prcomp(trait_data, scale. = TRUE)
summary(pca_result)
pca_scores <- as.data.frame(pca_result$x)


# Compute eigenvalues
eigenvalues <- (pca_result$sdev)^2

# Print eigenvalues
print(eigenvalues)

# Variance explained by each PC
explained_variance <- eigenvalues / sum(eigenvalues) * 100
print(explained_variance)

eigenvectors <- pca_result$rotation

print(eigenvectors)

# Create a summary table
summary_table <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = eigenvalues,
  Variance_Explained = explained_variance,
  Cumulative_Variance = cumsum(explained_variance)
)

# Print summary table
print(summary_table)

# Scree plot of eigenvalues
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 100))

install.packages("writexl")  # Run this if `writexl` is not installed
library(writexl)
# Compute eigenvalues
eigenvalues <- (pca_result$sdev)^2
explained_variance <- eigenvalues / sum(eigenvalues) * 100
cumulative_variance <- cumsum(explained_variance)

# Create eigenvalues table
eigenvalues_table <- data.frame(
  PC = paste0("PC", 1:length(eigenvalues)),
  Eigenvalue = eigenvalues,
  Variance_Explained = explained_variance,
  Cumulative_Variance = cumulative_variance
)

# Extract eigenvectors (loadings)
eigenvectors_table <- as.data.frame(pca_result$rotation)
eigenvectors_table <- cbind(Variable = rownames(eigenvectors_table), eigenvectors_table)
rownames(eigenvectors_table) <- NULL  # Reset row names

# Save to Excel
write_xlsx(
  list(
    Eigenvalues = eigenvalues_table,
    Eigenvectors = eigenvectors_table
  ),
  "PCA_Results.xlsx"
)

# Output file saved as "PCA_Results.xlsx" in your working directory











##############################################################################
cov(data_beat)
res.cov= cov(data_beat[,3:17])
round(res.cov,2)
eigen(res.cov)
# Combine PCA scores with genotype and treatment for plotting
pca_scores$Genotype <- data_beat$Genotype
pca_scores$Treatment <- data_beat$Treatment

# Create a PCA biplot using ggplot2 and color by Genotype and shape by Treatment
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Genotype, shape = Treatment)) +
  geom_point(size = 4) +  # Scatter plot with points
  theme_minimal() +       # Clean theme
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add dashed lines for the origin
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right") +  # Legend position
  scale_color_brewer(palette = "Set1")  # Color palette for Genotypes

fviz_pca_var(pca_result, col.var = "black", repel = TRUE)  # Arrows showing variable contributions







# Load necessary libraries
library(ggplot2)
library(factoextra)
library(ggrepel)     # For text label repelling (avoids overlap)
library(dplyr)       # For data manipulation
library(gridExtra)   # To combine multiple plots if needed

#Variance explained for the first two PCs
variance_explained <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)
xlab_text <- paste0("PC1 (", variance_explained[1], "% Variance)")
ylab_text <- paste0("PC2 (", variance_explained[2], "% Variance)")

# Create a publication-quality PCA biplot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  
  # Scatter plot with points (color for Genotype, shape for Treatment)
  geom_point(aes(color = Genotype, shape = Treatment), size = 4, alpha = 0.8) +
  
  # Add labels for samples (optional, can be omitted if too cluttered)
  geom_text_repel(aes(label = rownames(pca_scores)), 
                  size = 3, 
                  segment.color = 'grey50', 
                  max.overlaps = 10) +  # Limits overlaps
  
  # Customize the color and shape palettes
  scale_color_manual(values = c("red", "blue", "green", "purple","yellow"))+  # Add more colors if needed +  # Color for Genotype
  scale_shape_manual(values = c(16, 17, 15, 18)) +  # Custom shapes for Treatment
  
  # Add custom labels for axes with variance explained
  labs(title = "PCA Biplot of Genotypes and Treatments",
       x = xlab_text, y = ylab_text) +
  
  # Center the plot at the origin
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  
  # Customize the theme for a clean, publication-quality appearance
  theme_minimal(base_size = 15) +  # Base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove grid lines for cleaner look
    axis.line = element_line(color = "black")  # Add axis lines
  )

# Add the biplot arrows for variable contributions (traits)
biplot_arrows <- fviz_pca_var(pca_result, 
                              col.var = "black", 
                              repel = TRUE, 
                              labelsize = 5, 
                              arrowsize = 1.2) + 
  theme_minimal(base_size = 15) +  # Same clean theme
  labs(title = "Contribution of Traits") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )
biplot_arrows
# Combine the two plots: PCA plot and biplot arrows
grid.arrange(pca_plot, biplot_arrows, ncol = 4)




# Check the number of unique treatments
unique_treatments <- unique(pca_scores$Treatment)
length(unique_treatments)  # This will tell you how many shapes you need

# Add enough shape values based on the unique treatments
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  
  # Scatter plot with points (color for Genotype, shape for Treatment)
  geom_point(aes(color = Genotype, shape = Treatment), size = 4, alpha = 0.8) +
  
  # Add labels for samples (optional, can be omitted if too cluttered)
  geom_text_repel(aes(label = rownames(pca_scores)), 
                  size = 3, 
                  segment.color = 'grey50', 
                  max.overlaps = 10) +  # Limits overlaps
  
  # Customize the color and shape palettes
  scale_color_manual(values = c("red", "blue", "green", "purple", "yellow")) +  # Adjust colors as needed
  scale_shape_manual(values = c(16, 17, 15, 18, 19, 20)) +  # Add more shapes here based on the unique treatments
  
  # Add custom labels for axes with variance explained
  labs(title = "PCA Biplot of Genotypes and Treatments",
       x = xlab_text, y = ylab_text) +
  
  # Center the plot at the origin
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  
  # Customize the theme for a clean, publication-quality appearance
  theme_minimal(base_size = 15) +  # Base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "right",  # Adjust legend position
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove grid lines for cleaner look
    axis.line = element_line(color = "black")  # Add axis lines
  )

# Continue with the rest of your code...

pca_plot
# Install plotly if you don't have it
# install.packages("plotly")

# Load necessary libraries
library(ggplot2)
library(plotly)

# Assuming you already have a PCA result with at least 3 PCs
# Plot in 3D using plotly
pca_3d <- plot_ly(pca_scores, 
                  x = ~PC1, y = ~PC2, z = ~PC3, 
                  type = 'scatter3d', 
                  mode = 'markers', 
                  marker = list(size = 6, 
                                color = ~Genotype, 
                                colorscale = 'Viridis', 
                                symbol = ~Treatment, 
                                symbols = c('circle', 'square', 'diamond', 'cross', 'x'))) %>%
  
  # Add layout settings for labels and titles
  layout(title = "3D PCA Biplot of Genotypes and Treatments",
         scene = list(xaxis = list(title = xlab_text), 
                      yaxis = list(title = ylab_text), 
                      zaxis = list(title = "PC3")),
         legend = list(title = list(text = 'Genotype')))

# Display the plot
pca_3d

# Check the structure of the pca_scores data
str(pca_scores)

# Ensure you have PC1, PC2, and PC3
head(pca_scores)
# Load necessary libraries
library(plotly)

# Simplified 3D PCA plot
pca_3d <- plot_ly(
  data = pca_scores, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(size = 6)
)

# Display the plot
pca_3d
# Plot with color and symbols
pca_3d <- plot_ly(
  data = pca_scores, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(size = 6, color = ~Genotype, symbol = ~Treatment)
)

# Display the plot
pca_3d
# Plot with custom layout
pca_3d <- plot_ly(
  data = pca_scores, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(size = 6, color = ~Genotype, colorscale = 'Viridis', symbol = ~Treatment)
) %>%
  layout(
    title = "3D PCA Biplot of Genotypes and Treatments",
    scene = list(
      xaxis = list(title = xlab_text), 
      yaxis = list(title = ylab_text), 
      zaxis = list(title = "PC3")
    )
  )

####################################################################
# Ensure Genotype and Treatment are factors
pca_scores$Genotype <- as.factor(pca_scores$Genotype)
pca_scores$Treatment <- as.factor(pca_scores$Treatment)

# Create a custom color scale for Genotype (number of colors must match the number of Genotypes)
genotype_colors <- c("red", "blue", "green", "purple", "yellow")

# Create a custom symbol mapping for Treatment (number of symbols must match the number of Treatments)
treatment_symbols <- c("circle", "square", "diamond", "cross", "x", "square-open")
# Check if the number of levels in Treatment matches the number of symbols
if (length(unique(pca_scores$Treatment)) > length(treatment_symbols)) {
  stop("You need more symbols to match the number of Treatment levels!")
}
# Plotly 3D PCA plot with manual color and symbol mapping
pca_3d <- plot_ly(
  data = pca_scores, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(
    size = 6, 
    color = genotype_colors[as.numeric(pca_scores$Genotype)],  # Assign colors manually
    symbol = treatment_symbols[as.numeric(pca_scores$Treatment)]  # Assign symbols manually
  )
)%>%
  layout(
    title = "3D PCA Biplot of Genotypes and Treatments",
    scene = list(
      xaxis = list(title = xlab_text), 
      yaxis = list(title = ylab_text), 
      zaxis = list(title = "PC3")
    )
  )

# Display the plot
pca_3d

#------------------------------------------------------------
library(plotly)

# Ensure Genotype and Treatment are factors
pca_scores$Genotype <- as.factor(pca_scores$Genotype)
pca_scores$Treatment <- as.factor(pca_scores$Treatment)

# Create a custom color scale for Genotype
genotype_colors <- c("red", "blue", "green", "purple", "yellow")
genotype_levels <- levels(pca_scores$Genotype)

# Create a custom symbol mapping for Treatment
treatment_symbols <- c("circle", "square", "diamond", "cross", "x", "star")
treatment_levels <- levels(pca_scores$Treatment)

# Check if the number of levels in Treatment matches the number of symbols
if (length(treatment_levels) > length(treatment_symbols)) {
  stop("You need more symbols to match the number of Treatment levels!")
}

# Define color and symbol mapping functions
get_color <- function(genotype) {
  genotype_colors[as.numeric(genotype)]
}

get_symbol <- function(treatment) {
  treatment_symbols[as.numeric(treatment)]
}

# Plotly 3D PCA plot with manual color and symbol mapping
pca_3d <- plot_ly(
  data = pca_scores, 
  x = ~PC1, y = ~PC2, z = ~PC3, 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(
    size = 6, 
    color = ~Genotype,  # Use Genotype to map color
    colorscale = genotype_colors,
    symbol = ~Treatment,  # Use Treatment to map symbol
    symbol = sapply(pca_scores$Treatment, get_symbol)
  ),
  text = ~paste("Genotype:", Genotype, "<br>Treatment:", Treatment),  # Hover info
  hoverinfo = "text"
) %>%
  layout(
    title = "3D PCA Biplot of Genotypes and Treatments",
    scene = list(
      xaxis = list(title = xlab_text), 
      yaxis = list(title = ylab_text), 
      zaxis = list(title = "PC3")
    ),
    legend = list(
      title = list(text = "Genotypes"),
      orientation = 'h',  # Horizontal legend
      yanchor = 'bottom',
      y = 1.1,  # Adjust position if needed
      xanchor = 'right',
      x = 1
    )
  ) %>%
  add_trace(
    type = "scatter3d",
    mode = "markers",
    x = pca_scores$PC1,
    y = pca_scores$PC2,
    z = pca_scores$PC3,
    marker = list(
      size = 6,
      color = sapply(pca_scores$Genotype, get_color),
      symbol = sapply(pca_scores$Treatment, get_symbol)
    ),
    text = ~paste("Genotype:", Genotype, "<br>Treatment:", Treatment),
    hoverinfo = "text"
  )

# Display the plot
pca_3d


#############################################################
# Heatmap
# Load necessary libraries
library(pheatmap)    # For creating the heatmap with clustering
library(ggplot2)     # For additional aesthetic tweaks if necessary
library(dplyr)       # For data manipulation
library(RColorBrewer) # For color palettes

# Load your data (replace 'your_data.csv' with the actual file)
# Assuming the data has columns for 'Genotype', 'Treatment', and trait values
# your_data <- read.csv("your_data.csv")

# Ensure 'Genotype' and 'Treatment' are treated as factors
data_beat$Genotype <- as.factor(data_beat$Genotype)
data_beat$Treatment <- as.factor(data_beat$Treatment)

# Extract only the trait columns (ignoring Genotype and Treatment) for clustering
factor_columns <- c("Genotype", "Treatment")
trait_data <- data_beat %>% select(-one_of(factor_columns))

# Set row names as the combined Genotype and Treatment for better identification in heatmap
rownames(trait_data) <- paste(data_beat$Genotype, data_beat$Treatment, sep = "_")

# Scale the trait data (normalize the traits for heatmap)
scaled_trait_data <- scale(trait_data)

# Customize color palette for the heatmap
heatmap_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# Adjust the plot window size by setting R graphical parameters
par(mar = c(1, 1, 1, 1))  # Adjust the plot margins if needed

# Split the heatmap by Treatment
split_treatments <- data_beat$Treatment

# Generate the heatmap with splitting by treatment
pheatmap(scaled_trait_data,
         color = heatmap_colors,        # Color palette for the heatmap
         clustering_distance_rows = "euclidean",  # Clustering method for rows
         clustering_distance_cols = "euclidean",  # Clustering method for columns
         clustering_method = "complete",  # Hierarchical clustering method
         fontsize = 10,                  # Font size for heatmap
         cellwidth = 20,                 # Adjusted cell width
         cellheight = 10,                # Adjusted cell height
         treeheight_row = 50,            # Height of the row dendrogram
         treeheight_col = 50,            # Height of the column dendrogram
         show_rownames = TRUE,           # Show row names (Genotype_Treatment)
         show_colnames = TRUE,           # Show column names (traits)
         border_color = NA,              # Remove borders around cells
         main = "Heatmap of Traits by Treatment",
         split = split_treatments        # Split heatmap by Treatment
)

# Reset graphical parameters if necessary after plotting
par(mar = c(5.1, 4.1, 4.1, 2.1))  # Reset to default margins after the plot













data_beat$Genotype <- as.factor(data_beat$Genotype)
data_beat$Treatment <- as.factor(data_beat$Treatment)

# Extract only the trait columns (ignoring Genotype and Treatment) for clustering
factor_columns <- c("Genotype", "Treatment")
trait_data <- data_beat %>% select(-one_of(factor_columns))

# Set row names as the combined Genotype and Treatment for better identification in heatmap
rownames(trait_data) <- paste(data_beat$Genotype, data_beat$Treatment, sep = "_")

# Scale the trait data (normalize the traits for heatmap)
scaled_trait_data <- scale(trait_data)

# Transpose the data to flip the heatmap
transposed_trait_data <- t(scaled_trait_data)

# Customize color palette for the heatmap
heatmap_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# Adjust the plot window size by setting R graphical parameters
par(mar = c(1, 1, 1, 1))  # Adjust the plot margins if needed

# Split the heatmap by Treatment
split_treatments <- data_beat$Treatment

# Generate the transposed heatmap with splitting by treatment
pheatmap(transposed_trait_data,
         color = heatmap_colors,        # Color palette for the heatmap
         clustering_distance_rows = "euclidean",  # Clustering method for rows
         clustering_distance_cols = "euclidean",  # Clustering method for columns
         clustering_method = "complete",  # Hierarchical clustering method
         fontsize = 10,                  # Font size for heatmap
         cellwidth = 10,                 # Adjusted cell width for transposed data
         cellheight = 20,                # Adjusted cell height for transposed data
         treeheight_row = 50,            # Height of the row dendrogram
         treeheight_col = 50,            # Height of the column dendrogram
         show_rownames = TRUE,           # Show row names (Traits)
         show_colnames = TRUE,           # Show column names (Genotype_Treatment)
         border_color = NA,              # Remove borders around cells
         main = "Transposed Heatmap of Traits by Treatment",
         split = split_treatments        # Split heatmap by Treatment
)

# Reset graphical parameters if necessary after plotting
par(mar = c(5.1, 4.1, 4.1, 2.1))  # Reset to default margins after the plot



library(pheatmap)    # For creating the heatmap with clustering
library(ggplot2)     # For additional aesthetic tweaks if necessary
library(dplyr)       # For data manipulation
library(RColorBrewer) # For color palettes

# The mtcars dataset:
data_beat <- read.csv("Sugar beat algea_Means.csv",row.names = 1, header=T)

heat <- as.matrix(data_beat)

# Default Heatmap
heatmap(heat, scale="column")# No dendrogram nor reordering for neither column or row
heatmap(heat, scale="column",
        clustering_distance_rows = "euclidean",  # Clustering method for rows
        clustering_distance_cols = NA,  # Clustering method for columns
        clustering_method = "complete")

# The mtcars dataset:
data <- as.matrix(mtcars)

# Default Heatmap
heatmap(data)










# Correlatioin Analysis #######################################################
library(corrplot)
data_beat <- read.csv("Sugar beat algea.csv", header=T)
MATr = data_beat[76:90,3:18]

C = cor(MATr)

testRes = cor.mtest(MATr, conf.level = 0.95)
corrplot(C,
         method="circle",
         type = "upper",
         number.font = 0.2,
         na.label = "NA",
         number.cex = 0.4,
         cl.length = 11,
         tl.cex=0.6,
         cl.pos = "b",
         tl.col = "black",
         addgrid.col = "grey",
         p.mat = testRes$p, sig.level = 0.05)

testRes = cor.mtest(MATr, conf.level = 0.95)

corrplot(C, p.mat = testRes$p, method = 'color', diag = FALSE, type = 'upper',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'grey20', order = 'AOE')

corrplot(C, p.mat = testRes$p, sig.level = 0.05, order = 'AOE', addrect = 2,type = 'upper')

corrplot.mixed(C,
               lower = 'shade', upper = 'pie',
               sig.level=0.05,p.mat = testRes$p,insig='blank',
               addCoef.col ='black', number.cex = 0.8, order = 'AOE')


corrplot(C, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)




# Install the necessary packages
install.packages("pheatmap")
library(pheatmap)
library(reshape2)

# Load your data
data <- read.csv("your_file_path.csv")

# Select a trait to visualize, e.g., "RY(ton-fed)"
heatmap_data <- dcast(data_beat, Genotype ~ Treatment, value.var = "Chl")

# Set rownames as Genotypes
rownames(heatmap_data) <- heatmap_data$Genotype
heatmap_data <- heatmap_data[, -1]

# Create a heatmap
pheatmap(heatmap_data, 
         cluster_rows = TRUE,  
         cluster_cols = TRUE,  
         display_numbers = TRUE, 
         scale = "row",  # Normalize rows if needed
         fontsize_row = 12, 
         fontsize_col = 12, 
         main = "Heatmap of RY(ton-fed) for Genotypes and Treatments",
         color = colorRampPalette(c("blue", "white", "red"))(50))

#######################################################################3
# Bar 
#######################################################################
# Load libraries
# Load libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(ggpubr)
library(tidyr) 

# Set working directory and load the data
df <- read.csv("Sugar beat algea.csv", header = TRUE)

# Melt the data frame (selecting relevant traits)
df_melt <- df %>% 
  melt(id.vars = c("Genotype", "Treatment"), 
       measure.vars = c("RY.ton.fed.", "Root.Weight", "Leaves.Weight"))

# Calculate summary statistics
summary_df <- df_melt %>%
  group_by(Genotype, Treatment, variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()))

# Function to create a plot for each trait
create_plot <- function(trait) {
  ggplot(summary_df %>% filter(variable == trait), aes(x = Genotype, y = mean_value, fill = Treatment)) + 
    geom_col(position = position_dodge(0.8), width = 0.8, color = "black", size = 0.7) + # Increase dodge width for spacing
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                  width = 0.2, position = position_dodge(0.8)) +  # Adjust dodge width for error bars
    labs(title = paste("Barplot of", trait), 
         x = "Genotype", y = trait) + 
    theme_classic() +  # Classic theme
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, family = "Times New Roman", face = "bold", size = 12),
      axis.text.y = element_text(family = "Times New Roman", face = "bold"), 
      axis.title.x = element_text(family = "Times New Roman", face = "bold"),  
      axis.title.y = element_text(family = "Times New Roman", face = "bold"),  
      plot.title = element_text(family = "Times New Roman", face = "bold"),  
      legend.text = element_text(family = "Times New Roman", face = "bold"),  
      legend.title = element_text(family = "Times New Roman", face = "bold")  
    )
}

# Create plots for each trait
plot1 <- create_plot("RY.ton.fed")
plot2 <- create_plot("Root.Weight")
plot3 <- create_plot("Leaves.Weight")

# Arrange plots in a grid
grid.arrange(plot1, plot2, plot3, ncol = 1)







# Load libraries
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(ggpubr)
library(tidyr)

# Set working directory and load the data
df <- read.csv("Sugar beat algea.csv", header = TRUE)

# Melt the data frame (selecting relevant traits)
df_melt <- df %>% 
  melt(id.vars = c("Genotype", "Treatment"), 
       measure.vars = c("Chl", "TSS", "Sucrose."))

# Calculate summary statistics
summary_df <- df_melt %>%
  group_by(Treatment, Genotype, variable) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            se_value = sd(value, na.rm = TRUE) / sqrt(n()))

# Function to create a plot for each trait
create_plot <- function(trait) {
  ggplot(summary_df %>% filter(variable == trait), aes(x = Treatment, y = mean_value, fill = Genotype)) + 
    geom_col(position = position_dodge(0.8), width = 0.8, color = "black", size = 0.7) + # Dodge width for spacing
    geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value), 
                  width = 0.2, position = position_dodge(0.8)) +  # Error bars
    labs(title = paste("Barplot of", trait), 
         x = "Treatment", y = trait) + 
    theme_classic() +  # Classic theme
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, family = "Times New Roman", face = "bold", size = 12),
      axis.text.y = element_text(family = "Times New Roman", face = "bold"), 
      axis.title.x = element_text(family = "Times New Roman", face = "bold"),  
      axis.title.y = element_text(family = "Times New Roman", face = "bold"),  
      plot.title = element_text(family = "Times New Roman", face = "bold"),  
      legend.text = element_text(family = "Times New Roman", face = "bold"),  
      legend.title = element_text(family = "Times New Roman", face = "bold")  
    )
}

# Create plots for each trait
plot1 <- create_plot("Chl")
plot2 <- create_plot("TSS")
plot3 <- create_plot("Sucrose")

# Arrange plots in a grid
grid.arrange(plot1, plot2, plot3, ncol = 1)



###############################################################################
######## 
###############################################################################
# Install necessary packages if not already installed
if (!requireNamespace("ggbump", quietly = TRUE)) install.packages("ggbump")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggthemes", quietly = TRUE)) install.packages("ggthemes")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(ggbump)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(tidyr)
data_beat <- read.csv("Data.csv", header=T)

# Load your data here (replace 'data' with your actual dataset)
data <- read.csv("Data.csv")  # Make sure your data is loaded correctly

# Dynamically pivot data for easy plotting using ggbump
data_long <- data %>%
  pivot_longer(cols = -c(Treatments, Varieties), 
               names_to = "Trait", values_to = "Value")

# Create a function to plot each trait
plot_trait <- function(trait) {
  ggplot(data_long %>% filter(Trait == trait), aes(x = Treatments, y = Value, color = Varieties, group = Varieties)) +
    geom_bump(size = 1.5) +  # Line bump chart
    geom_point(size = 4) +   # Points
    labs(title = paste("Variation in", trait, "across Treatments"), 
         x = "Treatments", y = trait) +
    scale_color_manual(values = c("#2E86C1", "#E74C3C")) + # Customize colors for varieties
    theme_few() +  # Theme from ggthemes for publication
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
}

# Extract all unique traits from the data and plot each
unique_traits <- unique(data_long$Trait)

# Generate a plot for each trait and store in a list
plots <- lapply(unique_traits, plot_trait)

# Display all plots in a grid (optional)
library(gridExtra)
grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol to fit your layout preference



#############################################33333333
# Required libraries
# Required libraries
library(ggbump)
library(ggplot2)
library(dplyr)

# Load your data
data <- read.csv("Data.csv")

# Choose a trait to rank varieties based on the averaged values
selected_trait <- "Leaves.weight..g."  # Replace with any other column name if needed

# Average replicates and rank varieties within each treatment
data_avg_ranked <- data %>%
  group_by(Treatments, Varieties) %>%
  summarise(AverageValue = mean(get(selected_trait), na.rm = TRUE)) %>%  # Average replicates
  group_by(Treatments) %>%
  mutate(Rank = rank(-AverageValue)) %>%  # Rank in descending order
  ungroup()

# Bump chart with ggplot and ggbump
ggplot(data_avg_ranked, aes(x = Treatments, y = Rank, color = Varieties, group = Varieties)) +
  geom_bump(size = 2) +       # Thick lines for bump effect
  geom_point(size = 6) +      # Larger points to match the style
  scale_y_reverse(breaks = 1:max(data_avg_ranked$Rank)) +  # Reverse y-axis to show rank from top to bottom
  labs(title = paste("Ranking of Varieties by Average", selected_trait), x = "Treatments", y = "Rank") +
  theme_minimal() +           # Minimal theme for a clean look
  theme(
    legend.position = "right",               # Adjust legend position if needed
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.text.y = element_text(size = 12),   # Customize y-axis text
    axis.text.x = element_text(size = 12)
  ) +
  scale_color_manual(values = c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7")) # Customize colors
# Required libraries
library(ggbump)
library(ggplot2)
library(dplyr)
setwd("/home/raafat/Documents/Sugar_Beet")

# Load your data
data <- read.csv("Data.csv")
data$Sucrose.
# Choose a trait to rank varieties based on the averaged values
selected_trait <- "Sucrose."  # Replace with any other column name if needed

# Average replicates and rank varieties within each treatment
data_avg_ranked <- data %>%
  group_by(Treatments, Varieties) %>%
  summarise(AverageValue = mean(get(selected_trait), na.rm = TRUE)) %>%  # Average replicates
  group_by(Treatments) %>%
  mutate(Rank = rank(-AverageValue)) %>%  # Rank in descending order
  ungroup()

# Customize colors for the varieties
custom_colors <- c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7")

# Bump chart with ggplot and ggbump
ggplot(data_avg_ranked, aes(x = Treatments, y = Rank, color = Varieties, group = Varieties)) +
  geom_bump(size = 2) +       # Thick lines for bump effect
  geom_point(size = 6) +      # Larger points to match the style
  scale_y_reverse(breaks = 1:max(data_avg_ranked$Rank)) +  # Reverse y-axis to show rank from top to bottom
  labs(title = paste("(A)"), x = "Treatments", y = "Sucrose.") +
  theme_minimal(base_size = 14) + # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title centered and bold
    axis.title.x = element_text(size = 14, face = "bold"), # X-axis label formatting
    axis.title.y = element_text(size = 14, face = "bold"), # Y-axis label formatting
    axis.text = element_text(size = 12),                  # Axis text size
    legend.position = "right",                             # Remove legend for simplicity
    plot.margin = unit(c(1,1,1,1), "lines")               # Set margin to make it more square
  ) +
  scale_color_manual(values = custom_colors) +
  coord_fixed(ratio = 1) # Set aspect ratio to 1:1 for a square plot



?geom_bump


# Required libraries
library(ggbump)
library(ggplot2)
library(dplyr)

# Load your data
data <- read.csv("Data.csv")
str(data$Chl)
# Choose the trait to plot
selected_trait <- "Sucrose."  # Replace with any other column name if needed

# Average replicates for each treatment and variety
data_avg <- data %>%
  group_by(Treatments, Varieties) %>%
  summarise(AverageValue = mean(get(selected_trait), na.rm = TRUE)) %>%
  ungroup()


# Perform LSD test and add LSD grouping
lsd_results <- data %>%
  pivot_longer(cols = -c(Treatments, Varieties), names_to = "Trait", values_to = "Value") %>%
  group_by(Trait) %>%
  nest() %>%
  mutate(lsd_results = map(data, ~ {
    model <- aov(Value ~ Treatments, data = .x)
    lsd <- LSD.test(model, "Treatments", group = TRUE)$groups
    lsd$Treatment <- rownames(lsd)
    return(lsd)
  })) %>%
  unnest(cols = lsd_results) %>%
  rename(lsd_group = groups)
# Ensure Treatments is a factor to space out x-axis points
data_avg$Treatments <- factor(data_avg$Treatments, levels = unique(data_avg$Treatments))

# Customize colors for the varieties
custom_colors <- c("#FF6F61", "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251", "#B565A7")

# Line plot for trait values across treatments
ggplot(data_avg, aes(x = Treatments, y = AverageValue, color = Varieties, group = Varieties)) +
  geom_line(size = 1.5) +          # Lines connecting average values
  geom_point(size = 4) +           # Points representing average values
  labs(title = "(A)", x = "Treatments", y = "Sucrose.") +
  theme_minimal(base_size = 14) +
  # Base font size for publication
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title centered and bold
    axis.title.x = element_text(size = 14, face = "bold"),            # X-axis label formatting
    axis.title.y = element_text(size = 14, face = "bold"),            # Y-axis label formatting
    axis.text = element_text(size = 12),                              # Axis text size
    legend.position = "right"                                        # Position the legend
    #plot.margin = unit(c(1,1,1,1), "lines")                           # Set margin to make it more square
  ) +
  scale_color_manual(values = custom_colors)   # Adjust aspect ratio for better spacing on x-axis


# Required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(purrr)
setwd("/home/raafat/Documents/Sugar_Beet")
# Load your data
data <- read.csv("Data.csv")
data$T.S.S.
# Choose the trait to plot
selected_trait <- "T.S.S."  # Replace with any other column name if needed

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
    model <- aov(Value ~ Treatments+ Varieties, data = .x)
    lsd <- LSD.test(model, "Treatments", group = TRUE, alpha = 0.05)$groups
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
  labs(title = "(A)", x = "Treatments", y = "TSS%") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Title centered and bold
    axis.title.x = element_text(size = 14, face = "bold"),            # X-axis label formatting
    axis.title.y = element_text(size = 14, face = "bold"),            # Y-axis label formatting
    axis.text = element_text(size = 12),                              # Axis text size
    legend.position = "right"                                         # Position the legend
  ) +
  scale_color_manual(values = custom_colors)

