setwd("/home/raafat/Documents/Sugar_Beet")


data_beat <- read.csv("Data.csv", header=T)


# Varieties * Treatmemts
# Load necessary libraries
library(ggplot2)
library(factoextra)  # For PCA visualization
library(dplyr)       # For data manipulation
library(FactoMineR)
library("factoextra")
data_beat$Varieties <- as.factor(data_beat$Varieties)
data_beat$Treatments <- as.factor(data_beat$Treatments)

factor_columns <- c("Varieties", "Treatments")
trait_data <- data_beat %>% select(-one_of(factor_columns))
any(is.na(trait_data))       # Checks for NA values
any(is.infinite(trait_data)) # Checks for infinite values

pca_result <- prcomp(trait_data, scale. = TRUE)
summary(pca_result)
pca_scores <- as.data.frame(pca_result$x)

# Combine PCA scores with genotype and treatment for plotting
pca_scores$Varieties <- data_beat$Varieties
pca_scores$Treatments <- data_beat$Treatments

# Create a PCA biplot using ggplot2 and color by Genotype and shape by Treatment
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Varieties, shape = Treatments)) +
  geom_point(size = 4) +  # Scatter plot with points
  theme_minimal() +       # Clean theme
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add dashed lines for the origin
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right") +  # Legend position
  scale_color_brewer(palette = "Set1")  # Color palette for Genotypes


?stat_ellipse
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Varieties, shape = Treatments)) +
  geom_point(size = 4) +  # Scatter plot with points
  stat_ellipse(aes(group = Varieties), type = "t", level = 0.95, linetype = "solid") +  # Ellipses around each group
  theme_minimal() +       # Clean theme
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add dashed lines for the origin
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right") +  # Legend position
  scale_color_brewer(palette = "Set1")  # Color palette for Varieties

fviz_pca_var(pca_result, col.var = "black", repel = TRUE)  # Arrows showing variable contributions
fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_ind(pca_result,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = data_beat$Varieties, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#64e022", "#9e0be2", "#f13e3e"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
# Add confidence ellipses
fviz_pca_ind(pca_result, geom.ind = "point", col.ind = data_beat$Varieties, 
             palette = c("#00AFBB", "#E7B800", "#FC4E07", "#64e022", "#9e0be2", "#f13e3e"),
             addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Groups"
)



fviz_pca_biplot(pca_result, 
                col.ind = data_beat$Varieties, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species") 




fviz_pca_biplot(pca_result, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                fill.ind = data_beat$Varieties,
                col.ind = data_beat$Treatments,
                # Color variable by groups
               # col.var = factor(c("sepal", "sepal", "petal", "petal")),
               addEllipses = TRUE,
                
                legend.title = list(fill = "Species", color = "Clusters"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")      # Variable colors



# Load ggExtra if not already loaded
library(ggExtra)

# Create the base PCA plot
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Varieties, shape = Treatments)) +
  geom_point(size = 4) +  # Scatter plot with points
  stat_ellipse(aes(group = Varieties), type = "norm", level = 0.95, linetype = "solid") +  # Ellipses around each group
  theme_minimal() +       # Clean theme
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add dashed lines for the origin
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "right") +  # Legend position
  scale_color_brewer(palette = "Set1")  # Color palette for Varieties

# Add marginal density plots
ggMarginal(p, type = "density", margins = "both", groupColour = TRUE, groupFill = TRUE)

