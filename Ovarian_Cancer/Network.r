pathway_colors <- c(
  "mTOR" = "#D55E00",
  "Testosterone" = "#0072B2",
  "Androgen" = "#009E73",
  "Other" = "grey80"
)

V(G)$color <- pathway_colors[V(G)$pathway]

# Scale node size by gene membership
ME <- MEs[, paste0("ME", module_of_interest)]
kME <- cor(datExpr[, Co-exp_m_t_e], ME, use = "p")

V(G)$size <- scales::rescale(abs(kME), to = c(6, 16))

plot(
  G,
  layout = layout_with_fr,
  vertex.label = V(G)$label,
  vertex.label.cex = 0.75,
  vertex.label.color = "black",
  vertex.frame.color = NA,
  edge.width = E(G)$weight * 3,
  edge.color = "grey60",
  main = "Co-expression Network of mTOR, Testosterone and Androgen Pathways"
)
