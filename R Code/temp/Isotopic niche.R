# Plot the percent overlap between species pairs
labels <- c(expression(italic("Allolobophora \n   chlorotica")),
            expression(italic("Aporrectodea \n   caliginosa")),
            expression(italic("Aporrectodea \n  trapezoides")),
            expression(italic("Lumbricus \n   friendi")),
            expression(italic("Lumbricus \n  rubellus")))

ggplot(Isoniche_stats_df1, aes(x = Sp1, y = Sp2)) + 
  geom_point(aes(size = Percent_overlap_Sp2*100, color = as.factor(PERMANOVA_Pval_sig))) + 
  geom_abline(intercept = 0, slope = 1, size = 2, color = "grey") +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1, vjust = 0.75),
        axis.text.y = element_text(size = 10, color = "black", hjust = 1, vjust = 0.75),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.spacing.x = unit(-1, "cm"),
        legend.spacing.y = unit(0.15, "cm"),
        legend.key.width = unit(3.5, "cm"),
        legend.key.size = unit(1, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.text.align = 0.5,
        legend.box.just = "right",
        legend.justification = c(0.5, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) + 
  scale_x_discrete(labels = labels) + 
  scale_y_discrete(labels = labels) + 
  scale_size(range = c(2, 15), 
             name = "Percent of Sp2's niche \n overlapping with Sp1 (%)",
             breaks = c(5, 20, 50, 75, 100)) + 
  scale_color_manual(values = c("black", "red"), guide = F) +
  coord_fixed(ratio = 1)
  
ggsave("Output/Figures/Niche_overlap1.tiff", width = 8, height = 6, dpi = 600)


