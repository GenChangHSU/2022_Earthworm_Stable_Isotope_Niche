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


# Bayesian SEA
parms <- list()
parms$n.iter <- 2*10^4     # Number of iterations per chain
parms$n.burnin <- 1*10^3   # Burn-in
parms$n.thin <- 10         # Thinning
parms$n.chains <- 3        # Number of chains

priors <- list()
priors$R <- 1*diag(2)      # Inverse Wishart prior
priors$k <- 2
priors$tau.mu <- 0.001

SEA_posterior1 <- siberMVN(SIBER_data1, parms, priors) %>% 
  siberEllipses() %>%
  as.data.frame() %>% 
  `names<-`(SIBER_data1$all.groups)

SEA_Bayes1 <- rbind(
  apply(SEA_posterior1, 2, mean),
  apply(SEA_posterior1, 2, function(x) {hdr(x, prob = c(95))$hdr[1, 1]}),
  apply(SEA_posterior1, 2, function(x) {hdr(x, prob = c(95))$hdr[1, 2]})) %>%
  as.data.frame() %>%
  `row.names<-`(c("SEAb", "SEAb_HDR_lower", "SEAb_HDR_upper"))

SEA_posterior2 <- siberMVN(SIBER_data2, parms, priors) %>% 
  siberEllipses() %>%
  as.data.frame() %>% 
  `names<-`(SIBER_data1$all.groups)

SEA_Bayes2 <- rbind(
  apply(SEA_posterior2, 2, mean),
  apply(SEA_posterior2, 2, function(x) {hdr(x, prob = c(95))$hdr[1, 1]}),
  apply(SEA_posterior2, 2, function(x) {hdr(x, prob = c(95))$hdr[1, 2]})) %>%
  as.data.frame() %>%
  `row.names<-`(c("SEAb", "SEAb_HDR_lower", "SEAb_HDR_upper"))

# Plot the SEAb
SEA_Bayes_plot_df1 <- SEA_Bayes1 %>%
  t() %>% 
  as.data.frame() %>%
  mutate(Species = row.names(.)) %>%
  select(Species, 1:3)

ggplot(data = SEA_Bayes_plot_df1, aes(x = fct_rev(Species), y = SEAb, color = Species)) + 
  geom_pointrange(aes(ymin = SEAb_HDR_lower, ymax = SEAb_HDR_upper)) +
  coord_flip() + 
  labs(x = NULL, y = "Bayesian SEA (\u2030) (Mean ± 95% HDR)") +
  scale_x_discrete(label = rev(c(expression(italic("Allolobophora chlorotica")),
                                 expression(italic("Aporrectodea caliginosa")),
                                 expression(italic("Aporrectodea trapezoides")),
                                 expression(italic("Lumbricus friendi")),
                                 expression(italic("Lumbricus rubellus"))))) + 
  scale_y_continuous(limits = c(-1, 31), breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_color_npg(guide = F) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
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
        panel.grid.minor.y = element_blank())

ggsave("Output/Figures/SEAb1.tiff", width = 6, height = 5, dpi = 600)



# ### 2. Population-level Layman metrics
# SIBER_data2 <- Earthworm_data_corrected %>% 
#   select(iso1 = d13C_earthworm_correct, 
#          iso2 = d15N_earthworm_correct,
#          community = Species) %>%
#   arrange(community) %>%
#   dplyr::group_by(community) %>%
#   dplyr::mutate(group = 1:n()) %>%
#   select(1, 2, 4, 3) %>%
#   as.data.frame() %>%
#   createSiberObject()
# 
# # ML estimates
# Layman_metrics_ML <- communityMetricsML(SIBER_data2) %>% 
#   as.data.frame() %>%
#   `row.names<-`(c("NR", "CR", row.names(.)[3:6]))
# 
# Layman_ML_plot_df <- Layman_metrics_ML %>%
#   t() %>% 
#   as.data.frame() %>%
#   mutate(Species = row.names(.)) %>%
#   pivot_longer(cols = -Species, 
#                names_to = "Layman_mtcs", 
#                values_to = "Value") %>%
#   mutate(Layman_mtcs = fct_relevel(Layman_mtcs, "NR", "CR", "TA", "CD", "MNND", "SDNND"))
# 
# ggplot(data = Layman_ML_plot_df, 
#        aes(x = Species, y = Value, fill = Species)) + 
#   geom_bar(stat = "identity") +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
#   facet_wrap(~Layman_mtcs, ncol = 3, scales = "free") + 
#   labs(x = NULL, y = expression(paste("Value", " (\u2030)", sep = ""))) +
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 12, color = "black"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 15, margin = margin(r = 8)),
#         axis.ticks.x = element_blank(),
#         plot.title = element_text(hjust = 0.5, size = 18),
#         plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(colour = "transparent"),
#         panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         strip.background = element_rect(fill = "white"),
#         strip.text = element_text(size = 12),
#         legend.position = c(0.5, 1.08),
#         legend.direction = "horizontal",
#         legend.spacing.x = unit(0.1, "cm"),
#         legend.spacing.y = unit(0.15, "cm"),
#         legend.key.width = unit(0.5, "cm"),
#         legend.key.size = unit(1, "line"),
#         legend.key = element_blank(),
#         legend.text = element_text(size = 8),
#         legend.text.align = 0,
#         legend.box.just = "center",
#         legend.justification = c(0.5, 0.5),
#         legend.title.align = 0.5,
#         legend.background = element_rect(fill = "transparent", 
#                                          size = 0.5, 
#                                          linetype = "solid", 
#                                          colour = "transparent")) +
#   scale_fill_npg(name = NULL, 
#                  label = c(expression(italic("Allolobophora chlorotica    ")),
#                            expression(italic("Aporrectodea caliginosa    ")),
#                            expression(italic("Aporrectodea trapezoides    ")),
#                            expression(paste(italic("Diplocardia"), " sp.    ")),
#                            expression(italic("Lumbricus friendi    ")),
#                            expression(italic("Lumbricus rubellus    ")))) + 
#   guides(fill = guide_legend(nrow = 1))
# 
# ggsave("Output/Figures/Layman.tiff", width = 9, height = 6, dpi = 600)
# 
