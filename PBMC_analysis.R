library(scran)
library(scater)
library(ggrepel)
library(ggplot2)
library(ggsci)
library(ggbeeswarm)
library(ggpubr)
library(patchwork)

sc.sce.rna <- readRDS("SCE_Flo_FD1-6_integrated_annot.rds")

scater::plotDots(sc.sce.rna, 
        features = c("CD79A", "CD3E", "CD8A", "GZMB", "KLRD1",  "CD68", "CD86","CD36"  ), 
        group = "Monaco_cluster_main")

sc.t1 <- makePerCellDF(sc.sce.rna, features = c("JUN", "JUNB", "FOS", "CD163", "S100A8", "HLA-A", "B2M", 
                                                "LTB", "FTL" ,"IL17A", "IFNG", "IL10", "IL6", 
                                                "IL23A", "IL22", "IL12A", "CCL3", "CCL4", "GZMB",
                                                "CXCR3", "STAT3", "CSF2", "STAT4", "CXCL10", "CXCL9",
                                                "CCR6", "CCR4", "TNF", "IL26" , "IL17RA", "IL17RB", 
                                                "IL17RC", "IRF4", "RUNX1", "BATF", "CCL20", 
                                                "FOXP3", "IL2RA", "IL17A", "IL17F", "RORC",
                                                "STAT1", "FOS", "FKBP5", "STAT1", "STAT3", "LTA",
                                                "IL12A", "IL12B", "IL21"))
setDT(sc.t1)
sc.t1[, UMAP.1 := reducedDim(sc.sce.rna, "UMAP")[,1]]

sc.t1[, UMAP.2 := reducedDim(sc.sce.rna, "UMAP")[,2]]

sc.t1[, Monaco_cluster_fine2 := Monaco_cluster_fine]
sc.t1[is.na(Monaco_cluster_fine), Monaco_cluster_fine2 := "Unknown"]

sc.t1[Condition == "BL", Timepoint := "TP\U02080"]
sc.t1[Condition == "FU", Timepoint := "TP\U02081"]


ggplot(sc.t1, aes(x = UMAP.1, y = UMAP.2, color = Monaco_cluster_main)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_d3(name = "Cell type") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave("Supplementary_Figure_6A.pdf", width = 12, height = 8, device = cairo_pdf)

ggplot(sc.t1, aes(x = UMAP.1, y = UMAP.2, color = Monaco_cluster_fine2)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_igv(name = "Cell subtype") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3),ncol = 1)) +
  facet_grid(AE ~ Timepoint)
ggsave("Supplementary_Figure_6B.pdf", width = 12, height = 8, device = cairo_pdf)

ggplot(sc.t1, aes(x = Timepoint, fill = Monaco_cluster_fine2)) +
  geom_bar(position = "fill") +
  scale_fill_igv(name = "Cell type") +
  facet_wrap(~AE) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18)) +
  guides(fill = guide_legend(ncol = 1))

ggsave("Supplementary_Figure_6C.pdf", width = 12, height = 8, device = cairo_pdf)

sc.prop <- sc.t1[, .N, by = c("AE", "Timepoint", "Monaco_cluster_fine2")]
sc.prop[, Prop := round(N/sum(N) * 100, 1), by = c("AE", "Timepoint")]
sc.prop[, N_other := sum(N) - N, by = c("AE", "Timepoint")]

sc.prop[AE %in% c("no AE", "AE") & Timepoint == "TP\U2080", sample_class := "1"]

sc.prop[AE %in% c("AE") & Timepoint == "TP\U2081", sample_class := "2"]
sc.prop[AE %in% c("no AE") & Timepoint == "TP\U2081", sample_class := "3"]
sc.prop[Monaco_cluster_fine2 == "Naive CD4 T cells", Monaco_cluster_fine2 := "Naive CD4+ T cells"]
sc.prop[Monaco_cluster_fine2 == "Naive CD8 T cells", Monaco_cluster_fine2 := "Naive CD8+ T cells"]

chisq.test(sc.prop[Monaco_cluster_fine2 == "Intermediate monocytes" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Intermediate monocytes" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Myeloid dendritic cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Myeloid dendritic cells" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Naive CD4+ T cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Naive CD4+ T cells" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Naive CD8+ T cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Naive CD8+ T cells" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Plasmacytoid dendritic cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Plasmacytoid dendritic cells" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Th17 cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Th17 cells" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Non-switched memory B cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Non-switched memory B cells" & AE == "no AE", .(N, N_other)])

chisq.test(sc.prop[Monaco_cluster_fine2 == "Th1 cells" & AE == "AE", .(N, N_other)])
chisq.test(sc.prop[Monaco_cluster_fine2 == "Th1 cells" & AE == "no AE", .(N, N_other)])

sigdf <- data.table(Monaco_cluster_fine2 = c("Intermediate monocytes", "Intermediate monocytes", 
                                             "Naive CD4+ T cells","Naive CD4+ T cells",
                                             "Naive CD8+ T cells","Naive CD8+ T cells",
                                             "Th17 cells", "Th17 cells",
                                             "Myeloid dendritic cells", "Myeloid dendritic cells",
                                             "Th1 cells", "Th1 cells"),
                    AE = c("AE", "no AE", 
                           "AE", "no AE", 
                           "AE", "no AE",
                           "AE", "no AE",
                           "AE", "no AE",
                           "AE", "no AE"),
                    X1 = c("TP\U02080"),
                    Timepoint = c("TP\U02081"),
                    sample_class = "2",
                    y_pos = c(14, 6,
                              9,5,
                              4,2,
                              9,8,
                              3,3,
                              3,5),
                    value = c("2.2e-16", "0.000594",
                              "2.2e-16", "0.345",
                              "2.22e-16", "0.000392",
                              "2.2e-16", "5.22e-16",
                              "2.2e-16", "0.00887",
                              "1.62e-7", "2.2e-16"))
ggplot(sc.prop[Monaco_cluster_fine2 %in% c("Th17 cells", 
                                           "Intermediate monocytes", 
                                           "Naive CD4+ T cells",
                                           "Naive CD8+ T cells",
                                           "Myeloid dendritic cells",
                                           "Th1 cells")], 
       aes(x = Timepoint, y = Prop, fill = sample_class)) +
  geom_bar(stat = "identity") +
  ggsignif::geom_signif(data = sigdf,
                        aes(annotations = value, xmin = X1, xmax = Timepoint, y_position = y_pos), 
                        manual = TRUE) +
  scale_fill_npg() +
  facet_grid(Monaco_cluster_fine2 ~ AE) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none",
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16), 
        strip.text = element_text(size = 16),
        strip.text.y = element_text(angle = 0)) +
  ylab("Proportion (%)") +
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))
ggsave("Supplementary_Figure8.pdf", width = 12, height = 8, device = cairo_pdf)





vio.stat1 <- ggplot(sc.t1[Monaco_cluster_main %in% c("CD8+ T cells", "CD4+ T cells")], 
                    aes(x = Timepoint, y = STAT1, color = Monaco_cluster_fine2)) +
  geom_quasirandom() +
  scale_color_igv(name = "T cell subtype") +
  facet_grid(Monaco_cluster_fine2~AE, labeller = labeller(Monaco_cluster_fine2 = label_wrap_gen(width = 20))) +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16), 
        legend.position = "none",
        strip.text = element_text(size = 16),
        strip.text.y = element_text(angle = 0)) + 
  guides(color = guide_legend(override.aes = list(size = 3))) +
  
  scale_y_continuous(expand = expansion(mult = c(0,0.3))) +
  ylab("STAT1 (normalized counts)")
ggsave("Figure_4A.pdf", plot = vio.stat1, width = 12, height = 8, device = cairo_pdf)

vio.stat3 <- ggplot(sc.t1[Monaco_cluster_main %in% c("CD8+ T cells", "CD4+ T cells")], 
                    aes(x = Timepoint, y = STAT3, color = Monaco_cluster_fine2)) +
  geom_quasirandom() +
  scale_color_igv() +
  facet_grid(Monaco_cluster_fine2~AE, labeller = labeller(Monaco_cluster_fine2 = label_wrap_gen(width = 20))) +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16), 
        legend.position = "none",
        strip.text = element_text(size = 16),
        strip.text.y = element_text(angle = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.3))) +
  ylab("STAT3 (normalized counts)")
ggsave("Figure_4B.pdf", plot = vio.stat3, width = 12, height = 8, device = cairo_pdf)


umap.cxcl9 <- ggplot(sc.t1[order(CXCL9),], aes(x = UMAP.1, y = UMAP.2, color = CXCL9)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.cxcl9 <- ggplot(sc.t1, aes(x = Timepoint, y = CXCL9, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  facet_grid(Monaco_cluster_main~AE,labeller = labeller(Monaco_cluster_main = label_wrap_gen(width = 10))) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,2,4), expand = expansion(mult = c(0,0.3)))  +
  ylab("CXCL9 (normalized counts)")

umap.cxcl9 + vio.cxcl9
ggsave("Figure_4C.pdf", width = 12, height = 8, device = cairo_pdf)


umap.cxcl10 <- ggplot(sc.t1[order(CXCL10),], aes(x = UMAP.1, y = UMAP.2, color = CXCL10)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.cxcl10 <- ggplot(sc.t1, aes(x = Timepoint, y = CXCL10, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  facet_grid(Monaco_cluster_main~AE,labeller = labeller(Monaco_cluster_main = label_wrap_gen(width = 10))) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,2,4), expand = expansion(mult = c(0,0.3))) +
  ylab("CXCL10 (normalized counts)")

umap.cxcl10 + vio.cxcl10
ggsave("Figure_4D.pdf", width = 12, height = 8, device = cairo_pdf)

umap.ifng <- ggplot(sc.t1[order(IFNG),], aes(x = UMAP.1, y = UMAP.2, color = IFNG)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.ifng <- ggplot(sc.t1, aes(x = Timepoint, y = IFNG, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  facet_grid(Monaco_cluster_main~AE,labeller = labeller(Monaco_cluster_main = label_wrap_gen(width = 10))) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,2,4), expand = expansion(mult = c(0,0.3)))  +
  ylab("IFNG (normalized counts)")

umap.ifng + vio.ifng
ggsave("Figure_4E.pdf", width = 12, height = 8, device = cairo_pdf)

umap.il10 <- ggplot(sc.t1[order(IL10),], aes(x = UMAP.1, y = UMAP.2, color = IL10)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.il10 <- ggplot(sc.t1, aes(x = Timepoint, y = IL10, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  facet_grid(Monaco_cluster_main~AE,labeller = labeller(Monaco_cluster_main = label_wrap_gen(width = 10))) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,2,4), expand = expansion(mult = c(0,0.3)))  +
  ylab("IL10 (normalized counts)")

umap.il10 + vio.il10
ggsave("Figure_4F.pdf", width = 12, height = 8, device = cairo_pdf)

umap.tnf <- ggplot(sc.t1[order(TNF),], aes(x = UMAP.1, y = UMAP.2, color = TNF)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.tnf <- ggplot(sc.t1, aes(x = Timepoint, y = TNF, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  facet_grid(Monaco_cluster_main~AE,labeller = labeller(Monaco_cluster_main = label_wrap_gen(width = 10))) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,2,4), expand = expansion(mult = c(0,0.3)))  +
  ylab("TNF (normalized counts)")

umap.tnf + vio.tnf
ggsave("Figure_4G.pdf", width = 12, height = 8, device = cairo_pdf)

umap.lta <- ggplot(sc.t1[order(LTA),], aes(x = UMAP.1, y = UMAP.2, color = LTA)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.lta <- ggplot(sc.t1, aes(x = Timepoint, y = LTA, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("TP\U02080", "TP\U02081"))) +
  
  theme_bw() +
  facet_grid(Monaco_cluster_main~AE,labeller = labeller(Monaco_cluster_main = label_wrap_gen(width = 10))) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        strip.text.y = element_text(angle = 0),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,2,4), expand = expansion(mult = c(0,0.3)))  +
  ylab("LTA (normalized counts)")

umap.lta + vio.lta
ggsave("Figure_4H.pdf", width = 12, height = 8, device = cairo_pdf)
