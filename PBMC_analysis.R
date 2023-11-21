library(data.table)
library(BiocParallel)
library(scDblFinder)
library(SingleR)
library(Seurat)
library(future)
library(scran)
library(scater)
library(ggrepel)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(ggsci)
library(ggbeeswarm)
library(ggpubr)
library(dittoSeq)
library(Matrix.utils)
library(RColorBrewer)
library(ComplexHeatmap)
library(scales)

sc.sce.rna <- readRDS("FD1-6_integrated_annot.rds")
sc.sce.rna$Timepoint <- factor(sc.sce.rna$Condition, levels = c("BL", "FU"), labels = c("0", "1"))
sc.sce.rna$Timepoint2 <- as.numeric(sc.sce.rna$Timepoint) - 1
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

### Figure 4
lab_main <- sc.t1[, .(x = median(UMAP.1),
                      y = median(UMAP.2)),
                  by = Monaco_cluster_main]


umap.main <-ggplot(sc.t1, aes(x = UMAP.1, y = UMAP.2, color = Monaco_cluster_main)) +
  geom_point(size = 0.4, alpha = 0.8) +
  geom_label_repel(data = lab_main, aes(x = x, y = y, label = Monaco_cluster_main), color = "black", size = 6) +
  scale_color_npg(name = "Cell type") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

dotplot.cyto <- dittoDotPlot(sc.sce.rna, 
                             vars = c("CXCL9", "CXCL10", "IFNG", "IL10", "IL12A", "IL21"),
                             group.by = "Monaco_cluster_main",
                             split.by = c("AE", "Timepoint2"),
                             split.adjust = list(labeller = label_bquote(cols = TP[.(Timepoint2)])),
                             
                             max.color = "#EE0000FF",
                             scale = FALSE,
                             ylab = NULL,
                             theme =
                               theme(axis.text = element_text(size = 16),
                                     panel.background = element_rect(linewidth = 1, color = "black", fill = "white"),
                                     strip.background = element_rect(linewidth = 1, color = "black", fill = "white"),
                                     #text = element_text(family = "Arial"),
                                     strip.text = element_text(size = 14)
                               ))

umap.ifng <- ggplot(sc.t1[order(IFNG),], aes(x = UMAP.1, y = UMAP.2, color = IFNG)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint2, labeller = label_bquote(cols = TP[.(Timepoint2)])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.ifng <- ggplot(sc.t1, aes(x = Condition, y = IFNG, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 3, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
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
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")

umap.cxcl9 <- ggplot(sc.t1[order(CXCL9),], aes(x = UMAP.1, y = UMAP.2, color = CXCL9)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint2, labeller = label_bquote(cols = TP[.(Timepoint2)])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.cxcl9 <- ggplot(sc.t1, aes(x = Condition, y = CXCL9, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 3, groupOnX = TRUE) +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
  scale_color_d3(name = "Cell type") +
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
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1])))  +
  xlab("")

umap.cxcl10 <- ggplot(sc.t1[order(CXCL10),], aes(x = UMAP.1, y = UMAP.2, color = CXCL10)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint2, labeller = label_bquote(cols = TP[.(Timepoint2)])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom") 

vio.cxcl10 <- ggplot(sc.t1, aes(x = Condition, y = CXCL10, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 5, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
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
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")

umap.il10 <- ggplot(sc.t1[order(IL10),], aes(x = UMAP.1, y = UMAP.2, color = IL10)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint2, labeller = label_bquote(cols = TP[.(Timepoint2)])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.il10 <- ggplot(sc.t1, aes(x = Condition, y = IL10, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 3, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
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
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")

umap.il12a <- ggplot(sc.t1[order(IL12A),], aes(x = UMAP.1, y = UMAP.2, color = IL12A)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint2, labeller = label_bquote(cols = TP[.(Timepoint2)])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.il12a <- ggplot(sc.t1, aes(x = Condition, y = IL12A, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 3, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
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
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")

umap.il21 <- ggplot(sc.t1[order(IL21),], aes(x = UMAP.1, y = UMAP.2, color = IL21)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis_c(option = "turbo") +
  theme_bw() +
  facet_grid(AE~Timepoint2, labeller = label_bquote(cols = TP[.(Timepoint2)])) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18),
        legend.position = "bottom")

vio.il21 <- ggplot(sc.t1, aes(x = Condition, y = IL21, color = Monaco_cluster_main)) +
  geom_violin(size = 0.6, alpha = 0.8) +
  geom_quasirandom(bandwidth = 3, groupOnX = TRUE) +
  scale_color_d3(name = "Cell type") +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
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
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")

fig4 <- (umap.main | dotplot.cyto) / (umap.cxcl9 | vio.cxcl9 | umap.cxcl10 | vio.cxcl10) / (umap.ifng | vio.ifng | umap.il10 | vio.il10) / (umap.il12a | vio.il12a | umap.il21 | vio.il21) + plot_annotation(title = "Figure 4", tag_levels = "A") &  theme(plot.title = element_text(size = 24), plot.tag = element_text(size = 20))
ggsave("Figure4.pdf", plot = fig4, width = 20, height = 28)
ggsave("Figure4.png", plot = fig4, width = 20, height = 28)

#Figure 5
lab_coord <- sc.t1[Monaco_cluster_fine2 %in% c("Terminal effector CD8 T cells", "Naive CD8 T cells", "Th17 cells", 
                                               "MAIT cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", 
                                               "Naive CD4 T cells", "Vd2 gd T cells", "T regulatory cells", 
                                               "Non-Vd2 gd T cells", "Th1 cells"), 
                   .(x = median(UMAP.1),
                     y = median(UMAP.2)),
                   by = Monaco_cluster_fine2]

umap.tcell <- ggplot(sc.t1[Monaco_cluster_fine2 %in% c("Terminal effector CD8 T cells", "Naive CD8 T cells", "Th17 cells", 
                                                       "MAIT cells", "Central memory CD8 T cells", "Effector memory CD8 T cells", 
                                                       "Naive CD4 T cells", "Vd2 gd T cells", "T regulatory cells", 
                                                       "Non-Vd2 gd T cells", "Th1 cells")], 
                     aes(x = UMAP.1, y = UMAP.2, color = Monaco_cluster_fine2)) +
  geom_point(size = 0.4, alpha = 0.8) +
  geom_label_repel(data = lab_coord, aes(x = x, y = y, label = Monaco_cluster_fine2), color = "black", size = 6) +
  scale_color_d3("category20", name = "T Cell subtype") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18)) +
  guides(color = guide_legend(override.aes = list(size = 3),ncol = 1)) +
  coord_cartesian(xlim = c(-1, 13), ylim = c(-7, 4))

fine.ht <- lapply(unique(sc.t1[Monaco_cluster_main %in% c("CD8+ T cells", "CD4+ T cells"), Monaco_cluster_fine]), function(x) {
  sc.test <- sc.sce.rna[,sc.sce.rna$Monaco_cluster_fine %in% x]
  
  mark.test <- findMarkers(sc.test, group = sc.test$condition_ae)
  
  #select cytokines going up from BL to FU in AE
  mark.test.1 <- as.data.frame(mark.test[[1]])
  mark.test.1a <- rownames(mark.test.1[mark.test.1$logFC.FU.AE < -0.25 & mark.test.1$FDR < 0.05,])
  
  #select cytokines going up from BL to FU in no AE
  mark.test.2 <- as.data.frame(mark.test[[2]])
  mark.test.2a <- rownames(mark.test.2[mark.test.2$logFC.FU.no.AE < -0.25 & mark.test.2$FDR < 0.05,])
  
  ae.genes <- setdiff(mark.test.1a, mark.test.2a)
  noae.genes <- setdiff(mark.test.2a, mark.test.1a)
  mark.genes <- c(ae.genes[!grepl("[MR][TP][SL-]", ae.genes)], noae.genes[!grepl("[MR][TP][SL-]", noae.genes)])
  x2 <- gsub(" ", "_", x)
  #plotGroupedHeatmap(sc.test, features = mark.genes, group = "condition_ae", main = paste(x), center = TRUE, cluster_cols = FALSE,
  #                   filename = paste0("230202_analysis/231010_Heatmap_DEG_", x2, ".pdf"))
  
  sc.test.ht <- assay(summarizeAssayByGroup(sc.test, sc.test$condition_ae, statistics = "mean"))
  sc.test.ht2 <- sc.test.ht[mark.genes,]
  sc.test.ht3 <- apply(sc.test.ht2, 1, scale)
  rownames(sc.test.ht3) <- colnames(sc.test.ht2)
  
  ca <- columnAnnotation(foo = anno_mark(at = c(grep("STAT1", colnames(sc.test.ht3)),
                                                grep("STAT3", colnames(sc.test.ht3))),
                                         labels = c("STAT1", 
                                                    "STAT3"),
                                         which = "column", 
                                         side = "bottom"))
  
  x3 <- ifelse(x == "Terminal effector CD8 T cells", "Terminal effector\nCD8  T cells", 
               ifelse(x == "Naive CD8 T cells", "Naive CD8\n T cells", 
                      ifelse(x == "Central memory CD8 T cells", "Central memory\nCD8 T cells",
                             ifelse(x == "Th17 cells", "Th17\n cells", 
                                    ifelse(x == "Effector memory CD8 T cells", "Effector memory\nCD8 T cells", 
                                           ifelse(x == "Naive CD4 T cells", "Naive CD4\nT cells", 
                                                  ifelse(x == "T regulatory cells", "T regulatory\n cells", x)))))))
  if(x == "Terminal effector CD8 T cells") {
    Heatmap(sc.test.ht3,
            name = "expression",
            column_title = x3,
            show_column_names = FALSE,
            bottom_annotation = ca, 
            row_split = c("AE", "no AE", "AE", "no AE"),
            row_title_rot = 0,
            row_title_gp = gpar(fontsize = 20),
            row_labels = expression(TP[0],TP[0],TP[1],TP[1]),
            row_names_rot = 0,
            row_names_gp = gpar(fontsize = 20))
  } else {
    Heatmap(sc.test.ht3,
            column_title = x3,
            show_column_names = FALSE,
            bottom_annotation = ca, 
            
            row_split = c("AE", "no AE", "AE", "no AE"),
            row_title_rot = 0,
            row_title_gp = gpar(fontsize = 20),
            row_labels = expression(TP[0],TP[0],TP[1],TP[1]),
            row_names_rot = 0,
            row_names_gp = gpar(fontsize = 20),
            show_heatmap_legend = FALSE)
  }
  
})

#heatmap in row format
fine.ht2 <- fine.ht[[1]] + fine.ht[[2]] + fine.ht[[3]] + fine.ht[[4]] + fine.ht[[5]] + fine.ht[[6]] + fine.ht[[7]] + fine.ht[[8]]

fine.ht3 <- grid.grabExpr(draw(fine.ht2, ht_gap = unit(0.5, "cm")))

vio.stat1 <- ggplot(sc.t1[Monaco_cluster_main %in% c("CD8+ T cells", "CD4+ T cells")], 
                    aes(x = Condition, y = STAT1, color = Monaco_cluster_fine2)) +
  geom_quasirandom() +
  scale_color_igv(name = "T cell subtype") +
  facet_grid(Monaco_cluster_fine2~AE, labeller = labeller(Monaco_cluster_fine2 = label_wrap_gen(width = 20))) +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
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
  ylab("STAT1 (normalized counts)") +
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")
ggsave("230202_analysis/231010_Violin_Tcells_STAT1.pdf", plot = vio.stat1, width = 12, height = 8)

vio.stat3 <- ggplot(sc.t1[Monaco_cluster_main %in% c("CD8+ T cells", "CD4+ T cells")], 
                    aes(x = Condition, y = STAT3, color = Monaco_cluster_fine2)) +
  geom_quasirandom() +
  scale_color_igv() +
  facet_grid(Monaco_cluster_fine2~AE, labeller = labeller(Monaco_cluster_fine2 = label_wrap_gen(width = 20))) +
  stat_compare_means(comparisons = list(c("BL", "FU"))) +
  
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
  ylab("STAT3 (normalized counts)") +
  scale_x_discrete(labels = c("BL" = bquote(TP[0]), "FU" = bquote(TP[1]))) +
  xlab("")
ggsave("230202_analysis/230303_Violin_Tcells_STAT3.pdf", plot = vio.stat3, width = 12, height = 8)

fig5 <- (umap.tcell | plot_spacer()) / fine.ht3 / (vio.stat1 | vio.stat3) + plot_annotation(title = "Figure 5", tag_levels = "A") &  theme(plot.title = element_text(size = 24), plot.tag = element_text(size = 20))
ggsave("Figure5.pdf", plot = fig5, width = 24, height = 20)
ggsave("Figure5.png", plot = fig5, width = 24, height = 20)

#Supplemental Figure 4

umap.fine.ae <- ggplot(sc.t1, aes(x = UMAP.1, y = UMAP.2, color = Monaco_cluster_fine2)) +
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
  guides(color = guide_legend(override.aes = list(size = 3),ncol = 2)) +
  facet_grid(AE ~ Timepoint2,
             labeller = label_bquote(cols = TP[.(Timepoint2)]))
prop.fine <- ggplot(sc.t1, aes(x = Timepoint, fill = Monaco_cluster_fine2)) +
  geom_bar(position = "fill") +
  scale_fill_igv(name = "Cell type") +
  facet_wrap(~AE) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18), 
        strip.text = element_text(size = 18)) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_x_discrete(labels = c(bquote(TP[0]), bquote(TP[1]))) +
  ylab("Proportion")

sfig4 <- umap.fine.ae / prop.fine + plot_layout(guides = "collect") + plot_annotation(title = "Supplemental Figure 4", tag_levels = "A") &  theme(plot.title = element_text(size = 24), plot.tag = element_text(size = 20))
ggsave("Supplemental_Figure4.pdf", plot = sfig4, width = 12, height = 8)

