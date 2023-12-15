#Code for Figure 6 and 7
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ComplexHeatmap)
library(data.table)
library(ggbeeswarm)
library(ggpubr)
library(dittoSeq)
library(ggsci)

# Normal skin 2  Skin rash 6 
img.skin <- readPNG("skin_hires_image.png")
# colitis 1 colon 1 
img.col <- readPNG("colon_hires_image.png")

cell_cols <- pal_d3("category20")(11)
names(cell_cols) <- c("Fibroblasts and immune cells",
                      "Epithelial and immune cells",
                      "Epithelial cell",
                      "Eccrine cell",
                      "B cell",
                      "Fibroblast",
                      "Apoptotic cell",
                      "T and B cell",
                      "Muscle cell",
                      "Keratinocyte",
                      "Hair follicle cell")

s1b <- readRDS("AE_Visium_Seurat_harmony.rds")
rownames(s1b@meta.data) <- s1b$Cell
s1b <- RenameCells(s1b, new.names = s1b$Cell)

s1b.meta <- s1b@meta.data
setDT(s1b.meta)
s1b.meta[, UMAP.1 := s1b@reductions$umap@cell.embeddings[,1]]
s1b.meta[, UMAP.2 := s1b@reductions$umap@cell.embeddings[,2]]

s1.0.coord <- data.table(Cell = rownames(s1b@images$slice1@coordinates), s1b@images$slice1@coordinates)
s1.1.coord <- data.table(Cell = rownames(s1b@images$slice1.1@coordinates), s1b@images$slice1.1@coordinates)
s1.2.coord <- data.table(Cell = rownames(s1b@images$slice1.2@coordinates), s1b@images$slice1.2@coordinates)
s1.3.coord <- data.table(Cell = rownames(s1b@images$slice1.3@coordinates), s1b@images$slice1.3@coordinates)
s1.4.coord <- data.table(Cell = rownames(s1b@images$slice1.4@coordinates), s1b@images$slice1.4@coordinates)
s1.5.coord <- data.table(Cell = rownames(s1b@images$slice1.5@coordinates), s1b@images$slice1.5@coordinates)

s1.coord <- rbindlist(list(s1.0.coord, 
                           s1.1.coord, 
                           s1.2.coord, 
                           s1.3.coord, 
                           s1.4.coord, 
                           s1.5.coord))

identical(s1b.meta$Cell, s1.coord$Cell)

s1b.meta[, `:=` (imagerow = s1.coord$imagerow,
                 imagecol = s1.coord$imagecol,
                 row = s1.coord$row,
                 col = s1.coord$col)]

#make UMAP
lab.spa.coord <- s1b.meta[, .(X = median(UMAP.1), Y = median(UMAP.2)), by = cell_type]

umap.spat3 <- ggplot(s1b.meta, aes(x = UMAP.1, y = UMAP.2, color = cell_type)) +
  geom_point()+
  geom_label_repel(data = lab.spa.coord, aes(x = X, y = Y, fill = cell_type, label = cell_type), color = "black", size = 6) +
  scale_color_manual(values = cell_cols, name = "Cell type") +
  scale_fill_manual(values = cell_cols, name = "Cell type") +
  theme_bw() +
  theme(legend.title = element_text(size= 16),
        legend.text = element_text(size = 14))


marker.genes <- c("LCN2", "OLFM4", "MUC12", "PIGR", "SLC26A3", "PHGR1", "FABP1", 
                  "AQP8", "SLC9A3", "KRT1", "FLG", "KRT10", "KRTDAP", "KRT2", "IGHG1", 
                  "TCHH", "FADS2", "KRT25", "KRT71", "DCN", "COL1A1", "COL1A2", 
                  "COL3A1", "TNXB", "IGHM", "IGHA1", "IGKC", "JCHAIN", "DES", "TAGLN", 
                  "TPM2", "MYH11", "MYL9", "DCD", "MUCL1", "SCGB2A2", "SCGB1D2", 
                  "SLC12A2", "CCL19", "LCP1", "IL32", "CCL17", "CD74", "MT-ND5", 
                  "MT-ND1", "MT-ND3", "MT-ND6", "MT-CYB")

heat.cell <- s1b@assays$SCT@scale.data[marker.genes,]


top_ha <- HeatmapAnnotation(df = data.frame(Cell_type = s1b.meta$cell_type),
                            col = list(Cell_type = cell_cols))

col_fun <- circlize::colorRamp2(seq(-4, 5, by = 1), PurpleAndYellow(10))


dm1 <- Heatmap(heat.cell,
               col = col_fun(seq(-5,10)),
               name = "expression",
               show_column_names = FALSE,
               top_annotation = top_ha,
               column_split = s1b.meta$cell_type,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               column_title_rot = 90,
               row_names_gp = gpar(fontsize = 10))

dm2 <- grid.grabExpr(draw(dm1))

dot.spat <- dittoDotPlot(s1b, vars = c("CD19", "CD8A", "CD68", "KRT10", "PIGR", "COL1A1", "MYH11", "MT-ND5"), group.by = "cell_type",
                         scale = TRUE) +
  ylab("Cell type") +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

pro.spat <- ggplot(s1b.meta, aes(x = Paper_sample, fill = cell_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cell_cols, name = "Cell type") +
  theme_bw() +
  facet_wrap(~Condition, scales = "free_x") +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  xlab("") +
  ylab("Proportion")

fig6.design <- "AABBB
CCDDD"

fig6 <- wrap_elements(full = umap.spat3) + dm2 + dot.spat + pro.spat + plot_layout(design = fig6.design) + plot_annotation(title = "Figure 6", tag_levels = "A") &  theme(plot.title = element_text(size = 24), plot.tag = element_text(size = 20))
ggsave("Figure 6.pdf", plot = fig6, width = 24, height = 18)
ggsave("Figure 6.png", plot = fig6,width = 24, height = 18)

#Figure 7
skin.cell <- ggplot(s1b.meta[Paper_sample %in% c("Normal skin 2")],
                    aes(x = imagecol, y = imagerow, fill = cell_type)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = cell_cols, name = "Cell type") +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(200, 700), ylim = c(1400, 800)) +
  theme_void() +
  ggtitle("Normal skin") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16))

skin.il17 <- ggplot(s1b.meta[Paper_sample %in% c("Normal skin 2")],
                    aes(x = imagecol, y = imagerow, fill = WP_IL17_SIGNALING_PATHWAY)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(-1000, 5000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(200, 700), ylim = c(1400, 800)) +
  theme_void() +
  theme(legend.position = "none") 

skin.ifng <- ggplot(s1b.meta[Paper_sample %in% c("Normal skin 2")],
                    aes(x = imagecol, y = imagerow, fill = WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, 7000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(200, 700), ylim = c(1400, 800)) +
  theme_void() +
  theme(legend.position = "none")

skin.il6 <- ggplot(s1b.meta[Paper_sample %in% c("Normal skin 2")],
                   aes(x = imagecol, y = imagerow, fill = WP_IL6_SIGNALING_PATHWAY)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0,6000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(200, 700), ylim = c(1400, 800)) +
  theme_void() +
  theme(legend.position = "none")

rash.cell <- ggplot(s1b.meta[Paper_sample %in% c("Skin rash 6")],
                    aes(x = imagecol, y = imagerow, fill = cell_type)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = cell_cols, name = "Cell type") +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(875, 1375), ylim = c(1425, 825))+
  theme_void() +   
  ggtitle("Skin rash") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16))

rash.il17 <- ggplot(s1b.meta[Paper_sample %in% c("Skin rash 6")],
                    aes(x = imagecol, y = imagerow, fill = WP_IL17_SIGNALING_PATHWAY)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(-1000, 5000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(875, 1375), ylim = c(1425, 825))+
  theme_void() +
  theme(legend.position = "none")

rash.ifng <- ggplot(s1b.meta[Paper_sample %in% c("Skin rash 6")],
                    aes(x = imagecol, y = imagerow, fill = WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, 7000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(875, 1375), ylim = c(1425, 825))+
  theme_void() +
  theme(legend.position = "none")

rash.il6 <- ggplot(s1b.meta[Paper_sample %in% c("Skin rash 6")],
                   aes(x = imagecol, y = imagerow, fill = WP_IL6_SIGNALING_PATHWAY)) +
  annotation_raster(img.skin, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0,6000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(875, 1375), ylim = c(1425, 825))+
  theme_void() +
  theme(legend.position = "none")

colitis.cell <- ggplot(s1b.meta[Paper_sample %in% c("Colitis 1")],
                       aes(x = imagecol, y = imagerow, fill = cell_type)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = cell_cols, name = "Cell type", drop = FALSE) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(250, 750), ylim = c(1825, 1225)) +
  theme_void() +
  ggtitle("Colitis") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

colitis.il17 <- ggplot(s1b.meta[Paper_sample %in% c("Colitis 1")],
                       aes(x = imagecol, y = imagerow, fill = WP_IL17_SIGNALING_PATHWAY)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(-1000, 5000), name = "IL17 signaling\npathway") +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(250, 750), ylim = c(1825, 1250)) +
  theme_void() +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

colitis.ifng <- ggplot(s1b.meta[Paper_sample %in% c("Colitis 1")],
                       aes(x = imagecol, y = imagerow, fill = WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, 7000), name = "IFNG signaling\npathway") +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(250, 750), ylim = c(1825, 1250)) +
  theme_void() +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

colitis.il6 <- ggplot(s1b.meta[Paper_sample %in% c("Colitis 1")],
                      aes(x = imagecol, y = imagerow, fill = WP_IL6_SIGNALING_PATHWAY)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0,6000), name = "IL6 signaling\npathway") +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(250, 750), ylim = c(1825, 1250)) +
  theme_void() +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

colon.cell <- ggplot(s1b.meta[Paper_sample %in% c("Colon 1")],
                     aes(x = imagecol, y = imagerow, fill = cell_type)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_manual(values = cell_cols, name = "Cell type", drop = FALSE) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(1000, 1500), ylim = c(2100, 1500)) +
  theme_void() +
  ggtitle("Normal colon")  +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 16))

colon.il17 <- ggplot(s1b.meta[Paper_sample %in% c("Colon 1")],
                     aes(x = imagecol, y = imagerow, fill = WP_IL17_SIGNALING_PATHWAY)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(-1000, 5000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(1000, 1500), ylim = c(2100, 1500)) +
  theme_void() +
  theme(legend.position = "none")

colon.ifng <- ggplot(s1b.meta[Paper_sample %in% c("Colon 1")],
                     aes(x = imagecol, y = imagerow, fill = WP_TYPE_II_INTERFERON_SIGNALING_IFNG)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0, 7000), name = "WP_TYPE_II_INTERFERON\nSIGNALING_IFNG") +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(1000, 1500), ylim = c(2100, 1500)) +
  theme_void() +
  theme(legend.position = "none")

colon.il6 <- ggplot(s1b.meta[Paper_sample %in% c("Colon 1")],
                    aes(x = imagecol, y = imagerow, fill = WP_IL6_SIGNALING_PATHWAY)) +
  annotation_raster(img.col, xmin = 0, xmax = 3000, ymin = -3000, ymax = 0) +
  geom_point(pch = 21, size = 3, alpha = 0.8) +
  scale_fill_viridis_c(option = "turbo", limits = c(0,6000)) +
  theme_bw() +
  scale_y_reverse() +
  coord_cartesian(xlim = c(1000, 1500), ylim = c(2100, 1500)) +
  theme_void() +
  theme(legend.position = "none")

vios.il17 <- ggplot(s1b.meta, aes(x = Condition, y = WP_IL17_SIGNALING_PATHWAY, color = Condition)) +
  geom_violin() +
  geom_quasirandom() +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF","#E64B35FF", "#4DBBD5FF"),
                     name = "Condition") +
  stat_compare_means(comparisons = list(c("Normal skin", "Skin rash"), c("Normal colon", "Colitis")),
                     size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  xlab("") +
  ylab("IL17  Signaling pathway") +
  scale_y_continuous(expand = expansion(mult = c(0,0.2))) 

vios.ifng <- ggplot(s1b.meta, aes(x = Condition, y = WP_TYPE_II_INTERFERON_SIGNALING_IFNG, color = Condition)) +
  geom_violin() +
  geom_quasirandom() +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF","#E64B35FF", "#4DBBD5FF"),
                     name = "Condition") +
  stat_compare_means(comparisons = list(c("Normal skin", "Skin rash"), c("Normal colon", "Colitis")),
                     size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  xlab("") +
  ylab("IFNG Signaling pathway") +
  scale_y_continuous(expand = expansion(mult = c(0,0.2))) 

vios.il6 <- ggplot(s1b.meta, aes(x = Condition, y = WP_IL6_SIGNALING_PATHWAY, color = Condition)) +
  geom_violin() +
  geom_quasirandom() +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF","#E64B35FF", "#4DBBD5FF"),
                     name = "Condition") +
  stat_compare_means(comparisons = list(c("Normal skin", "Skin rash"), c("Normal colon", "Colitis")),
                     size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  xlab("") +
  ylab("IL6 Signaling pathway") +
  scale_y_continuous(expand = expansion(mult = c(0,0.2))) 

fig7 <- (skin.cell | rash.cell | colon.cell | colitis.cell) /
  (skin.il17 | rash.il17 | colon.il17 | colitis.il17) / 
  (skin.ifng | rash.ifng | colon.ifng | colitis.ifng) / 
  (skin.il6 | rash.il6 | colon.il6 | colitis.il6) / 
  (vios.il17 | vios.ifng | vios.il6) + 
  plot_annotation(title = "Figure 7", tag_levels = "A") &  
  theme(plot.title = element_text(size = 24), plot.tag = element_text(size = 20))

ggsave("Figure 7.pdf", plot = fig6.new3,  width = 20, height = 18)
ggsave("Figure 7.png", plot = fig6.new3, width = 20, height = 18)
