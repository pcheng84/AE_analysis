#visium with Seurat
library(Seurat)
library(ggplot2)
library(patchwork)
library(harmony)
library(dplyr)
library(ggsci)
library(UCell)
library(SignatuR)
library(msigdbr)
library(escape)
library(dittoSeq)
library(data.table)
library(ggbeeswarm)
library(ggpubr)

cols11 <- pal_d3("category20")(11)
names(cols11) <- unique(s1a$cell_simp)

s1a <- readRDS("Visium_Seurat_harmony_GEO.rds")

Idents(s1a) <- s1a$cell_simp

DimPlot(s1a, group.by = "cell_simp", label = TRUE, cols = cols11, label.size = 6, label.box = FALSE, repel = TRUE)
ggsave("Figure6A.pdf", width = 12, height = 8)

cluster.markers.simp <- FindAllMarkers(s1a, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster.markers.simp %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
data.table::fwrite(cluster.markers.simp, "230530_Visium/Seurat_FindAllMarkers_cell_simp.txt", sep = "\t")
DoHeatmap(s1a, features = top10$gene, group.colors = cols11, angle = 90) + NoLegend()
ggsave("Figure6B.pdf", width = 12, height = 10)

SpatialDimPlot(s1a, group.by = "cell_simp", images = "slice1",
               cols = cols11) + guides(fill = guide_legend(title = "Cell type", override.aes = list(size = 5)))
ggsave("Figure6E_right_cell.pdf", width = 12, height = 8)

SpatialFeaturePlot(s1a, features = "WP_IL17_SIGNALING_PATHWAY",
                   images = "slice1", min.cutoff = -1000, max.cutoff = 4000) + theme(legend.position = "right")
ggsave("Figure6E_right_IL17.pdf", width = 12, height = 8)

SpatialFeaturePlot(s1a, features = "WP_TYPE_II_INTERFERON_SIGNALING_IFNG",
                   images = "slice1", min.cutoff = 0, max.cutoff = 6000) + theme(legend.position = "right")
ggsave("Figure6E_right_IFNG.pdf", width = 12, height = 8)

SpatialDimPlot(s1a, group.by = "cell_simp", images = "slice1.5",
               cols = cols11) + guides(fill = guide_legend(title = "Cell type", override.aes = list(size = 5)))
ggsave("Figure6E_left_cell.pdf", width = 12, height = 8)

SpatialFeaturePlot(s1a, features = "WP_IL17_SIGNALING_PATHWAY",
                   images = "slice1.5", min.cutoff = -1000, max.cutoff = 4000) + theme(legend.position = "right")
ggsave("Figure6E_left_IL17.pdf", width = 12, height = 8)

SpatialFeaturePlot(s1a, features = "WP_TYPE_II_INTERFERON_SIGNALING_IFNG",
                   images = "slice1.5", min.cutoff = 0, max.cutoff = 6000) + theme(legend.position = "right")
ggsave("Figure6E_left_IFNG.pdf", width = 12, height = 8)

SpatialFeaturePlot(s1a, features = "WP_IL17_SIGNALING_PATHWAY",
                   images = "slice1.6")

dittoDotPlot(s1a, c("CD19", "CD8A", "CD68", "IL17A", "IFNG", 
                    "LAG3", "TIGIT", "PDCD1", "CTLA4", "FOXP3"
                    ),
             group.by = c("cell_simp")) +
  ylab("Cell type") +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Figure6C", width = 12, height = 10)

s1a.meta <- s1a@meta.data

ggplot(s1a.meta[!is.na(Condition2)], aes(x = Paper_sample2, fill = cell_simp)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = cols11, name = "Cell type") +
  theme_bw() +
  facet_wrap(~Condition3, scales = "free_x") +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  xlab("") +
  ylab("Proportion")
ggsave("Figure6D.pdf", width = 12, height = 10)

s1a.meta[, Condition3 := factor(Condition2, levels = c("Normal skin", "Rash", "Normal colon", "Colitis"), labels = c("Normal skin", "Skin rash", "Normal colon", "Colitis"))]

ggplot(s1a.meta[!is.na(Condition2)], aes(x = Condition3, y = WP_IL17_SIGNALING_PATHWAY, color = Condition3)) +
  geom_violin() +
  geom_quasirandom() +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF","#E64B35FF", "#4DBBD5FF"),
                     name = "Condition") +
  stat_compare_means(comparisons = list(c("Normal skin", "Skin rash"), c("Normal colon", "Colitis")),
                     size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  xlab("") +
  ylab("IL17  Signaling pathway")
ggsave("Figure6F_top.pdf", width = 12, height = 8)

ggplot(s1a.meta[!is.na(Condition2)], aes(x = Condition3, y = WP_TYPE_II_INTERFERON_SIGNALING_IFNG, color = Condition3)) +
  geom_violin() +
  geom_quasirandom() +
  scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF","#E64B35FF", "#4DBBD5FF"),
                     name = "Condition") +
  stat_compare_means(comparisons = list(c("Normal skin", "Skin rash"), c("Normal colon", "Colitis")),
                     size = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  xlab("") +
  ylab("IFNG Signaling pathway")
ggsave("Figure6F_bottom.pdf", width = 12, height = 8)




