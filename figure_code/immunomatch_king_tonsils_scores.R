# Analyse King et al. CSR vs SHM vs Pairing
library(Seurat)
king <- "/Volumes/Seagate4TB/GLT_datasets/King_Tonsil/SEURAT_OBJECTS/HumanTonsil_BCells_scRNA_IGHC_metacell_SeuratObject.rds"
king <- readRDS(king)
king <- UpdateSeuratObject(king)

# merge pairing score to the Seurat object
king_pairing <- read.csv(
  "Documents/Dongjun_Guo/King_Tonsil_allBcells_paired_immunomatch.csv",
  stringsAsFactors = FALSE
)
rownames(king_pairing) <- king_pairing$barcode
king_pairing$pairing_scores2 <- apply(
  king_pairing[, c("locus", "pairing_scores_k", "pairing_scores_l")], MARGIN = 1, function(x){
    if(x[1] == "IGK") return(as.numeric(x[2]))
    if(x[1] == "IGL") return(as.numeric(x[3]))
  }
)
king <- Seurat::AddMetaData(king, king_pairing[, "pairing_scores2", drop = FALSE])

# random control - pseudonegative generated for each cell
# by random shuffling of L chains between cells 
# of same cell type
king_pairing_random_celltype <- read.csv(
  "Documents/Dongjun_Guo/King_Tonsil_allBcells_paired_pairing_scores_random.csv",
  stringsAsFactors = FALSE
)
king_pairing_random_celltype <- king_pairing_random_celltype[king_pairing_random_celltype$label == 0, ]
rownames(king_pairing_random_celltype) <- king_pairing_random_celltype$old_H_ab_id
colnames(king_pairing_random_celltype)[9] <- "pairing_score_random"
king <- Seurat::AddMetaData(king, king_pairing_random_celltype[, "pairing_score_random", drop = FALSE])

# random control2 - same as the first random control but
# this time not restricting shuffling within same
# cell type
king_pairing_random <- read.csv(
  "Documents/Dongjun_Guo/King_Tonsil_allBcells_paired_pairing_scores_random2.csv",
  stringsAsFactors = FALSE
)
king_pairing_random <- king_pairing_random[king_pairing_random$label == 0, ]
rownames(king_pairing_random) <- king_pairing_random$old_H_ab_id
colnames(king_pairing_random)[9] <- "pairing_score_random2"
king <- Seurat::AddMetaData(king, king_pairing_random[, "pairing_score_random2", drop = FALSE])

# only retain those with pairing_scores2
king <- subset(king, cells=Cells(king)[which(!is.na(king$pairing_scores2) &
                                               !is.na(king$pairing_score_random) &
                                               !is.na(king$pairing_score_random2))])
library(ggplot2)
king$ISOTYPE2 <- factor(
  king$ISOTYPE, levels = c("IgM", "IgD", "IgG1", "IgG2",
                            "IgG3", "IgG4", "IgA1", "IgA2"),
  labels = c("IgM", "IgD", "IgG", "IgG", "IgG", "IgG", 
             "IgA", "IgA")
)
king$ISOTYPE2 <- factor(king$ISOTYPE2, levels = c("IgA", "IgG", "IgD", "IgM"))
king$SHM_CAT <- ggplot2::cut_interval(
  as.numeric(king$IGH_MU_FREQ), n = 4,
  breaks = c(0, 0.01, 0.05, 0.1, 0.3), labels = FALSE
)
king$SHM_CAT <- factor(
  king$SHM_CAT, levels = 1:4,
  labels = c("99-100", "95-99", "90-95", "<90")
)
king$SHM_CAT <- factor(
  king$SHM_CAT, levels = c("<90", "90-95", "95-99", "99-100")
)
g <- ggplot(reshape2::melt(king@meta.data[which(!is.na(king$ISOTYPE2)), ], 
                      id.vars = "ISOTYPE2", 
                      measure.vars = c("pairing_scores2", 
                                       "pairing_score_random2")), 
       aes(y = ISOTYPE2, x = value, color = variable)) +
  geom_boxplot(position = position_dodge(width = 1)) + 
  scale_color_manual(values = c("black", "grey"),
                     labels = c("Observed H-L pairs",
                                "H paired with randomly\nshuffled L chains")) + 
  cowplot::theme_cowplot() + ylab("BCR isotype") +
  xlab("ImmunoMatch pairing score")
ggsave(plot = g, width = 5.6, height = 2.4,
       filename = "Documents/Dongjun_Guo/King_Tonsil_allBcells_csr.svg")


g <- ggplot(reshape2::melt(king@meta.data[which(!is.na(king$SHM_CAT)), ], 
                      id.vars = "SHM_CAT", 
                      measure.vars = c("pairing_scores2", 
                                       "pairing_score_random2")), 
       aes(y = SHM_CAT, x = value, color = variable)) +
  geom_boxplot(position = position_dodge(width = 1)) + 
  scale_color_manual(values = c("black", "grey"),
                     labels = c("Observed H-L pairs",
                                "H paired with randomly\nshuffled L chains")) + 
  cowplot::theme_cowplot() + ylab("H chain identity to germline (%)") +
  xlab("ImmunoMatch pairing score")
ggsave(plot = g, width = 5.6, height = 2.4,
       filename = "Documents/Dongjun_Guo/King_Tonsil_allBcells_shm.svg")

# scatter plot of MU_FREQ vs pairing_scores2 and group by cell-type
# for simplicity here I try to simplify the labels
king$CellType2 <- factor(
  king$CellType,
  levels = levels(king$CellType),
  labels = c("Naive", "Activated", "preGC", "GC", "GC", "GC", "GC",
             "(pre-)Plasmablast", "(pre-)Plasmablast", "Memory", "Memory", "Cycling")
)

p <- lapply(levels(king$CellType2), function(x){
  color_list <- scales::hue_pal()(length(levels(king$CellType2)))
  names(color_list) <- levels(king$CellType2)
  o <- ggplot(king@meta.data[which(!is.na(king$SHM_CAT) &
                                king$CellType2 == x), ], 
         aes(x = IGH_MU_FREQ, y = pairing_scores2, colour = CellType2)) + 
    geom_point(size = 0.7) + ggtitle(NULL, subtitle = x) +
    scale_color_manual(values = color_list) +
    scale_y_continuous("ImmunoMatch\npairing score", limits = c(0, 1)) +
    scale_x_continuous(labels = c("100%", "95%", "90%", "85%", "80%"),
                       breaks = c(0, 0.05, 0.1, 0.15, 0.2),
                       name = "H chain identity to germline (%)",
                       limits = c(0, 0.22)) +
    theme_bw() + 
    theme(legend.position = "none")
  ggExtra::ggMarginal(o, type = "density", margins = "y", groupColour = TRUE)
})
ggsave("Documents/Dongjun_Guo/shm_vs_pairing_celltype.svg",
       width = 8.5, height = 7, plot = cowplot::plot_grid(plotlist = p))
# statistics
summary(
  lmerTest::lmer(pairing_scores2 ~ IGH_MU_FREQ + (1 | CellType2),
                 data = king@meta.data[which(!is.na(king$SHM_CAT)), ])
)
