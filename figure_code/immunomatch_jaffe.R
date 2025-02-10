# map clonotype clustering information to jaffe data

jaffe_clonotypes <- list(
  "Documents/Dongjun_Guo/donor_1_filtered_contig_annotations_IGH_clone-pass.tsv",
  "Documents/Dongjun_Guo/donor_2_filtered_contig_annotations_IGH_clone-pass.tsv",
  "Documents/Dongjun_Guo/donor_4_filtered_contig_annotations_IGH_clone-pass.tsv"
)
jaffe_clonotypes <- lapply(jaffe_clonotypes, read.table,
                           stringsAsFactors = FALSE, sep = "\t",
                           header = TRUE)
jaffe_clonotypes <- do.call("rbind", jaffe_clonotypes)

jaffe_immunomatch <- read.csv(
  "Documents/Dongjun_Guo/jaffe_2022_except_donor3_paired_all_immunomatch.csv",
  stringsAsFactors = FALSE
)

library(plyr)
clone_info <- ddply(
  jaffe_clonotypes, c("donor_id", "clone_id"), summarise,
  clone_size = length(sequence_id),
  celltypes_in_clone = paste(sort(unique(celltype)), collapse = ",")
)
jaffe_immunomatch <- merge(
  jaffe_immunomatch,
  jaffe_clonotypes[, c("donor_id", "file_id", "celltype", "barcode", "sequence_id", "clone_id")],
  by.x = c("donor_id_h", "file_id_h", "cell_type", "barcode_h", "contig_id_h"),
  by.y = c("donor_id", "file_id", "celltype", "barcode", "sequence_id"),
  all.x = TRUE, all.y = FALSE, sort = FALSE
)
jaffe_immunomatch <- merge(
  jaffe_immunomatch,
  clone_info,
  by.x = c("donor_id_h", "clone_id"),
  by.y = c("donor_id", "clone_id"),
  all.x = TRUE, all.y = FALSE, sort = FALSE
)
clone_info2 <- ddply(
  jaffe_immunomatch, c("donor_id_h", "clone_id"), summarise,
  n_cells_in_clone = length(unique(paste(file_id_h, "_", barcode_h))),
  clone_ltype = paste(sort(unique(chain_l)), collapse = ",")
) 
jaffe_immunomatch <- merge(
  jaffe_immunomatch,
  clone_info2,
  by = c("donor_id_h", "clone_id"),
  all.x = TRUE, all.y = FALSE, sort = FALSE
)
write.csv(jaffe_immunomatch, col.names = TRUE, row.names = FALSE,
          "Documents/Dongjun_Guo/jaffe_2022_except_donor3_paired_all_immunomatch_clonotypes.csv")

# stats: for each clone with both naive and memory,
# is the memory always higher than naive in the score?
jaffe_immunomatch$pairing_scores2 <- apply(
  jaffe_immunomatch[, c("chain_l", "pairing_scores_k", "pairing_scores_l")], MARGIN = 1, function(x){
    if(x[1] == "IGK") return(as.numeric(x[2]))
    if(x[1] == "IGL") return(as.numeric(x[3]))
  }
)
jaffe_stats <- ddply(
  jaffe_immunomatch[which(grepl("naive,", jaffe_immunomatch$celltypes_in_clone)), ],
  c("donor_id_h", "clone_id", "clone_size", "n_cells_in_clone",
    "celltypes_in_clone"), summarise,
  naive_score = median(pairing_scores2[which(cell_type == "naive")]),
  memory_score = median(pairing_scores2[which(cell_type != "naive")])
)
jaffe_stats$score_diff <- jaffe_stats$memory_score - jaffe_stats$naive_score

# Probably want to prioritise for building clonotype trees:
# * only 1 ltype
# * with naive + either/both (un)switched mem
# * 10-100 H sequences?
jaffe_immunomatch_priority <- jaffe_immunomatch[
  which(jaffe_immunomatch$clone_ltype != "IGK,IGL" &
          grepl("naive,", jaffe_immunomatch$celltypes_in_clone) &
          jaffe_immunomatch$clone_size > 10 &
          jaffe_immunomatch$clone_size <= 100), 
]
write.csv(jaffe_immunomatch_priority, col.names = TRUE, row.names = FALSE,
          "Documents/Dongjun_Guo/jaffe_2022_except_donor3_paired_all_immunomatch_clonotypes_prioritised.csv")

# BrepPhylo clonotype trees on the selected clones
jaffe_immunomatch_priority <- read.csv(
  "Documents/Dongjun_Guo/jaffe_2022_except_donor3_paired_all_immunomatch_clonotypes_prioritised.csv",
  stringsAsFactors = FALSE
)
# generate output folder
outputFolder <- path.expand( "Documents/Dongjun_Guo/jaffe_2022_BrepPhylo_trees" )
dir.create( outputFolder, showWarnings = FALSE )
dnapars_executable <- "/Users/josefng/phylip-3.695/exe/dnapars.app/Contents/MacOS/dnapars"

jaffe_immunomatch_priority$vh_seq <- apply(
  jaffe_immunomatch_priority[, c("fwr1_nt_h", "cdr1_nt_h",
                                 "fwr2_nt_h", "cdr2_nt_h",
                                 "fwr3_nt_h")], 
  MARGIN = 1, 
  function(x){
    paste(x, collapse = "")
  }
)
jaffe_immunomatch_priority$clone_id2 <- apply(
  jaffe_immunomatch_priority[, c("donor_id_h", "clone_id")], MARGIN =1 , function(x){
    paste(x, collapse = "_clone")
  }
)

library(BrepPhylo)
batch_results <- doBatchCloneAnalysis( 
  jaffe_immunomatch_priority,
  outputFolder = outputFolder,
  species = "Homo sapiens",
  sequence_column = "vh_seq",
  cloneID_column = "clone_id2",
  IGHVgeneandallele_column = "v_gene_h",
  label_column = "barcode_h",
  colourColumn = "cell_type",
  colourList = list("naive" = "black",
                    "switched_memory" = "red",
                    "unswitched_memory" = "red"),
  plotFormat = "pdf",
  phyloTreeType = "dnapars",
  germlineSet = "Documents/Dongjun_Guo/imgt_human_vgenes_20241231.fasta",
  phyloTreeOptions = list( "executable" = dnapars_executable ),
  makeArboTree = FALSE,
  useTempDir = FALSE 
)
saveRDS(batch_results, paste0(outputFolder, '/jaffe_batch_analysis.rds'))

batch_results <- readRDS(paste0(outputFolder, '/jaffe_batch_analysis.rds'))
germline_distance <- lapply(names(batch_results), function(x) {
  tb <- batch_results[[x]]$distances
  if(!is.null(tb)){
    tb$donor_id <- unlist(strsplit(x, split = "_clone"))[1]
    tb$clone_id <- gsub(" ", "", unlist(strsplit(x, split = "_clone"))[2])
    tb
  } else return(NULL)
})
germline_distance <- do.call("rbind", germline_distance)
germline_distance <- merge(
  germline_distance, 
  jaffe_immunomatch_priority[, c("donor_id_h", "clone_id", "barcode_h", "pairing_scores2", "cell_type")],
  by.x = c("donor_id", "clone_id", "SeqID"),
  by.y = c("donor_id_h", "clone_id", "barcode_h"),
  all.x = TRUE, all.y = FALSE, sort = FALSE
)
library(plyr)
germline_distance_stats <- ddply(
  germline_distance, c("donor_id", "clone_id"), summarise,
  r = cor(distFromGermline, pairing_scores2, method = "pearson"),
  clone_size = length(unique(SeqID)),
  n_naive = sum(cell_type == "naive"),
  n_memory = sum(cell_type != "naive"),
  min_score = min(pairing_scores2, na.rm = TRUE),
  max_score = max(pairing_scores2, na.rm = TRUE)
)
germline_distance_stats <- germline_distance_stats[
  which(!is.na(germline_distance_stats$r)), 
]
germline_distance_stats$prop_naive <- germline_distance_stats$n_naive / germline_distance_stats$clone_size
germline_distance_stats$prop_memory <- germline_distance_stats$n_memory / germline_distance_stats$clone_size
germline_distance_stats <- germline_distance_stats[
  which(germline_distance_stats$prop_naive >= 0.3 &
          germline_distance_stats$prop_naive <= 0.7 &
          germline_distance_stats$min_score > 0.5), 
]

library(treeio)
library(ggplot2)
library(ggtree)
clone100526 <- jaffe_immunomatch_priority[
  which(jaffe_immunomatch_priority$donor_id_h == "donor_1" &
          jaffe_immunomatch_priority$clone_id == "100526"),
]
tree_clone100526 <- "Documents/Dongjun_Guo/jaffe_2022_BrepPhylo_trees/clone_donor_1_clone100526/clone_donor_1_clone100526.tree"
tree_clone100526 <- read.tree(tree_clone100526)
clone100526$node <- as.numeric(factor(
  clone100526$barcode_h, 
  levels = c(tree_clone100526$tip.label,
             tree_clone100526$node.label)
))
p <- ggtree(tree_clone100526) %<+% clone100526 +
  geom_tippoint(aes(color = pairing_scores2)) +
  theme_tree() +
  scale_color_viridis_c(option = "B", name = "ImmunoMatch\npairing\nscore")
ggsave(plot = p, width = 3.4, height = 5.4,
       filename = "Documents/Dongjun_Guo/jaffe_donor1_clone100526_tree.svg")

clone124219 <- jaffe_immunomatch_priority[
  which(jaffe_immunomatch_priority$donor_id_h == "donor_1" &
          jaffe_immunomatch_priority$clone_id == "124219"),
]
tree_clone124219 <- "Documents/Dongjun_Guo/jaffe_2022_BrepPhylo_trees/clone_donor_1_clone124219/clone_donor_1_clone124219.tree"
tree_clone124219 <- read.tree(tree_clone124219)
clone124219$node <- as.numeric(factor(
  clone124219$barcode_h, 
  levels = c(tree_clone124219$tip.label,
             tree_clone124219$node.label)
))

ggtree(tree_clone124219) %<+% clone124219 +
  geom_tippoint(aes(color = pairing_scores2)) +
  theme_tree() +
  scale_color_viridis_c(option = "B", name = "ImmunoMatch\npairing\nscore")
