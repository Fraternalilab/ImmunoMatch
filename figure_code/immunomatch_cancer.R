# compare ImmunoMatch score distribution in different datasets
load("Documents/Dongjun_Guo/immunomatch_scored_data.RData")

# merge annotation of 
# 1. CLL and Lymphoma samples - IGHV gremlin identity
cll_anno <- "Documents/Dongjun_Guo/cll_sequence_aa_with_info.txt"
cll_anno <- read.table(cll_anno, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(cll_igblast)[which(colnames(cll_igblast) == "v_identity")] <- "videntity" 
cll_anno <- merge(
  cll_anno, cll_igblast[, c("sequence_id", "locus", "videntity")],
  by = "sequence_id", all.x = TRUE, all.y = FALSE, sort = FALSE
)
cll <- merge(
  cll, cll_anno[which(cll_anno$chain_type == "VL"), c("ab_id", "locus", "videntity")],
  by = "ab_id", all.x = TRUE, all.y = FALSE, sort = FALSE
)
cll <- merge(
  cll, cll_anno[which(cll_anno$chain_type == "VH"), c("ab_id", "locus", "videntity")],
  by = "ab_id", all.x = TRUE, all.y = FALSE, sort = FALSE, suffixes = c("_VL", "_VH")
)

lymphoma_anno <- "Documents/Dongjun_Guo/lymphoma_sequence_aa_with_info.txt"
lymphoma_anno <- read.table(lymphoma_anno, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(lymphoma_igblast)[which(colnames(lymphoma_igblast) == "v_identity")] <- "videntity" 
lymphoma_anno <- merge(
  lymphoma_anno, lymphoma_igblast[, c("sequence_id", "locus", "videntity")],
  by = "sequence_id", all.x = TRUE, all.y = FALSE, sort = FALSE
)
lymphoma <- merge(
  lymphoma, lymphoma_anno[which(lymphoma_anno$chain_type == "VL"), c("ab_id", "locus", "videntity")],
  by = "ab_id", all.x = TRUE, all.y = FALSE, sort = FALSE
)
lymphoma <- merge(
  lymphoma, lymphoma_anno[which(lymphoma_anno$chain_type == "VH"), c("ab_id", "locus", "videntity")],
  by = "ab_id", all.x = TRUE, all.y = FALSE, sort = FALSE, suffixes = c("_VL", "_VH")
)

# 2. Tonsil GC
tonsil_anno <- read.csv("Documents/Dongjun_Guo/King_Tonsil_GC_paired.csv", 
                        stringsAsFactors = FALSE)
tonsil <- merge(
  tonsil, tonsil_anno[, c("sample_name", "barcode", "locus")],
  by = c("sample_name", "barcode"), all.x = TRUE, all.y = FALSE, sort = FALSE
)
tonsil$donor <- sapply(tonsil$sample_name, function(x){
  unlist(strsplit(x, split = "_"))[1]
})
tonsil <- tonsil[, c("donor", "sample_name", "pairing_scores_original", 
                     "pairing_scores_k", "pairing_scores_l", "locus")]

# 3. Jaffe (i.e. Naive and memory) - here nothing else is needed - just subset to necessary columns
jaffe <- jaffe[, c("donor_id_h", "cell_type", "pairing_scores_original", 
                   "pairing_scores_k", "pairing_scores_l", "chain_l")]

# 4. Cancer cell lines - merge IGHV gremlin identity
celllines_vident <- read.csv(
  "Documents/Dongjun_Guo/ccle_celllines_filtered_mixcr_results_videntity.csv",
  stringsAsFactors = FALSE
)
celllines <- merge(
  celllines, unique(celllines_vident[, c("Cell_line", "VH", "videntity_VH")]),
  by = c("Cell_line", "VH"), all.x = TRUE, all.y = FALSE, sort = FALSE
)
celllines <- merge(
  celllines, unique(celllines_vident[, c("Cell_line", "VL", "videntity_VL")]),
  by = c("Cell_line", "VL"), all.x = TRUE, all.y = FALSE, sort = FALSE
)
celllines$locus <- sapply(celllines$CL, function(x){
  if(grepl("IGL", x)) return("IGL")
  if(grepl("IGK", x)) return("IGK")
  return(NA)
})
celllines[which(is.na(celllines$locus)), "locus"] <- "IGK" # I checked this sequence manually!
celllines2 <- celllines[, c("cancer_type", "Sample_Name", "pairing_scores_original",
                            "pairing_scores_k", "pairing_scores_l", "locus")]

cll$sample_type <- "leukaemia"
lymphoma$sample_type <- "lymphoma"
cll2 <- cll[, c("sample_type", "ab_id", 
                "pairing_scores_original", "pairing_scores_k", 
                "pairing_scores_l", "locus_VL")]
lymphoma2 <- lymphoma[, c("sample_type", "ab_id", 
                          "pairing_scores_original", "pairing_scores_k", 
                          "pairing_scores_l", "locus_VL")]

# Merge all the datasets together
library(plyr)
all_data <- list(
  jaffe,
  celllines2,
  cll2,
  lymphoma2,
  tonsil
)
all_data[[5]][, 2] <- "GC"
colnames(all_data[[5]])[2] <- "sample_type"
colnames(all_data[[1]]) <- c("sample_name", "sample_type", "pairing_score_original", "pairing_score_k", "pairing_score_l", "chain_l")
colnames(all_data[[2]]) <- c("sample_type", "sample_name", "pairing_score_original", "pairing_score_k", "pairing_score_l", "chain_l")
colnames(all_data[[3]]) <- c("sample_type", "sample_name", "pairing_score_original", "pairing_score_k", "pairing_score_l", "chain_l")
colnames(all_data[[4]]) <- c("sample_type", "sample_name", "pairing_score_original", "pairing_score_k", "pairing_score_l", "chain_l")
colnames(all_data[[5]]) <- c("sample_name", "sample_type", "pairing_score_original", "pairing_score_k", "pairing_score_l", "chain_l")
all_data <- do.call("rbind", lapply(all_data, function(tb){
  tb[, c("sample_name", "sample_type", "pairing_score_original", "pairing_score_k", "pairing_score_l", "chain_l")]
}))

# use scores from the lambda model for lambda sequences & scores from the kappa model for kappa sequences
all_data$pairing_score2 <- apply(all_data[, c("pairing_score_k", "pairing_score_l", "chain_l")], MARGIN = 1, function(x){
  if(is.na(x[3])) return(NA)
  if(x[3] == "IGK") return(as.numeric(x[1]))
  if(x[3] == "IGL") return(as.numeric(x[2]))
  return(NA)
})
all_data$sample_type <- factor(
  all_data$sample_type,
  levels = rev(c("naive", "GC", "unswitched_memory", 
             "switched_memory", "leukaemia", "lymphoma")),
  labels = rev(c("Naive", "GC", "Memory", 
             "Memory", "Leukaemia", "Lymphoma"))
)

# Take average per donor
library(plyr)
all_data_avg <- ddply(all_data[which(!is.na(all_data$chain_l)), ], 
                      c("sample_type", "chain_l", "sample_name"),
                      summarise, pairing_score = mean(pairing_score2))

#------------------------------------
# Plotting
library(ggplot2)

# 1. Naive vs GC vs Memory
p <- ggplot(all_data_avg[which(! all_data_avg$sample_type %in% c("Lymphoma", "Leukaemia")), ], 
       aes(y = sample_type, x = pairing_score)) + xlim(0.4, 0.9) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_boxplot() + geom_point() + cowplot::theme_cowplot() +
  facet_wrap(~ chain_l, ncol = 1, scales = "free_x") + 
  scale_y_discrete(name = "") + 
  xlab("ImmunoMatch Pairing Score")
ggsave(plot = p, "Documents/Dongjun_Guo/normal_immunomatch.svg",
       width = 3.5, height=  5)

#------------------------------------
# 2. Cancers - compare scores by cancer subtypes
# manually curated leukaemia/lymphoma subtypes
# now plot scores by subtypes
cancers <- read.csv("Documents/Dongjun_Guo/immunomatch_scored_datasets_cancers.csv",
                    header = TRUE, stringsAsFactors = FALSE)
cancers$cancer_subtype <- factor(
  cancers$cancer_subtype,
  levels = c("Pre-B-ALL", "unmutated CLL", "mantle cell lymphoma", 
             "Burkitt lymphoma", "follicular lymphoma",
             "GCB-DLBCL", "ABC-DLBCL", "mutated CLL",
             "primary effusion lymphoma",
             "lymphoplasmacytic lymphoma"),
  labels = c("Pre-B-ALL", "Unmutated CLL", "Mantle cell lymphoma", 
             "Burkitt lymphoma", "Follicular lymphoma",
             "GCB-DLBCL", "ABC-DLBCL", "Mutated CLL",
             "Primary effusion lymphoma",
             "Lymphoplasmacytic lymphoma")
)
library(ggplot2)
p <- ggplot(cancers, aes(y = cancer_subtype, x = pairing_score2)) +
  geom_boxplot(outlier.shape = NA) + geom_point() + 
  scale_y_discrete(limits = rev, name = "") +
  cowplot::theme_cowplot() + xlab("ImmunoMatch Pairing Score")
ggsave(plot = p, "Documents/Dongjun_Guo/cancers_immunomatch.svg",
       width = 5, height = 5.6)

