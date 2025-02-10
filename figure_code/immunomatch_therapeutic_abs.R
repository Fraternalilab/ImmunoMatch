thera <- read.csv("Documents/Dongjun_Guo/immunomatch_thera_abs/TheraSAbDab_SeqStruc_OnlineDownload_20241205_human_monospecific_approved-active_scored.csv",
                  stringsAsFactors = FALSE)
thera$pairing_scores2 <- apply(thera[, c("VD.LC", "pairing_scores_k", "pairing_scores_l")], MARGIN = 1, function(x){
  if(x[1] == "Lambda") return(as.numeric(x[3])) else return(as.numeric(x[2]))
})
thera_anarci <- read.csv("Documents/Dongjun_Guo/immunomatch_thera_abs/TheraSAbDab_SeqStruc_OnlineDownload_20241205_human_monospecific_approved-active_anarci.csv",
                         stringsAsFactors = FALSE)
thera <- merge(
  thera, thera_anarci[
    , c("Therapeutic", "species_h", "species_l", "videntity_h",
        "videntity_l", "videntity_human_h", "videntity_human_l")
  ],
  by = "Therapeutic", all.x = TRUE, all.y = FALSE, sort = FALSE
)
thera <- thera[which(!is.na(thera$species_h)), ]
# scored random V-J pairings for therapeutic antibodies
randomH <- read.csv("Documents/Dongjun_Guo/immunomatch_thera_abs/TheraSAbDab_SeqStruc_OnlineDownload_20241205_human_monospecific_approved-active_randomHeavy_scored_pident.csv",
                    stringsAsFactors = FALSE)
randomH$pairing_scores2 <- apply(randomH[, c("VD.LC", "pairing_scores_k", "pairing_scores_l")], MARGIN = 1, function(x){
  if(x[1] == "Lambda") return(as.numeric(x[3])) else return(as.numeric(x[2]))
})
randomH <- merge(randomH, thera[, c("Therapeutic", "pairing_scores2")],
                 by = "Therapeutic", all.x = TRUE, all.y = FALSE,
                 sort = FALSE, suffixes = c("", "_WT"))
randomH$pident_group <- ggplot2::cut_width(
  randomH$pident, width = 10, labels = c("<50", "50-60", "60-70", "70-80", "80-90", "90-100"),
  breaks = c(0, 50, 60, 70, 80, 90, 100)
)
randomH$score_diff <- randomH$pairing_scores2 - randomH$pairing_scores2_WT

library(ggplot2)
ggplot(randomH[which(randomH$Therapeutic == "Trastuzumab"), ],
       aes(x = pident, y = pairing_scores2)) + geom_point() +
  geom_hline(aes(yintercept = unique(pairing_scores2_WT))) +
  cowplot::theme_cowplot() + ylim(0, 1)

library(plyr)
pair_stats <- ddply(
  randomH, c("Therapeutic", "pairing_scores2_WT", "pident_group"), summarise, 
  prop_pair = mean(pairing_scores2 >= 0.5, na.rm = TRUE), 
  n_pair = sum(pairing_scores2 >= 0.5), 
  n = length(VH_random) 
)

best_random <- ddply(
  randomH, c("Therapeutic", "pairing_scores2_WT"), summarise, 
  best_random = VH_random[which.max(pairing_scores2)],
  max_score = max(pairing_scores2),
  best_pident = pident[which.max(pairing_scores2)],
  n_pair = sum(pairing_scores2 >= 0.5),
  total_random = length(pairing_scores2)
)

highest_pident_random <- ddply(
  randomH, c("Therapeutic", "HeavySequence", "X100..SI.Structure",
             "X99..SI.Structure", "X95.98..SI.Structure"), 
  summarise,
  score_wt = unique(pairing_scores2_WT),
  highest_pident = VH_random[which.max(pident)],
  score = pairing_scores2[which.max(pident)],
  pident = pident[which.max(pident)],
  n_pair = sum(pairing_scores2 >= 0.5),
  total_random = length(pairing_scores2)
)
highest_pident_random$score_diff <- highest_pident_random$score - highest_pident_random$score_wt
highest_pident_random$pident_group <- ggplot2::cut_width(
  highest_pident_random$pident, width = 10, 
  labels = c("60-70", "70-80", "80-90", "90-100"),
  breaks = c(60, 70, 80, 90, 100)
)
highest_pident_random$pident_group <- factor(
  highest_pident_random$pident_group,
  levels = rev(levels(highest_pident_random$pident_group))
)

ggplot(highest_pident_random, aes(x = pident, y = abs(score_diff))) + 
  geom_point()

library(ggplot2)
library(ggdist)
p <- ggplot(highest_pident_random, 
       aes(x = pident_group, y = abs(score_diff))) + 
  geom_boxplot(width = 0.15) + 
  stat_halfeye(adjust = 0.5, justification = -0.2,
               .width = 0, point_colour = NA, width=0.8) +
  cowplot::theme_cowplot() +
  ylab("|Score_random - Score_WT|") + 
  xlab("Sequence identity WT vs random (%)")
ggsave(plot = p, width = 3.7, height = 3.1,
       "Documents/Dongjun_Guo/immunomatch_thera_abs/seq_id-vs_score-diff.svg")
# random kappa
randomK <- read.csv("Documents/Dongjun_Guo/immunomatch_thera_abs/TheraSAbDab_SeqStruc_OnlineDownload_20241205_human_monospecific_approved-active_randomKappa_scored_pident.csv",
                    stringsAsFactors = FALSE)
randomK$pident_group <- ggplot2::cut_width(
  randomK$pident, width = 10, labels = c("<50", "50-60", "60-70", "70-80", "80-90", "90-100"),
  breaks = c(0, 50, 60, 70, 80, 90, 100)
)
pair_stats <- ddply(
  randomK, c("Therapeutic", "VD.LC", "pident_group"), summarise, 
  prop_pair = mean(pairing_scores_k >= 0.5, na.rm = TRUE), 
  n_pair = sum(pairing_scores_k >= 0.5), 
  n = length(VKappa_random) 
)
best_random <- ddply(
  randomK, c("Therapeutic", "VD.LC"), summarise, 
  best_random = VKappa_random[which.max(pairing_scores_k)],
  max_score = max(pairing_scores_k),
  best_pident = pident[which.max(pairing_scores_k)],
  n_pair = sum(pairing_scores_k >= 0.5),
  total_random = length(pairing_scores_k)
)
highest_pident_random <- ddply(
  randomK, c("Therapeutic", "VD.LC"), summarise, 
  highest_pident = VKappa_random[which.max(pident)],
  score = pairing_scores_k[which.max(pident)],
  pident = pident[which.max(pident)],
  n_pair = sum(pairing_scores_k >= 0.5),
  total_random = length(pairing_scores_k)
)

randomL <- read.csv("Documents/Dongjun_Guo/immunomatch_thera_abs/TheraSAbDab_SeqStruc_OnlineDownload_20241205_human_monospecific_approved-active_randomLambda_scored_pident.csv",
                    stringsAsFactors = FALSE)
randomL$pident_group <- ggplot2::cut_width(
  randomL$pident, width = 10, labels = c("<50", "50-60", "60-70", "70-80", "80-90", "90-100"),
  breaks = c(0, 50, 60, 70, 80, 90, 100)
)
pair_stats <- ddply(
  randomL, c("Therapeutic", "VD.LC", "pident_group"), summarise, 
  prop_pair = mean(pairing_scores_l >= 0.5, na.rm = TRUE), 
  n_pair = sum(pairing_scores_l >= 0.5), 
  n = length(VLambda_random) 
)
best_random <- ddply(
  randomL, c("Therapeutic", "VD.LC"), summarise, 
  best_random = VLambda_random[which.max(pairing_scores_l)],
  max_score = max(pairing_scores_l),
  best_pident = pident[which.max(pairing_scores_l)],
  n_pair = sum(pairing_scores_l >= 0.5),
  total_random = length(pairing_scores_l)
)
highest_pident_random <- ddply(
  randomL, c("Therapeutic", "VD.LC"), summarise, 
  highest_pident = VLambda_random[which.max(pident)],
  score = pairing_scores_l[which.max(pident)],
  pident = pident[which.max(pident)],
  n_pair = sum(pairing_scores_l >= 0.5),
  total_random = length(pairing_scores_l)
)


