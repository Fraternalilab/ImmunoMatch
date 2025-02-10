# plot popscomp statistics
popscomp_all_results <- "Documents/Dongjun_Guo/vcab_human_for_pairing_interface_positions_imgt.csv"
popscomp_all_results <- read.table(popscomp_all_results, sep = ",", stringsAsFactors = FALSE, header = TRUE)
library(plyr)
popscomp_counts <- ddply(
  popscomp_all_results, c("chain_type", "output_numbering_IMGT"),
  summarise, n_interface = length(D_SASA.A.2)
)
popscomp_counts$prop_interface <- popscomp_counts$n_interface / length(unique(popscomp_all_results$iden_code))

popscomp <- 'Documents/Dongjun_Guo/vcab_human_for_pairing_interface_stats.csv'
popscomp <- read.csv(popscomp, stringsAsFactors = FALSE)
popscomp <- reshape2::melt(popscomp, id.vars = c("iden_code"),
                           measure.vars = c("VH_CDR1", "VH_CDR2",
                                            "VH_CDR3", "VH_FWR",
                                            "VL_CDR1", "VL_CDR2",
                                            "VL_CDR3", "VL_FWR"))
popscomp$chain <- factor(
  popscomp$variable, levels = levels(popscomp$variable),
  labels = c("VH", "VH", "VH", "VH", "VL", "VL", "VL", "VL")
)

popscomp$region <- factor(
  popscomp$variable, levels = levels(popscomp$variable),
  labels = c("CDR1", "CDR2", "CDR3", "FWR", "CDR1", "CDR2", "CDR3", "FWR")
)

library(ggplot2)
p <- ggplot(popscomp, aes(x = region, y = value)) + geom_boxplot() +
  facet_wrap(~ chain) + cowplot::theme_cowplot() +
  ylab("Proportion of VH-VL interface\nsurface area (A^2)")
ggsave(width =5.6, height = 3.2, plot = p,
       "Documents/Dongjun_Guo/vcab_human_popscomp_interface-sa.svg")
