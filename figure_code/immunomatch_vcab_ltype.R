# plot sequence comparison of kappa and lambda from vcab
# calculated using mmseqs2
# I extracted the ones labelled 'Homo_Sapiens' in vcab
# and use column VL_Seq
# mmseqs easy-search final_vcab_human_VL.fasta final_vcab_human_VL.fasta final_vcab_human_VL_mmseqs2.m8 tmp --alignment-mode 3 --exhaustive-search --filter-msa 0 --rescore-mode 3

mmseqs_result <- read.table(
  "Documents/Dongjun_Guo/final_vcab_human_VL_mmseqs2.m8",
  sep = "\t", header = FALSE, stringsAsFactors = FALSE
)
final_vcab <- read.csv(
  "Documents/Dongjun_Guo/final_vcab_with_V_coor.csv",
  stringsAsFactors = FALSE
)
final_vcab <- final_vcab[which(final_vcab['Species'] == "Homo_Sapiens"), ]
final_vcab$Ltype2 <- sapply(final_vcab$Ltype, function(x){
  unlist(strsplit(x, split = "(", fixed = TRUE))[1]
})
final_vcab <- final_vcab[
  order(final_vcab$Ltype2, final_vcab$iden_code),
]
final_vcab$iden_code <- factor(
  final_vcab$iden_code, levels = final_vcab$iden_code
)
mmseqs_result[, 1] <- factor(
  mmseqs_result[, 1], levels = levels(final_vcab$iden_code)
)
mmseqs_result[, 2] <- factor(
  mmseqs_result[, 2], levels = levels(final_vcab$iden_code)
)
mmseqs_result <- reshape2::acast(
  V1 ~ V2, value.var = "V3", data = mmseqs_result
)
mmseqs_result <- mmseqs_result[levels(final_vcab$iden_code),
                               levels(final_vcab$iden_code)]
mmseqs_result[upper.tri(mmseqs_result, diag = TRUE)] <- NA
svg("Documents/Dongjun_Guo/vcab_ltype_seqid_heatmap.svg", 
    width = 4.5, height = 4)
fields::image.plot(mmseqs_result, useRaster=  TRUE, col = hcl.colors(100, palette = "YlOrRd", rev = TRUE))
dev.off()

k_k_pident <- c(mmseqs_result[final_vcab[final_vcab$Ltype2 == "kappa", "iden_code"],
                              final_vcab[final_vcab$Ltype2 == "kappa", "iden_code"]])
l_l_pident <- c(mmseqs_result[final_vcab[final_vcab$Ltype2 == "lambda", "iden_code"],
                              final_vcab[final_vcab$Ltype2 == "lambda", "iden_code"]])
k_l_pident <- c(mmseqs_result[final_vcab[final_vcab$Ltype2 == "lambda", "iden_code"],
                              final_vcab[final_vcab$Ltype2 == "kappa", "iden_code"]])
k_k_pident <- mean(k_k_pident, na.rm = TRUE)
l_l_pident <- mean(l_l_pident, na.rm = TRUE)
k_l_pident <- mean(k_l_pident, na.rm = TRUE)
