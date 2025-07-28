#downloading a package with MCPcounter (cell deconvolution tool)
library(devtools)
library(dplyr)
library(tidyverse)
library(MCPcounter)
library(rstatix)

#loading metadata from the patients (as in serrated_RNA.R document)
metadata_serrated <- as.data.frame(read_csv("~/Desktop/serrated seq/RNAseq final/RNAmeta_serrated.csv"))

#removing rows with samples that have no RNAseq samples available (patient ID 6341) or have undefined histology (946, 6894) 
rows_to_remove <- which(metadata_serrated$`Patient ID` == 946 | 
                          metadata_serrated$`Patient ID` == 6341 | 
                          metadata_serrated$`Patient ID` == 6894)
metadata_serrated <- metadata_serrated[-rows_to_remove, ]

#extracting only interesting metadata 
metadata_serrated <- data.frame(
  names = metadata_serrated$`RNA Name`,
  condition = metadata_serrated$samples,
  histology = metadata_serrated$histology,
  subject = metadata_serrated$`pair ID`,
  batch = metadata_serrated$`Sample batch`,
  stringsAsFactors = FALSE
)
view(metadata_serrated)

#loading batch-adjusted expression matrix, created in serrated_RNA.R document 
#serrated_matrix <- read_csv("~/Desktop/serrated seq/RNAseq final/batchadjusted_serrated_expression.csv")
serrated_matrix <- read_csv("~/Desktop/serrated seq/RNAseq final/normalized_counts_matrix.csv")
matrix_genes <- serrated_matrix$...1
serrated_matrix <- serrated_matrix[ ,-1]
rownames(serrated_matrix) <- matrix_genes

#running MCPcounter
MCP_serrated <- MCPcounter.estimate(expression = serrated_matrix, 
                    featuresType = "HUGO_symbols",
                    )
#prepare the MCP data to join with metadata
MCP_serrated <- as.data.frame(t(MCP_serrated))
MCP_serrated$names <- rownames(MCP_serrated)
str(MCP_serrated)

#join MCP and metadata
mcp_counter <- left_join(metadata_serrated, MCP_serrated, by = "names")

#converting to long format for visualisation
mcp_counter_long <- mcp_counter |> pivot_longer(-c(names, condition, histology, subject, batch), names_to='cell_type', values_to='cell_estimate')

#generating a boxplot for immune cell type deconvolution 
ggplot(mcp_counter_long, aes(x=histology, y=cell_estimate, fill=condition)) +
  geom_boxplot(outlier.shape=NA) + # Boxplot with transparency
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.7) +
  facet_wrap(~cell_type, nrow = 2,  scales = "free_y") +  # One boxplot per gene
  labs(title="MCP counter cell type deconvolution",
       y="MCP counter cell abundance estimate") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) +
  theme_bw() +
  scale_fill_brewer(palette="Set2")

#stats for MCP deconvolution ADENOMA normal vs polyp MCP counter
mcp_counter_adenoma <- mcp_counter_long |> filter(histology == "adenoma") |> arrange(cell_type, subject)
stats_ade <- tibble()
for (cell in unique(mcp_counter_adenoma$cell_type)) {
  cell_data <- mcp_counter_adenoma |> filter(cell_type == cell)
  wilcox_result <- wilcox_test(cell_estimate ~ condition, data = cell_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(cell_type = cell)
  stats_ade <- bind_rows(stats_ade, wilcox_result)
}
stats_ade$padj <- p.adjust(stats_ade$p, method = "BH")
view(stats_ade)

#stats for MCP deconvolution SERRATED normal vs polyp
mcp_counter_serrated <- mcp_counter_long |> filter(histology == "serrated") |> arrange(cell_type, subject)
stats_ser <- tibble()
for (cell in unique(mcp_counter_serrated$cell_type)) {
  cell_data <- mcp_counter_serrated |> filter(cell_type == cell)
  wilcox_result <- wilcox_test(cell_estimate ~ condition, data = cell_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(cell_type = cell)
  stats_ser <- bind_rows(stats_ser, wilcox_result)
}
stats_ser$padj <- p.adjust(stats_ser$p, method = "BH")
view(stats_ser)

#stats for MCP deconvolution NORMAL adeneoma vs serrated - no significant differences between normal tissuess from adenomas and form serrated polyps 
mcp_counter_normal <- mcp_counter_long |> filter(condition == "normal") |> arrange(cell_type, subject)
stats_normal <- tibble()
for (cell in unique(mcp_counter_normal$cell_type)) {
  cell_data <- mcp_counter_normal |> filter(cell_type == cell)
  wilcox_result <- wilcox_test(cell_estimate ~ histology, data = cell_data)
  wilcox_result <- wilcox_result |> mutate(cell_type = cell)
  stats_normal <- bind_rows(stats_normal, wilcox_result)
}
stats_normal$padj <- p.adjust(stats_normal$p, method = "BH")
view(stats_normal)

#stats for MCP deconvolution adeneoma polyps vs serrated polyps
mcp_counter_polyp <- mcp_counter_long |> filter(condition == "polyp") |> arrange(cell_type, subject)
stats_polyp <- tibble()
for (cell in unique(mcp_counter_polyp$cell_type)) {
  cell_data <- mcp_counter_polyp |> filter(cell_type == cell)
  wilcox_result <- wilcox_test(cell_estimate ~ histology, data = cell_data)
  wilcox_result <- wilcox_result |> mutate(cell_type = cell)
  stats_polyp <- bind_rows(stats_polyp, wilcox_result)
}
stats_polyp$padj <- p.adjust(stats_polyp$p, method = "BH")
view(stats_polyp)


#runnig single sample gene set enrichemnt analysis using GSVA pachages. This allows to generate enrichement score for each sample and to then investigate correlations between enrichemnt score and other scores of interest  ------
library(GSVA)
library(GSEABase)

#maing sure that the expression matrix is numeric and investigate the resulting matrix
norm_count_matrix <- as.matrix(serrated_matrix)
str(norm_count_matrix)

# Load gene sets for KEGG gene sets and hallmarks gene sets (previously downloaded form MSigDb)
hallmark <- getGmt("~/Desktop/serrated seq/RNAseq final/GMTs/h.all.v2024.1.Hs.symbols.gmt")
kegg <- getGmt("~/Desktop/serrated seq/RNAseq final/GMTs/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")
kegg_med <- getGmt("~/Desktop/serrated seq/RNAseq final/GMTs/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt")

# Combine all KEGG gene sets and hallmarks of cancer gene set
all_pathways <- c(hallmark, kegg, kegg_med)

#running GSVE to get per-sample enrichment scores 
gsvaPar_h <- gsvaParam(norm_count_matrix, hallmark)
ssgsea_h <- gsva(gsvaPar_h, verbose=FALSE)
str(ssgsea_h)

gsvaPar_kegg <- gsvaParam(norm_count_matrix, kegg)
ssgsea_kegg <- gsva(gsvaPar_kegg, verbose=FALSE)
str(ssgsea_kegg)

gsvaPar_kegg_med <- gsvaParam(norm_count_matrix, kegg_med)
ssgsea_kegg_med <- gsva(gsvaPar_kegg_med, verbose=FALSE)
str(ssgsea_kegg_med)

# Row-bind gene sets for enrichemnts in hallmarks of cancer and KEGG databases
ssgsea_combined <- rbind(ssgsea_h, ssgsea_kegg, ssgsea_kegg_med)

# Extract the pathways of interest 
pathways_of_interest <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                          "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                          "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
                          "KEGG_MEDICUS_REFERENCE_CD80_CD86_CTLA4_PP2A_SIGNALING_PATHWAY",
                          "KEGG_MEDICUS_REFERENCE_PDL_PD1_SHP_PI3K_SIGNALING_PATHWAY")
cor_df <- t(ssgsea_combined[pathways_of_interest , ])

#merging cor_df and MCP_serrated 
cor_df <- as.data.frame(cor_df)
cor_df$names <- rownames(cor_df)
cor_combined_df <- merge(cor_df, MCP_serrated, by = "names")
rownames(cor_combined_df) <- cor_combined_df$names
cor_combined_df$names <- NULL

# Correlation matrix formed 
cor_matrix <- cor(cor_combined_df, method = "pearson")
view(cor_matrix)
write.csv(cor_matrix, "ssgsea_cell_correlations.csv")

#repeating after subsetting the matricx for patients with adenoma and serrated polyps, separately
serrated_samples <- c("S22_4732_Nor", "S21_4732_Pol", "S24_4732_Nor", "S23_4732_Pol",
                      "S12_6344_Nor", "S11_6344_Pol", "S6_6411_Norm", "S5_6411_Poly",
                      "S18_6794_Nor", "S17_6794_Pol", "S16_6794_Nor", "S15_6794_Pol",
                      "S38_7713_Nor", "S37_7713_Pol", "S42_7984_Nor", "S41_7984_Nor",
                      "S40_7984_Pol", "S39_7984_Pol")

#subsetting samples by histology 
cor_df_serrated <- cor_combined_df[rownames(cor_df) %in% serrated_samples, ]
cor_df_adenoma <- cor_combined_df[!(rownames(cor_df) %in% serrated_samples), ]

# generation correlation matrices for each histology 
cor_mat_ser <- cor(cor_df_serrated, method = "pearson")
cor_mat_ade <- cor(cor_df_adenoma, method = "pearson")

#top get p values corresponding to each correlation, matrices with corresponding p values in place of pearson coefficients were generated. Here the p value matrices are initialized 
p_mat_ser <- matrix(NA, ncol = ncol(cor_df_serrated), nrow = ncol(cor_df_serrated))
colnames(p_mat_ser) <- rownames(p_mat_ser) <- colnames(cor_df_serrated)

p_mat_ade <- matrix(NA, ncol = ncol(cor_df_adenoma), nrow = ncol(cor_df_adenoma))
colnames(p_mat_ade) <- rownames(p_mat_ade) <- colnames(cor_df_adenoma)

#here, a loop is made to generate p values for each pearson coefficinet in the correlation matrix 
for (i in 1:ncol(cor_df_serrated)) {
  for (j in 1:ncol(cor_df_serrated)) {
    test <- cor.test(cor_df_serrated[, i], cor_df_serrated[, j])
    p_mat_ser[i, j] <- test$p.value
  }
}

for (i in 1:ncol(cor_df_adenoma)) {
  for (j in 1:ncol(cor_df_adenoma)) {
    test <- cor.test(cor_df_adenoma[, i], cor_df_adenoma[, j])
    p_mat_ade[i, j] <- test$p.value
  }
}

#Here I am converting the p-value matrices to character matrices and converting all p values above 0.05 to NS (for easier visualisation of significant correlation)
p_mat_ser_char <- ifelse(p_mat_ser > 0.05, "NS", formatC(p_mat_ser, format = "e", digits = 2))
p_mat_ade_char <- ifelse(p_mat_ade > 0.05, "NS", formatC(p_mat_ade, format = "e", digits = 2))

#save matrices with pearson coefficients and with corresponiding p values
write.csv(cor_matrix_serrated, "ssgsea_cells_correlations_serrated.csv")
write.csv(cor_matrix_adenoma, "ssgsea_cells_correlations_adenoma.csv")

write.csv(p_mat_ser_char, "pvalues_correlations_serrated.csv")
write.csv(p_mat_ade_char, "pvalues_cells_correlations_adenoma.csv")

#statistical tests for correlations of interest
ade_ifn_mhc <- cor.test(cor_df_adenoma$HALLMARK_INTERFERON_GAMMA_RESPONSE, cor_df_adenoma$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION)
ade_tnf_mhc <- cor.test()
