# install Bioconductor packages from BiocManager. ----------------------------------------------------
#NOTE: CRAN pacakges can be downloaded with install() function, but packages from bioconductor are instaled with Biocmanager::install() function 
library(clusterProfiler) #package nessessary for GSEA 
library(AnnotationDbi) #Provides a framework to access annotation databases. Required for org.Hs.eg.db
library(org.Hs.eg.db) #Maps gene IDs and provides gene metadata (chr position, GOs). species-specific
library(msigdbr)
library(biomaRt)
library(tidyverse)
library(dplyr)
library(ggplot2) # package for graphs 
library(ggbreak)
library(ggrepel) #pacakge for better spreadout labels 
library(cowplot)
library(enrichplot)
library(RColorBrewer) #package for colour palletes for graphs
library(pheatmap) #package for heatmaps 
library(magick) #package for image processing 

#load the csv files with results of the DESeq2 analysis - output of DESEQ2 includes log2 fold changes in gene expression
res_adenoma <- read_csv("~/Desktop/serrated seq/RNAseq final/adenoma_polyp_vs_normal_final_results.csv")
res_serrated <- read_csv("~/Desktop/serrated seq/RNAseq final/serrated_polyp_vs_normal_final_results.csv")
res_comparison <- read_csv("~/Desktop/serrated seq/RNAseq final/change_from_normal_to_ser_vs_ade_results.csv")
res_nvsn <- read_csv("~/Desktop/serrated seq/RNAseq final/normal_serrated_vs_adenoma_final_results.csv")
res_pvsp <- read_csv("~/Desktop/serrated seq/RNAseq final/polyp_serrated_vs_adenoma_final_results.csv")

#for GSEA analysis create a table including fold changes sorted in decreasing order and gene_ids. Make sure that the are no duplicates in gene names ---------
##for res_adenoma 
res_adenoma <- res_adenoma |> rename(gene_id = ...1) |> as.data.frame()
res_adenoma_list <- res_adenoma |> dplyr::select(gene_id, log2FoldChange) |> arrange(desc(log2FoldChange)) 
geneList_adenoma <- setNames(res_adenoma_list$log2FoldChange, res_adenoma_list$gene_id)

##for res_serrated 
res_serrated <- res_serrated |> rename(gene_id = ...1)
res_serrated_list <- res_serrated |> dplyr::select(gene_id, log2FoldChange) |> arrange(desc(log2FoldChange))
geneList_serrated <- setNames(res_serrated_list$log2FoldChange, res_serrated_list$gene_id)

##for res_comparison
res_comparison <- res_comparison |> rename(gene_id = ...1)
res_comparison_list <- res_comparison |> dplyr::select(log2FoldChange, gene_id) |> arrange(desc(log2FoldChange))
geneList_comparison <- setNames(res_comparison_list$log2FoldChange, res_comparison_list$gene_id)

##for res_nvsn
res_nvsn <- res_nvsn |> rename(gene_id = ...1)
res_nvsn_list <- res_nvsn |> dplyr::select(log2FoldChange, gene_id) |> arrange(desc(log2FoldChange))
geneList_nvsn <- setNames(res_nvsn_list$log2FoldChange, res_nvsn_list$gene_id)

##for res_pvsp
res_pvsp <- res_pvsp |> rename(gene_id = ...1)
res_pvsp_list <- res_pvsp |> dplyr::select(log2FoldChange, gene_id) |> arrange(desc(log2FoldChange))
geneList_pvsp <- setNames(res_pvsp_list$log2FoldChange, res_pvsp_list$gene_id)

#create gene enrichment results using MSigDb. From this database extract the category of gene set you are interested in. For example C5 includes gene ontology information, H includes clancer hallmarks gene sets. Extract the gene names from the same category that you have in gene lists. For currated gene sets C2 (includes kegg). C8 are cell types ---------
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") |>
  dplyr::select(gs_name, gene_symbol) |>
  mutate(gs_name = sub("^HALLMARK_", "", gs_name)) |>
  mutate(gs_name = gsub("_", " ", gs_name))
head(H_t2g)

C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") |> 
  dplyr::select(gs_name, gene_symbol) 

head(C2_t2g)
#filter KEGG pathways form currated signature C2 
KEGG_t2g <- C2_t2g |> 
  filter(grepl(pattern = "^KEGG", gs_name)) |>
  mutate(gs_name = sub("^KEGG_", "", gs_name)) |>
  mutate(gs_name = sub("^MEDICUS_REFERENCE_", "", gs_name)) |>
  mutate(gs_name = sub("^MEDICUS_", "", gs_name)) |>
  mutate(gs_name = gsub("_", " ", gs_name))

#comparisons in hallmarks ------
#ade vs normal
set.seed(97)
gse_hallmarks_adenoma <- GSEA(geneList_adenoma, TERM2GENE = H_t2g)
gse_hallmarks_adenoma_dotplot <- dotplot(gse_hallmarks_adenoma, x = "NES", color = "p.adjust", showCategory=25) + 
  theme_light() +
  lims(x=c(-3, 3)) +
  ggtitle("adenoma/ polyp vs. normal hallmarks")

#serrated vs normal
set.seed(256)
gse_hallmarks_serrated <- GSEA(geneList_serrated, TERM2GENE = H_t2g)
gse_hallmarks_serrated_dotplot <- dotplot(gse_hallmarks_serrated, x = "NES", showCategory=20) + 
  theme_light() +
  lims(x=c(-3, 3)) +
  ggtitle("serrated/ polyp vs. normal hallmarks") 

#normal to serrated vs polyp 
gse_hallmarks_comparison <- GSEA(geneList_comparison, TERM2GENE = H_t2g)
gse_hallmarks_comparison_dotplot <- dotplot(gse_hallmarks_comparison,  x = "NES", showCategory=25) + 
  ggtitle("normal to serrated vs. adenoma hallmarks") +
  theme_light() +
  lims(x=c(-3, 3))

#shared doplot for serrated and adenomas 
dfH_adenoma <- as.data.frame(gse_hallmarks_adenoma_dotplot[["data"]])
dfH_serrated <- as.data.frame(gse_hallmarks_serrated_dotplot[["data"]])
dfH_comparison <- as.data.frame(gse_hallmarks_comparison_dotplot[["data"]])
dfH_adenoma$group <- "Adenoma vs normal"
dfH_serrated$group <- "Serrated vs normal"
dfH_comparison$group <- "Serrated vs adenoma"
combined_gsea <- rbind(dfH_adenoma, dfH_serrated, dfH_comparison)
combined_gsea_sig <- combined_gsea |>
  filter(p.adjust < 0.05) 
combined_gsea_sig$group <- as.factor(combined_gsea_sig$group)
combined_gsea_sig$group <- relevel(combined_gsea_sig$group, ref = "Serrated vs normal")

ggplot(combined_gsea_sig, aes(x = group, y = Description, size = GeneRatio, color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0 ) +
  theme_minimal() +
  labs(title = "Cancer hallmarks GSEA",
       x = "Comparison",
       y = "Hallmark Pathway",
       size = "GeneRatio",
       color = "NES")

combined_gsea_comparison <- combined_gsea_sig |> filter(combined_gsea_sig$group == "Serrated vs adenoma")
ggplot(combined_gsea_comparison, aes(x = NES, y = Description, size = -p.adjust, colour = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",  midpoint = 0 ) +
  theme_minimal(base_size = 12) +
  labs(title = "Cancer hallmarks serrated vs adenoma",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       color = "NES",
       size = "Adjusted p-value") +
  theme(strip.text = element_text(face = "bold"))


#comparison in KEGG --------
#adenoma vs normal
gse_KEGG_adenoma <- GSEA(geneList_adenoma, TERM2GENE = KEGG_t2g)
dotplot(gse_KEGG_adenoma, x = "NES", color = "p.adjust", showCategory=20) + 
  theme_light() +
  lims(x=c(-3, 3)) +
  ggtitle("adenoma/ polyp vs. normal KEGG") 

#serrated vs normal
set.seed(256)
gse_KEGG_serrated <- GSEA(geneList_serrated, TERM2GENE = KEGG_t2g)
dotplot(gse_KEGG_serrated, x = "NES", showCategory=20) + 
  theme_light() +
  lims(x=c(-3, 3)) +
  ggtitle("serrated/ polyp vs. normal KEGG") 

#normal to serrated vs polyp 
gse_KEGG_comparison <- GSEA(geneList_comparison, TERM2GENE = KEGG_t2g)
dotplot(gse_KEGG_comparison,  x = "NES", showCategory=20) + 
  ggtitle("normal to serrated vs. adenoma hallmarks") +
  theme_light() +
  lims(x=c(-3, 3))

#normal serrated vs adenoma 
gse_KEGG_nvsn <- GSEA(geneList_nvsn, TERM2GENE = KEGG_t2g)
dotplot(gse_KEGG_nvsn,  x = "NES", showCategory=20) + 
  ggtitle("normal in serrated vs. adenoma KEGG") +
  theme_light() +
  lims(x=c(-3, 3))

#polyp serrated vs adenoma 
gse_KEGG_pvsp <- GSEA(geneList_pvsp, TERM2GENE = KEGG_t2g)
dotplot(gse_KEGG_pvsp,  x = "NES", showCategory=20) + 
  ggtitle("serrated vs. adenoma polyp KEGG") +
  theme_light() +
  lims(x=c(-3, 3))

#shared doplot for serrated and adenomas 
dfKEGG_adenoma <- as.data.frame(gse_KEGG_adenoma@result)
dfKEGG_serrated <- as.data.frame(gse_KEGG_serrated@result)
dfKEGG_comparison <- as.data.frame(gse_KEGG_comparison@result)
dfKEGG_adenoma$group <- "Adenoma vs normal"
dfKEGG_serrated$group <- "Serrated vs normal"
dfKEGG_comparison$group <- "Serrated vs adenoma"
combined_KEGG <- rbind(dfKEGG_adenoma, dfKEGG_serrated, dfKEGG_comparison)
combined_KEGG_sig <- combined_KEGG |>
  filter(p.adjust < 0.05)
combined_KEGG_sig$group <- as.factor(combined_KEGG_sig$group)
combined_KEGG_sig$group <- relevel(combined_KEGG_sig$group, ref = "Serrated vs normal")

ggplot(combined_KEGG_sig, aes(x = group, y = Description, size = -p.adjust, color = NES)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
  theme_minimal() +
  labs(title = "KEGG Adenoma vs Serrated",
       x = "Group",
       y = "KEGG Pathway",
       size = "Adjusted p-value",
       color = "NES")

combined_KEGG_comparison <- combined_KEGG_sig |> filter(combined_KEGG_sig$group == "Serrated vs adenoma")
ggplot(combined_KEGG_comparison, aes(x = NES, y = Description, color = p.adjust, size = abs(NES))) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue", trans = "reverse") +
  theme_minimal(base_size = 12) +
  labs(title = "Serrated vs adenoma GSEA",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       color = "Adjusted p-value",
       size = "Count") +
  theme(strip.text = element_text(face = "bold"))

#verifying GSEA for antigen presentation --------
antigen_t2g <- KEGG_t2g |>
  filter(grepl(pattern = "^ANTIGEN PROCESSING AND PRESENTATION$", gs_name))

#GSEA for antigen processing for adenomas
hla_adenoma <- GSEA(geneList_adenoma, 
                    TERM2GENE = antigen_t2g, 
                    minGSSize = 1,
                    pvalueCutoff = 1)
gseaplot2(
  hla_adenoma,
  geneSetID = 1,
  color = c("#E495A5", "#7DB0DD"),
)

hla_adenoma_res <- dotplot(hla_adenoma)

#for serrated lesions 
hla_serrated <- GSEA(geneList_serrated, 
                    TERM2GENE = antigen_t2g, 
                    minGSSize = 1,
                    pvalueCutoff = 1)
gseaplot2(
  hla_serrated,
  geneSetID = 1,
  color = c("#E495A5", "#7DB0DD"),
)
hla_serrated_res <- dotplot(hla_serrated)

#for serrated vs adenoma lesions 
hla_ser_vs_ade <- GSEA(geneList_comparison, 
                       TERM2GENE = antigen_t2g, 
                       minGSSize = 1,
                       pvalueCutoff = 1)
gseaplot2(
  hla_ser_vs_ade,
  geneSetID = 1,
  color = c("#E495A5", "#7DB0DD"),
)

hla_comparison_res <- dotplot(hla_ser_vs_ade)

#making a single graph for hla 
hla_adenoma <- as.data.frame(hla_adenoma_res[["data"]])
hla_serrated <- as.data.frame(hla_serrated_res[["data"]])
hla_ser_vs_ade <- as.data.frame(hla_comparison_res[["data"]])
hla_adenoma$group <- "Adenoma vs normal"
hla_serrated$group <- "Serrated vs normal"
hla_ser_vs_ade$group <- "Serrated vs adenoma"
hla <- rbind(hla_adenoma, hla_serrated, hla_ser_vs_ade)
hla$group <- as.factor(hla$group)
hla$group <- relevel(hla$group, ref = "Serrated vs normal")

ggplot(hla, aes(x = NES, y = group, size = GeneRatio, color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white",  high = "red",  midpoint = 0 ) +
  theme_minimal(base_size = 12) +
  labs(
    x     = "NES for Antigen Processing and Presentation",
    y     = NULL,
    size  = "Gene Ratio",
    color = "NES") +
  theme(strip.text = element_text(face = "bold"))

ggplot(hla, aes(x = group, y = Description, size = GeneRatio, color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white",  high = "red",  midpoint = 0 ) +
  theme_minimal(base_size = 12) +
  labs(
    x     = "Comparison",
    y     = NULL,
    size  = "Gene Ratio",
    color = "NES") +
  theme(strip.text = element_text(face = "bold"))

#verifying GSEA for PDL1/PD and CTLA4 -------
icb_t2g <- KEGG_t2g |>
  filter(gs_name == "CD80 CD86 CTLA4 PP2A SIGNALING PATHWAY" | gs_name == "PDL PD1 SHP PI3K SIGNALING PATHWAY")

#GSEA for PDL1/PD and CTLA4  for adenomas
icb_adenoma <- GSEA(geneList_adenoma, 
                    TERM2GENE = icb_t2g, 
                    minGSSize = 1,
                    pvalueCutoff = 1)
gseaplot2(
  icb_adenoma,
  geneSetID = 1:2,
  color = c("#E495A5", "#7DB0DD"),
)

icb_adenoma_res <- dotplot(icb_adenoma)

#PDL1/PD and CTLA4 enrichemnt for serrated lesions 
icb_serrated <- GSEA(geneList_serrated, 
                     TERM2GENE = icb_t2g, 
                     minGSSize = 1,
                     pvalueCutoff = 1)
gseaplot2(
  icb_serrated,
  geneSetID = 1:2,
  color = c("#E495A5", "#7DB0DD"),
)

icb_serrated_res <- dotplot(icb_serrated)

#PDL1/PD and CTLA4 enrichment for serrated vs adenoma lesions 
icb_ser_vs_ade <- GSEA(geneList_comparison, 
                       TERM2GENE = icb_t2g, 
                       minGSSize = 1,
                       pvalueCutoff = 1)
gseaplot2(
  icb_ser_vs_ade,
  geneSetID = 1:2,
  color = c("#E495A5", "#7DB0DD"),
)

icb_comparison_res <- dotplot(icb_ser_vs_ade)

#making a single graph for hla 
icb_adenoma <- as.data.frame(icb_adenoma_res[["data"]])
icb_serrated <- as.data.frame(icb_serrated_res[["data"]])
icb_comparison <- as.data.frame(icb_comparison_res[["data"]])
icb_adenoma$group <- "Adenoma vs normal"
icb_serrated$group <- "Serrated vs normal"
icb_comparison$group <- "Serrated vs adenoma"
icb <- rbind(icb_adenoma, icb_serrated, icb_comparison)
icb$group <- as.factor(icb$group)
icb$group <- relevel(icb$group, ref = "Serrated vs normal")

ggplot(icb, aes(x = NES, y = group, size = GeneRatio, color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white",  high = "red",  midpoint = 0 ) +
  theme_minimal(base_size = 12) +
  labs(
    x     = "NES",
    y     = NULL,
    size  = "Gene Ratio",
    color = "NES") +
  theme(strip.text = element_text(face = "bold")) +
  facet_wrap(~ID, nrow = 2)

ggplot(icb, aes(x = group, y = Description, size = GeneRatio, color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white",  high = "red",  midpoint = 0 ) +
  theme_minimal(base_size = 12) +
  labs(
    x     = "Comparison",
    y     = NULL,
    size  = "Gene Ratio",
    color = "NES") +
  theme(strip.text = element_text(face = "bold"))
