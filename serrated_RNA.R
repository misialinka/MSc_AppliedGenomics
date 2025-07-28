# install Bioconductor packages from BiocManager. ----------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(tidyverse) #package for dataanalysis using ggplot2 and dplyr packages
library(tidyr) #picoting and nesting datasets 
library(ggplot2)
library(ggpubr) #costumizing ggplot plots 
library(ggrepel) #pacakge for better spreadout labels (geom_text_repel() and geom_label_repel())
library("EnhancedVolcano")
library(tximeta) #this package is needed for generating SummarizedExperiment file.
library(DESeq2) # package for DEG analysis 
library(rstatix) #package for statistical analysis 
library(ashr) #package for shirnking DESEQ results 
library(biomaRt) #changing a gene naming system 
library(sva) #package for removing batch effects
library(RColorBrewer)

#loading metadata (dataframe containing clinical about samples included in the analysis) 
metadata_serrated <- as.data.frame(read_csv("~/Desktop/serrated seq/RNAseq final/RNAmeta_serrated.csv"))

#cleaning the data - removing rows with samples that have no RNAseq samples available (patient ID 6341) or have undefined histology (946, 6894) 
rows_to_remove <- which(metadata_serrated$`Patient ID` == 946 | 
                          metadata_serrated$`Patient ID` == 6341 | 
                          metadata_serrated$`Patient ID` == 6894)

metadata_serrated <- metadata_serrated[-rows_to_remove, ]
view(metadata_serrated)

#selecting all clinical features needed for analysis from metadata - preparing information for the SummarizedExperiment object
coldata_serrated <- data.frame(
  names = metadata_serrated$`RNA Name`,  #sample names - needed to identify samples
  condition = metadata_serrated$samples,  #condition specifies if the samples comes form "polyp" or "normal" samples
  histology = metadata_serrated$histology, #histology defines if patient has serrated or adenoma polyps
  subject = metadata_serrated$`pair ID`,  #subject defines normal-polyp pairs from the same patient
  batch = metadata_serrated$`Sample batch`,  #defines from which batch of sequjencing the samples comes from 
  side = metadata_serrated$`Location side`,  #specifies "left" or "right" side of the colon, from which sample was taken 
  allhistologies = metadata_serrated$`histology all`, #includes categories "normal" "adenoma" "serrated"  
  stringsAsFactors = FALSE
)

coldata_serrated$histology <- as.factor(coldata_serrated$histology) 
coldata_serrated$subject <- as.factor(coldata_serrated$subject) 
coldata_serrated$condition <- as.factor(coldata_serrated$condition)
coldata_serrated$batch <- as.factor(coldata_serrated$batch)

#assigning sample names as rownames of coldata 
row.names(coldata_serrated) <- coldata_serrated$names

#loading file with expression matrix including gene counts from 1st and 2nd batch of sequencing (make sure that gene IDs are rows and sample names are columns)
gene_counts_serrated1 <- as.matrix(read.csv("~/Desktop/serrated seq/RNAseq final/counts_part1.csv", row.names = 1)) 

#removing row counts for patients with undefined histology (patient ID 946, 6894)
gene_counts_serrated1 <- gene_counts_serrated1[ , !(colnames(gene_counts_serrated1) %in% c("S36_946_Norm", "S35_946_Poly", "S30_6894_Nor", "S32_6894_Nor", "S31_6894_Pol", "S29_6894_Pol"))]

#ensuring that matrix contains integers (required for DESEQ function)
gene_counts_serrated1 <- round(gene_counts_serrated1)
str(gene_counts_serrated1)   #57773 genes present in the matrix

#loading file with with expression matrix including gene counts from the last batch of sequencing. I am using dataframe as an intermediate to prevent matrix switching from numeric to character
gene_counts_serrated2 <- as.data.frame(read_csv("~/Desktop/serrated seq/RNAseq final/counts_part2.csv"))
rownames(gene_counts_serrated2) <- gene_counts_serrated2$gene_id
gene_counts_serrated2 <- gene_counts_serrated2[ , !(colnames(gene_counts_serrated2) %in% c("gene_name", "gene_id"))]
gene_counts_serrated2 <- as.matrix(gene_counts_serrated2)
gene_counts_serrated2 <- round(gene_counts_serrated2)
str(gene_counts_serrated2)     #29744 genes

#CONVERTING ENSEMBL IDS TO HUGO. Gene names in the 1st batch are ENSEMBL IDs, while in the second are HUGO symbols. To merge the two batches into a single matrix and run DESEQ function on the batches together, the names need to be converted to the same type using biomart.  -----------------------------------------------------
#creating mart object for ENSEMBL IDs  
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#extracting all ENSEMBL gene names in the expression matrix from the 1st batch of samples
ensembl_genes <- rownames(gene_counts_serrated1)

#inspect number of genes in hugo_genes list. Here, there are 57773 genes
str(ensembl_genes)

#creating dataframe with HUGO names and matched ENSEMBL IDs 
conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = ensembl_genes,
  mart = ensembl
)

#loading an intrinsically generated file with ENSEMBL IDs and HUGO symbol gene names for conversion. This file includes older HUGO system (e.g., DEC1 gene symbol instead of DELEC1), which was used in the annotation for the last batch of RNA seq data. 
gene_names <- as.data.frame(read_csv("~/Desktop/serrated seq/RNAseq final/gene.names.csv"))

#merging external database with internal naming system, to help with more accurate mapping with the older HUGO symbols naming
gene_names <- gene_names |> 
  dplyr::rename(ensembl_gene_id = gene_id) 
conversion <- full_join(conversion, gene_names, by = "ensembl_gene_id") #object contains 57774 genes

#filtering protein coding genes to ensure that the most relevant match is kept when duplicated names will be removed (later int the code)
conversion <- conversion |>
  dplyr::filter(gene_biotype == "protein_coding") |>
  dplyr::filter(hgnc_symbol != "")

#inspect number of genes in conversion list. Here, there are 18926 genes
str(conversion)

#removing duplicates. NOTE: removing duplicates is necessary for further analysis, but may accidentaly remove relevant ENSEMBL IDs. 
conversion <- conversion[!duplicated(conversion$ensembl_gene_id), ]

#selecting genes in the 1st batch that map to the selected protein-coding genes
mapped_genes <- ensembl_genes[ensembl_genes %in% conversion$ensembl_gene_id]
gene_counts_serrated1 <- gene_counts_serrated1[mapped_genes, ]

#adding the new gene names as rownames
hugo_ids <- conversion$gene_name[match(rownames(gene_counts_serrated1), conversion$ensembl_gene_id)]
rownames(gene_counts_serrated1) <- hugo_ids

#merging the two files -----------------------------
#adding a name of a column with the gene names, to allow merging by the gene names. 
gene_counts_serrated1_t <- as_tibble(gene_counts_serrated1, rownames = "hugo_id")
gene_counts_serrated2_t <- as_tibble(gene_counts_serrated2, rownames = "hugo_id")

#merging the two files with matching by the gene names. 30117 genes here. 
serrated_counts <- full_join(gene_counts_serrated1_t, gene_counts_serrated2_t, by = "hugo_id")

#inspect how many NA values are found in the matrix. Here there are 406526. These could be non-protein coding genes present in the 2nd, unfiltered batch e.g., SEPT5-GP1BB (NAs are present in the added batch and not the protein-selected batch)
sum(is.na(serrated_counts))

#filter all rows with gene counts with value of 0 and NA values. 17705 left 
serrated_counts <- serrated_counts %>%
  filter(rowSums(dplyr::select(., -hugo_id) != 0) > 0)
serrated_counts <- serrated_counts |> drop_na()

#investigate any duplicates that are still present 
serrated_counts$hugo_id[duplicated(serrated_counts$hugo_id)] #NPIPA7 ZNF286A RPP14 SOX7 FAM47E-STBD1 TM4SF19 

#duplicates need to be removed for DESEQ. 
serrated_counts <- serrated_counts %>%
  distinct(hugo_id, .keep_all = TRUE) 

#assigning gene names to rownames 
rownames(serrated_counts) <- serrated_counts$hugo_id

#converting back to expression matrix
serrated_counts <- as.matrix(serrated_counts)

#removing the row "esembl_id" as it was only needed for matching the two files for merging
serrated_counts <- serrated_counts[, -1]
serrated_counts_numeric <- apply(serrated_counts, 2, as.numeric)
rownames(serrated_counts_numeric) <- rownames(serrated_counts)

#investigate if the expression matrix looks correct 
str(serrated_counts_numeric)  

##pre-filtering to exclude genes that are poorly detected, based on number of groups---------
#see what is the number of samples in your smallest group (following recommendation by the DESEQ2 author M. Love (https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) ) 
coldata_serrated |> group_by(histology, condition) |> summarise(n()) 
smallestGroupSize <- 9

#keeping rows in which at least 9 samples have at least 10 counts. 15087 genes remaining
keep <- rowSums(serrated_counts_numeric >= 10) >= smallestGroupSize
serrated_counts_numeric <- serrated_counts_numeric[keep,]

#make sure that rownames in coldata and colnames in gene_counts_MSI are the same. This is crucial for the SummarizedExperiemnt object to be able to link expression matrix to the samples ------------------
all(rownames(coldata_serrated) %in% colnames(serrated_counts_numeric))

#make sure that the order of rownames in coldata and colnames in gene_counts_MSI are the same. 
all(rownames(coldata_serrated) == colnames(serrated_counts_numeric))

#If not, reorder them to match and assing names of rows and column correctly
coldata_serrated <- coldata_serrated |>
  arrange(desc(histology), subject, condition)
serrated_counts_numeric <- serrated_counts_numeric[ ,coldata_serrated$names]

#adjusting for the batch effect, since PCA has previously shown strong batch effect--------
##assigning which samples in the matrix come from different batches
batch <- c(rep(1, 18), rep(2, 10), rep(1, 18)) 

##assigning matrices with biological effects of interest (covariates). These lines were used for testing of which model has best batch correction based on the PCA.
cov_condition <- c(rep(c(0, 1), 23))
cov_histology <- c(rep(0, 18), rep(1, 28))
covar_mat <- cbind(cov_condition, cov_histology)

##running batch correction. To look at PCA for different methods of batch corrections I have also tested "group = cov_condition", to adjust for the polyp vs. normal variation and "covar_mod=covar_mat" to adjust for both the histology and polyp vsn normal variation. These additions were not improving the batch effect in the PCA plot. 
serrated_counts_adjusted <- ComBat_seq(serrated_counts_numeric, batch = batch, group = NULL) 

#saving the merged, batch-adjusted matrix for deconvolution analysis 
write.csv(as.data.frame(serrated_counts_adjusted), 
          file="batchadjusted_serrated_expression.csv") 

# Create the SummarizedExperiment object with raw counts and coldata------------------------
se_serrated <- SummarizedExperiment(
  assays = list(counts = serrated_counts_adjusted),
  colData = coldata_serrated
)

# To account for hypermutated status and tissue type in a paired fashion and to avoid redundancy error hindering DESeq paired analysis, custom model matrix must be provided. -----------------------
tmp <- colData(se_serrated) |> 
  as_tibble() |> #tibble for easier manipulation
  arrange(desc(histology), subject, condition) |> 
  mutate(tmp_id=c(rep(seq(1, 9), each=2), rep(seq(1, 14), each=2))) |>
  dplyr::select(names, histology, tmp_id, condition, subject)

# Sort in the same way as colData(se). These steps arrange tmp_id which was created in the previous step, select names, arrange tmp by names and save only tmp_id from the table. 
tmp_id <- colData(se_serrated) |> 
  as.tibble() |> 
  dplyr::select(names) |> #
  left_join(tmp, by='names') |> 
  pull(tmp_id)

# Once tmp_id has been made and sorted in the correct order, it can be added
colData(se_serrated)$tmp_id <- as.factor(tmp_id)

#view the colData with the new tmp_id column to ensure that the order is correct 
view(colData(se_serrated))

#create model matrix for a design, following recommendation by the DESEQ2 author M. Love (https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)  -----------------------
mm <- model.matrix(~ histology + histology:tmp_id + histology:condition, colData(se_serrated))
unname(mm)
all.zero <- apply(mm, 2, function(x) all(x==0))
all.zero
idx <- which(all.zero)
mm <- mm[,-idx]
unname(mm)

#running DESEQ   -----------------------
#generation of a DESeqDataSet from a RangedSummarizedExperiment
dds_serrated <- DESeqDataSet(se_serrated, design = ~ subject + condition)

#sets "normal" as the reference for comparison between conditions
dds_serrated$condition <- relevel(dds_serrated$condition, ref = "normal") 

#to ensure if your experimental groups are read as levels. If not present, then it is neccessary to add them manually
levels(dds_serrated$condition)
levels(dds_serrated$histology)

#derive log2 fold change and Wald test p value by running DESeq(). NOTE: DESeq() automatically normalized
dds_serrated <- DESeq(dds_serrated, full = mm)
DEG_serrated <- results(dds_serrated)
resultsNames(dds_serrated)

#extract results relevant for my work, which is adenoma vs. normal, adenoma vs. serrated, serrated vs. normal  -----------------------
res_adenoma <- results(dds_serrated, alpha=0.05, name='histologyadenoma.conditionpolyp') #paired adenoma vs. normal
res_serrated <- results(dds_serrated, alpha = 0.05, name = 'histologyserrated.conditionpolyp') #paired serrated vs. normal

#extract comparison of changes from normal to adenoma polyp vs serrated polyp 
contrast_list <- list(c("histologyserrated.conditionpolyp"), c("histologyadenoma.conditionpolyp"))
res_comparison <- results(dds_serrated, contrast = contrast_list)

#export results as csv for later GSEA (seperate script for GSEA and cell deconvolution)
write.csv(as.data.frame(res_adenoma), 
          file="adenoma_polyp_vs_normal_final_results.csv")
write.csv(as.data.frame(res_serrated), 
          file="serrated_polyp_vs_normal_final_results.csv")
write.csv(as.data.frame(res_comparison),
          file = "change_from_normal_to_ser_vs_ade_results.csv")

#plot PCA for visual control of batch effect ---------
vsd_serrated <- vst(dds_serrated,blind=FALSE)

pcaData <- plotPCA(vsd_serrated, intgroup=c("allhistologies", "batch", "side", "tmp_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, shape=allhistologies)) +
  geom_point(size=4, aes(color = batch)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  stat_ellipse(aes(linetype = allhistologies), color = "gray60") +
  labs(title="PCA without ComBat-seq correction")

#Conduct log fold change shrinkage to improve visualization and ranking. To see which conditions to compare you can use resultsNames(dds) function.   -----------------------
res_adenoma_shrunk <- lfcShrink(dds_serrated, coef = 'histologyadenoma.conditionpolyp', type = "ashr") 
res_serrated_shrunk <- lfcShrink(dds_serrated, coef = 'histologyserrated.conditionpolyp', type = "ashr")
res_comparison_shrunk <- lfcShrink(dds_serrated, contrast = contrast_list, type = "ashr")

# Prettify results for plotting
res_adenoma_prettified <- res_adenoma_shrunk |> 
  as_tibble(rownames='gene_names') |> 
  mutate(gene_name=DEG_serrated$gene_name, .before=2)

res_serrated_prettified <- res_serrated_shrunk |> 
  as_tibble(rownames='gene_names') |> 
  mutate(gene_name=DEG_serrated$gene_name, .before=2)

res_comparison_prettified <- res_comparison_shrunk |>
  as_tibble(rownames='gene_names') |>
  mutate(gene_name=DEG_serrated$gene_name, .before=2)

#extracting normalized counts and creating counts per million ------------
normalized_counts <- counts(dds_serrated, normalized=TRUE)
write.csv(as.data.frame(normalized_counts), 
          file="normalized_counts_matrix.csv")

#creating counts per million for visualisation
serratedls_cpm <- t(t(normalized_counts)/colSums(normalized_counts))*1E6

#investigating HLA expression   -----------------------
#creating my own signature for HLA class I and II 
hla1_genes <- c(
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "TAP1",
    "TAP2",
    "IRF1"
    )

hla2_genes <- c(
      "CIITA",
      "HLA-DPA1",
      "HLA-DPB1",
      "HLA-DQA1",
      "HLA-DQA2",
      "HLA-DQB1",
      "HLA-DQB2",
      "HLA-DRA",
      "HLA-DMA",
      "HLA-DMB",
      "HLA-DOA",
      "HLA-DOB"
    )

hla_genes <- unique(c(hla1_genes, hla2_genes))

#creating a tibble with all gene names, gene IDs, fold changes and p values for HLA genes for all 3 comparisons 
hla_genes_adenoma <- tibble(gene_names=hla_genes) |> left_join(res_adenoma_prettified, by='gene_names') |> arrange(padj)
hla_genes_serrated <- tibble(gene_names=hla_genes) |> left_join(res_serrated_prettified, by='gene_names') |> arrange(padj)
hla_genes_comparison <- tibble(gene_names=hla_genes) |> left_join(res_comparison_prettified, by='gene_names') |> arrange(padj)

#selecting HLA genes and the counts per million trasncirpt count matrices 
serratedls_cpm_HLA1 <- t(serratedls_cpm[hla1_genes, ])
serratedls_cpm_HLA2 <- t(serratedls_cpm[hla2_genes, ])

#adding a column with names of the samples (which was previously colnames), as they will be needed for merging datasets
serratedls_cpm_HLA1 <- as.data.frame(serratedls_cpm_HLA1)
serratedls_cpm_HLA1$names <- rownames(serratedls_cpm_HLA1)

serratedls_cpm_HLA2 <- as.data.frame(serratedls_cpm_HLA2)
serratedls_cpm_HLA2$names <- rownames(serratedls_cpm_HLA2)

#saving coldata infromation for merging with cpm (counts per million)
HLA_coldata <- as.data.frame(colData(dds_serrated))

#merging cpm data and metadata information
serratedls_cpm_HLA1 <- left_join(serratedls_cpm_HLA1, HLA_coldata, by = 'names')
serratedls_cpm_HLA2 <- left_join(serratedls_cpm_HLA2, HLA_coldata, by = 'names')

#converting to long format for visualisation
all_hla1_data <- serratedls_cpm_HLA1 |> pivot_longer(-c(names, condition, histology, subject, tmp_id, batch, allhistologies, side, sizeFactor), names_to='gene', values_to='CPM')

all_hla2_data <- serratedls_cpm_HLA2 |> pivot_longer(-c(names, condition, histology, subject, tmp_id, batch, allhistologies, side, sizeFactor), names_to='gene', values_to='CPM')

#generating a boxplots and stats for HLA expression -------
ggplot(all_hla1_data, aes(x=histology, y=CPM, fill=condition)) +
  geom_boxplot(outlier.shape=NA) + # Boxplot with transparency
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.7) +
  scale_fill_brewer(palette="Set2") +
  facet_wrap(~gene, scales = "free_y") +  # One boxplot per gene
  labs(title="Expression of HLA Genes",
       y="Normalized CPM") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) +
  theme_bw()   

ggplot(all_hla2_data, aes(x=histology, y=CPM, fill=condition)) +
  geom_boxplot(outlier.shape=NA) + # Boxplot with transparency
  scale_fill_brewer(palette="Set2") +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.7) +
  facet_wrap(~gene, scales = "free_y") +  # One boxplot per gene
  labs(title="Expression of HLA Genes",
       y="Normalized CPM") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) +
  theme_bw()  

#runnign statistics for all HLA gene counts comparisons, using Wilcox test with BH correction for multiple testing
#comparing serrated polyp vs paired normal expression of HLA class I 
all_hla1_serrated <- all_hla1_data |> filter(histology == "serrated") |> arrange(gene, tmp_id)
hla1_stats_ser <- tibble()
for (marker in unique(all_hla1_serrated$gene)) {
  genes_data <- all_hla1_serrated |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = genes_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  hla1_stats_ser <- bind_rows(hla1_stats_ser, wilcox_result)
}

hla1_stats_ser$padj <- p.adjust(hla1_stats_ser$p, method = "BH")
view(hla1_stats_ser)
  
#comparing serrated polyp vs paired normal expression of HLA class II
all_hla2_serrated <- all_hla2_data |> filter(histology == "serrated") |> arrange(gene, tmp_id)
hla2_stats_ser <- tibble()
for (marker in unique(all_hla2_serrated$gene)) {
  genes_data <- all_hla2_serrated |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = genes_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  hla2_stats_ser <- bind_rows(hla2_stats_ser, wilcox_result)
}
hla2_stats_ser$padj <- p.adjust(hla2_stats_ser$p, method = "BH")
view(hla2_stats_ser)

#comparing adenoma polyp vs paired normal expression of HLA class I
all_hla1_adenoma <- all_hla1_data |> filter(histology == "adenoma") |> arrange(gene, tmp_id)
hla1_stats_ade <- tibble()
for (marker in unique(all_hla1_adenoma$gene)) {
  genes_data <- all_hla1_adenoma |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = genes_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  hla1_stats_ade <- bind_rows(hla1_stats_ade, wilcox_result)
}
hla1_stats_ade$padj <- p.adjust(hla1_stats_ade$p, method = "BH")
view(hla1_stats_ade)

#comparing adenoma polyp vs paired normal expression of HLA class II
all_hla2_adenoma <- all_hla2_data |> filter(histology == "adenoma") |> arrange(gene, tmp_id)
hla2_stats_ade <- tibble()
for (marker in unique(all_hla2_adenoma$gene)) {
  genes_data <- all_hla2_adenoma |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = genes_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  hla2_stats_ade <- bind_rows(hla2_stats_ade, wilcox_result)
}
hla2_stats_ade$padj <- p.adjust(hla2_stats_ade$p, method = "BH")
view(hla2_stats_ade)

#generating a boxplots and stats for cell marker expression -------
#extracting main cell markers and comparing their expression across categories from serratedls_cpm (a file containing all cell counts )
cell_markers <- res_adenoma_prettified |> filter(gene_names=="CD4" | #Th cells
                                                   gene_names=="CD8A" | #Tcyt cells
                                                   gene_names == "FOXP3" | #Treg cells
                                                   gene_names == "CD68" | #macrophages
                                                   gene_names=="MS4A1" | #CD20 B cell marker
                                                   gene_names == "NCAM1"  #CD56 gene, NK cells
) |> 
  pull(gene_names)

serratedls_cpm_cellmarkers <- t(serratedls_cpm[cell_markers, ]) 
serratedls_cpm_cellmarkers <- as.data.frame(serratedls_cpm_cellmarkers)
serratedls_cpm_cellmarkers$names <- rownames(serratedls_cpm_cellmarkers)
cells_coldata <- as.data.frame(colData(dds_serrated))
serratedls_cpm_cells <- left_join(serratedls_cpm_cellmarkers, cells_coldata, by = 'names')
allls_cells <- serratedls_cpm_cells |> 
  pivot_longer(-c(names, condition, histology, subject, tmp_id, side, batch, sizeFactor, allhistologies), names_to='gene', values_to='CPM')

#plotting cell markers
ggplot(allls_cells, aes(x=histology, y=CPM, fill=condition)) +
  geom_boxplot(outlier.shape=NA) + # Boxplot with transparency
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.7) +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  facet_wrap(~gene, scales = "free_y") +  # One boxplot per gene
  labs(title="Expression of immune cell markers",
       y="Normalized CPM") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) 


#comparing adenoma polyp vs paired normal expression of immune cell marker genes
allls_cells_adenoma <- allls_cells |> filter(histology == "adenoma") |> arrange(gene, tmp_id)
cells_stats_ade <- tibble()
for (marker in unique(allls_cells_adenoma$gene)) {
  markers_data <- allls_cells_adenoma |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = markers_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  cells_stats_ade <- bind_rows(cells_stats_ade, wilcox_result)
}
cells_stats_ade$padj <- p.adjust(cells_stats_ade$p, method = "BH")
view(cells_stats_ade)

#comparing serrated polyp vs paired normal expression of immune cell marker genes
allls_cells_serrated <- allls_cells |> filter(histology == "serrated") |> arrange(gene, tmp_id)
cells_stats_ser <- tibble()
for (marker in unique(allls_cells_serrated$gene)) {
  markers_data <- allls_cells_serrated |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = markers_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  cells_stats_ser <- bind_rows(cells_stats_ser, wilcox_result)
}
cells_stats_ser$padj <- p.adjust(cells_stats_ser$p, method = "BH")
view(cells_stats_ser)

#checking expression of immune checkpoints ------
##extracting main immune checkpoints and comparing their expression across categories from serratedls_cpm (a file containing all cell counts )
icb <- res_adenoma_prettified |> filter(gene_names=="CD274" | #PDL1 
                                                   gene_names=="CTLA4" | 
                                                   gene_names == "PDCD1" | #PD1
                                                   gene_names == "HAVCR2" |#for T cell exhaustion TIM3
                                                   gene_names == "TIGIT" |
                                                   gene_names == "LAG3"
                                        
                                                  
) |> 
  pull(gene_names)

serratedls_cpm_icb <- t(serratedls_cpm[icb, ]) 
serratedls_cpm_icb <- as.data.frame(serratedls_cpm_icb)
serratedls_cpm_icb$names <- rownames(serratedls_cpm_icb)
icb_coldata <- as.data.frame(colData(dds_serrated))
serratedls_cpm_icb <- left_join(serratedls_cpm_icb, icb_coldata, by = 'names')
allls_icb <- serratedls_cpm_icb |> 
  pivot_longer(-c(names, condition, histology, subject, tmp_id, side, batch, sizeFactor, allhistologies), names_to='gene', values_to='CPM')

#plotting immune checkpoin inhibitor expression across histologies and normal
ggplot(allls_icb, aes(x=histology, y=CPM, fill=condition)) +
  geom_boxplot(outlier.shape=NA) + # Boxplot with transparency
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), alpha = 0.7) +
  scale_fill_brewer(palette="Set2") +
  theme_bw() +
  facet_wrap(~gene, scales = "free_y") +  # One boxplot per gene
  labs(title="Expression of immune checkpoints",
       y="Normalized CPM") +
  scale_y_continuous(expand=expansion(mult=c(0.1, 0.4))) 

#comparing adenoma polyp vs paired normal expression of immune checkpoin inhibitor genes
allls_icb_adenoma <- allls_icb |> filter(histology == "adenoma") |> arrange(gene, tmp_id)
icb_stats_ade <- tibble()
for (marker in unique(allls_icb_adenoma$gene)) {
  markers_data <- allls_icb_adenoma |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = markers_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  icb_stats_ade <- bind_rows(icb_stats_ade, wilcox_result)
}
icb_stats_ade$padj <- p.adjust(icb_stats_ade$p, method = "BH")
view(icb_stats_ade)

#comparing serrarated polyp vs paired normal expression of immune checkpoin inhibitor genes
allls_icb_serrated <- allls_icb |> filter(histology == "serrated") |> arrange(gene, tmp_id)
icb_stats_ser <- tibble()
for (marker in unique(allls_icb_serrated$gene)) {
  markers_data <- allls_icb_serrated |> filter(gene == marker)
  wilcox_result <- wilcox_test(CPM ~ condition, data = markers_data, paired = TRUE)
  wilcox_result <- wilcox_result |> mutate(gene = marker)
  icb_stats_ser <- bind_rows(icb_stats_ser, wilcox_result)
}
icb_stats_ser$padj <- p.adjust(icb_stats_ser$p, method = "BH")
view(icb_stats_ser)