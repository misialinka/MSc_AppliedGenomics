#loading pachages for analysis of maf files 
if (!require("BiocManager"))
  install.packages("BiocManager")
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)
library(maftools)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(VennDiagram)
library(ggVennDiagram)


#creating signle maf object from different sequencing batches -------
#loading a MAF file with 1st and 2nd batches of sequencing
maf_serrated1 <- read.maf("~/Desktop/serrated seq/DNAseq final/maf_ls_maftools_batch12.maf")

#loading a MAF file with 3rd batch of sequencing  
maf_serrated2 <- read.maf("~/Desktop/serrated seq/DNAseq final/data_mutations_batch4_5.maf")

#this file includes data from polyps and colorectal cancers. I am filtering out the cancer samples to keep only polyps. 
cancer_samples <- c("10766_tissue_normal_sigmoid", "10766_tissue_colon_sigmoid_vs_10766_tissue_normal_", "11520_tissue_normal_ascending",
                    "11520_tissue_colon_ascending_vs_11520_tissue_norma", "11521_tissue_normal_ascending", "11521_tissue_colon_ascending_vs_11521_tissue_norma", "11545_tissue_normal_descending", "11545_tissue_colon_descending_vs_11545_tissue_norm", "11547_tissue_normal_hepatic_flexure", "11547_tissue_colon_hepatic_flexure_vs_11547_tissue")

maf_serrated2 <- subsetMaf(maf_serrated2, tsb = setdiff(maf_serrated2@data$Tumor_Sample_Barcode, cancer_samples))

#merging the maf files 
maf_ls <- merge_mafs(list(maf_serrated1, maf_serrated2))

#loading metadata including clinical information about each samples
metadata_serrated <- read_csv("~/Desktop/serrated seq/DNAseq final/MAFmetadata_serrated.csv")

#extracting relevant information from the metadata clinical data 
clinical_serrated <- data.frame(
  Tumor_Sample_Barcode = metadata_serrated$`MAF name`,  #MAF sample names 
  condition = metadata_serrated$samples,  #condition (samples) specifies "polyp" or "normal" samples
  histology = metadata_serrated$histology, #histology defines if patient has serrated or adenoma polyps
  subject = metadata_serrated$`pair ID`,  #subject defines normal-polyp pairs from the same patient
  batch = metadata_serrated$`Sample batch`,  
  location = metadata_serrated$Location, #location in the gut from which the sample was taken
  side = metadata_serrated$`Location side`,  #specifies "left" or "right" side of the colon
  allhistologies = metadata_serrated$`histology all`, #includes categories "normal" "adenoma" "serrated"  
  method = metadata_serrated$Genomics,
  stringsAsFactors = FALSE
)

rownames(clinical_serrated) <- clinical_serrated$Tumor_Sample_Barcode

#checking if Tumour_Sample_Barcode in maf and clinical_serrated are the same 
all(rownames(clinical_serrated) %in% maf_ls@clinical.data$Tumor_Sample_Barcode) #TRUE

#matching the sequence of Tumour_Sample_Barcode in clinical data to match the maf file
clinical_serrated <- clinical_serrated[maf_ls@clinical.data$Tumor_Sample_Barcode ,]

#checking if the sequence of Tumour_Sample_Barcode is the same
all(rownames(clinical_serrated) == maf_ls@clinical.data$Tumor_Sample_Barcode) #TRUE

#adding clinical data to the maf files 
maf_ls <- read.maf(maf_ls@data, clinicalData = clinical_serrated)

#excluding samples with unknown histology, to focus only on adenoma and serrated samples
maf_ls <- subsetMaf(maf = maf_ls, clinQuery = "histology == 'adenoma' | histology == 'serrated'")

#excluding normal samples to keep only somatic mutations for further analysis 
maf_ls <- subsetMaf(maf = maf_ls, tsb = setdiff(maf_ls@data$Tumor_Sample_Barcode, c("10436N", "10842N", "11316N", "11316N", "11452N", "1624N")))

#creating subsetted maf files, seperate for each histology, for comparisons in later analysis 
maf_serrated <- subsetMaf(maf = maf_ls, clinQuery = "histology == 'serrated'")
maf_adenoma <- subsetMaf(maf = maf_ls, clinQuery = "histology == 'adenoma'")

#checking how many variants have been described in COSMIC
has_cosmic <- grepl("COS", maf_ls@data$Existing_variation) 
table(has_cosmic) #60,917 FALSE, 26,927 TRUE

#subsetting only variants that have been described in COSMIC, to keep more reliable variants 
maf_cosmic <- maf_ls
maf_cosmic@data <- maf_ls@data[grepl("COS", maf_ls@data$Existing_variation), ]

#keeping only HIGH and MODERATE impact variants
maf_cosmic <- subsetMaf(maf_cosmic, query = "IMPACT == 'HIGH' | IMPACT == 'MODERATE'")

#create subsets of histologies for comparisons in later analysis for cosmic mutations
maf_cosmic_serrated <- subsetMaf(maf = maf_cosmic, clinQuery = "histology == 'serrated'")
maf_cosmic_adenoma <-  subsetMaf(maf = maf_cosmic, clinQuery = "histology == 'adenoma'")

#summaries of mutations ------
#Writes maf summary to an output file with basename maf_ls
write.mafSummary(maf = maf_cosmic, basename = 'maf_ls')

#getting summary of all mutations in samples 
cosmic_summary <- getSampleSummary(maf_cosmic)
view(cosmic_summary)

#Tumour mutational burden measures  ---------
#plotting tumour mutational burden (tmb) for all samples. The graph shows that samples sequenced with WES form outliers in tmb analysis
tmb_WES <- tmb(maf_cosmic, logScale = FALSE)

#selecting samples which were sequenced with WGS (and not WES) for TMB calculations, since WES samples form outliers 
maf_ls_WGS <- subsetMaf(maf = maf_cosmic, clinQuery = "method == 'WGS'")

#subsetting mafs for each histology for tumour mutational burden calculations
maf_serrated_WGS <- subsetMaf(maf = maf_ls_WGS, clinQuery = "histology == 'serrated'")
maf_adenoma_WGS <- subsetMaf(maf = maf_ls_WGS, clinQuery = "histology == 'adenoma'")

#calculating tmb 
tmb_serrated <- tmb(maf = maf_serrated_WGS, captureSize = 50)
mean(tmb_serrated$total_perMB)
sd(tmb_serrated$total_perMB)
tmb_adenoma <- tmb(maf = maf_adenoma_WGS, captureSize = 50)
mean(tmb_adenoma$total_perMB)
sd(tmb_adenoma$total_perMB)
tmb_all <- tmb(maf = maf_ls_WGS, captureSize = 50)

# Combine TMB data into a single dataframe for plotting
tmb_comparison <- data.frame(
  Sample = c(tmb_serrated$Tumor_Sample_Barcode, tmb_adenoma$Tumor_Sample_Barcode),
  TMB = c(tmb_serrated$total_perMB, tmb_adenoma$total_perMB),
  Cohort = c(rep("Serrated", nrow(tmb_serrated)), rep("Adenoma", nrow(tmb_adenoma)))
)

#runnign statistical analysis (Wilcox) for tmb in adenoma and serrated polyps
wilcox.test(tmb_serrated$total_perMB, tmb_adenoma$total_perMB)

#Making a plot to compare TMB
ggplot(tmb_comparison, aes(x = Cohort, y = TMB, fill = Cohort)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +  # Add individual points
  theme_minimal() +
  stat_compare_means(comparison = list(c("Serrated", "Adenoma")), 
                     method = "wilcox.test", 
                     label = "p.format") +  # Adjust Y position
  labs(title = "Tumor Mutation Burden Comparison", y = "TMB (mut/Mb)")

#subsetting data to calculate frameshift mutation burden. Frameshift mutations are known to have more serious effects then the more common single nucleotide variants, and are also genomic signature of Lynch Syndrome.  
maf_WGS_fs <- subsetMaf(maf = maf_ls_WGS, query = "Variant_Classification == 'Frame_Shift_Del' | Variant_Classification == 'Frame_Shift_Ins'")
maf_serrated_fs <- subsetMaf(maf = maf_serrated_WGS, query = "Variant_Classification == 'Frame_Shift_Del' | Variant_Classification == 'Frame_Shift_Ins'")
maf_adenoma_fs <- subsetMaf(maf = maf_adenoma_WGS, query = "Variant_Classification == 'Frame_Shift_Del' | Variant_Classification == 'Frame_Shift_Ins'")

#calculating tmb for frameshifts
fsmb_serrated <- tmb(maf = maf_serrated_fs, captureSize = 50)
mean(fsmb_serrated$total_perMB)
sd(fsmb_serrated$total_perMB)
fsmb_adenoma <- tmb(maf = maf_adenoma_fs, captureSize = 50)
mean(fsmb_adenoma$total_perMB)
sd(fsmb_adenoma$total_perMB)
fsmb_all <- tmb(maf = maf_WGS_fs, captureSize = 50)

# Combine TMB data into a single dataframe
fsmb_comparison <- data.frame(
  Sample = c(fsmb_serrated$Tumor_Sample_Barcode, fsmb_adenoma$Tumor_Sample_Barcode),
  TMB = c(fsmb_serrated$total_perMB, fsmb_adenoma$total_perMB),
  Cohort = c(rep("Serrated", nrow(tmb_serrated)), rep("Adenoma", nrow(tmb_adenoma)))
)

#runnign statistical analysis (Wilcox) for tmb in adenoma and serrated polyps
wilcox.test(fsmb_serrated$total_perMB, fsmb_adenoma$total_perMB)

# Boxplot to compare TMB
ggplot(fsmb_comparison, aes(x = Cohort, y = TMB, fill = Cohort)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +  # Add individual points
  theme_minimal() +
  stat_compare_means(comparison = list(c("Serrated", "Adenoma")), 
                     method = "wilcox.test", 
                     label = "p.format") +  # Adjust Y position
  labs(title = "Tumor Frame Shift Mutation Burden Comparison", y = "FSMB (mut/Mb)")

#oncoplots for all samples ------
#oncoplot for top 50 mutated genes - I could see HLA-A and HLA-DRB1 and HLA-DRB5 among these genes. 
oncoplot(
  maf = maf_cosmic, top = 50,
  draw_titv = TRUE, 
  clinicalFeatures = c('histology', 'location'),
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = TRUE
)

#oncoplot including pathways. smgbp for biological pathways. sigpw for oncogenic pathways. check out ?oncoplot for more functions, susch as collapsing pathways. The two main mutated genes in the immune signalling pathway were HLA-A and HLA-B. 
dev.off()
oncoplot(maf = maf_cosmic, 
         draw_titv = TRUE,
         pathways = 'smgbp', topPathways = 5, 
         collapsePathway = FALSE,
         sortByAnnotation = TRUE,
         clinicalFeatures = c("histology", "location")
)

#plots for antigen presentation ------
#creating antigen presentation signature based on "The genomic landscape of 2,023 colorectal cancers", Cornish et al, 2024
hla_genes <- data_frame(
  Genes = c("HLA-A", "HLA-B", "HLA-C", "HLA_E", "HLA-F", "HLA-G", "B2M", "TAP1", "TAP2", "ERAP1", "ERAP2", "PDIA3", "CALR",  "PSME1", "PSME2", "PSME3", "PSMA7", "IRF1", "CANX", "NLRC5", "CIITA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"
  ), 
  Pathway = rep(c(
    "HLA class I", "Class I antigen processing", "HLA class II"
  ), 
  c(7, 13, 15)),
  stringsAsFactors = FALSE
)

#creating an oncoplot for my own signature for antigen presentation genes, for serrated and adenomatous
dev.off()
oncoplot(maf = maf_cosmic_adenoma, 
         draw_titv = TRUE,
         pathways = hla_genes, 
         gene_mar = 8.5,
         #collapsePathway = TRUE,
         sortByAnnotation = TRUE,
         clinicalFeatures = c("histology", "location"))
dev.off()
oncoplot(maf = maf_cosmic_serrated, 
         draw_titv = TRUE,
         pathways = hla_genes, 
         gene_mar = 8.5,
         #collapsePathway = TRUE,
         sortByAnnotation = TRUE,
         clinicalFeatures = c("histology", "location"))

#investigating if there are mutations in the immune checkpoin inhibitors 
icb_genes <- c("CD274", "CTLA4", "PDCD1", "HAVCR2", "TIGIT", "LAG3")

dev.off()
oncoplot(maf = maf_cosmic, 
         draw_titv = TRUE,
         genes = icb_genes,
         sortByAnnotation = TRUE,
         clinicalFeatures = c("histology", "location"))
