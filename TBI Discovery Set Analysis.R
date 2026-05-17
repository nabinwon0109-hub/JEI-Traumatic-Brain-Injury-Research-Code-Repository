#Discovery Set Analysis 
# 1. Install and Load Required Libraries  ###########################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
source("/Users/nabinwon/Desktop/RNA-seq/renv/activate.R") #가상환경 설정
library(DESeq2)



# 2. Read the Count Matrix  ###########################

file_path <- "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/GSE173975_TBI_GENE_MATRIX_revised.txt"

# Read the count matrix
count_data <- read.table(file_path, sep = "\t", header = TRUE)

rownames(count_data) <- count_data$Gene

head(count_data)
count_data$Gene <- NULL


colnames(count_data)

ordered_cols <- c(
  "Sham.1.rep1", "Sham.1.rep2",       # Sham Day1
  "TBI.1.rep1", "TBI.1.rep2", "TBI.1.rep3", "TBI.1.rep4",  # TBI Day1
  "Sham.14.rep1", "Sham.14.rep2", "Sham.14.rep3", "Sham.14.rep4",  # Sham Day14
  "TBI.14.rep1", "TBI.14.rep2", "TBI.14.rep3", "TBI.14.rep4"  # TBI Day14
)

count_data_ordered <- count_data[, c(ordered_cols)]
colnames(count_data_ordered) <- gsub("\\.", "", ordered_cols)
colnames(count_data_ordered)
# 3. Prepare Metadata (Coldata) ###########################

condition <- c("Sham", "Sham",  # Sham Day1
               "TBI", "TBI", "TBI", "TBI",  # TBI Day1
               "Sham", "Sham", "Sham", "Sham",  # Sham Day14
               "TBI", "TBI", "TBI", "TBI")  # TBI Day14

# Define the sample days (Day1 vs. Day14) based on the reordered columns
day <- c("Day1", "Day1",  # Sham Day1
         "Day1", "Day1", "Day1", "Day1",  # TBI Day1
         "Day14", "Day14", "Day14", "Day14",  # Sham Day14
         "Day14", "Day14", "Day14", "Day14")  # TBI Day14

# Combine the two factors into a combined variable
conditionday <- paste(condition, day, sep = "")

# Combine into a metadata dataframe
coldata <- data.frame(condition = factor(condition), 
                      day = factor(day), 
                      conditionday = factor(conditionday))

# Set the row names to match the column names of the reordered count data
rownames(coldata) <- colnames(count_data_ordered)

# Check the colData structure
coldata



# 4. Create DESeq2 Dataset ###########################


# Proceed with DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = count_data_ordered, colData = coldata, design = ~ conditionday)
dds

smallestGroupSize <- ncol(count_data_ordered) * 0.5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds





# 5. PCA & Transformation ###########################
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)



plotPCA(vsd, intgroup=c("condition", "day"))




# 6. Run DESeq2 Analysis ###########################

# dds$conditionday <- factor(dds$conditionday, levels = c("ShamDay1",  "TBIDay1","ShamDay14", "TBIDay14"))
dds$conditionday <- relevel(dds$conditionday, ref = "ShamDay14")

# By setting ref = "ShamDay1", ShamDay1 is used as the reference and all other groups will be compared against it.


# Proceed with DESeq2 analysis
dds <- DESeq(dds)

# Get the results
res_day14_shamVSTBI <- results(dds, contrast = c("conditionday", "TBIDay14", "ShamDay14"))
head(res_day14_shamVSTBI)




## Log fold change shrinkage for visualization and ranking #######################
resultsNames(dds)

# The name used for the coef argument in the lfcShrink() function must match the exact name provided by resultsNames(dds).

res_sham14vs_TBI14 <- lfcShrink(dds, 
                                coef = "conditionday_TBIDay14_vs_ShamDay14", 
                                type = "apeglm")

res_sham14vs_Sham1 <- lfcShrink(dds, 
                                coef = "conditionday_ShamDay14_vs_ShamDay1", 
                                type = "apeglm")

res_TBI1vs_TBI14 <- lfcShrink(dds, 
                              coef = "conditionday_TBIDay1_vs_ShamDay1", 
                              type = "apeglm")



head(res_sham14vs_TBI14)

UpDEG <- subset(res_sham14vs_TBI14, log2FoldChange > 0.5 & padj < 0.05)
DnDEG <- subset(res_sham14vs_TBI14, log2FoldChange < -0.5 & padj < 0.05)


UpDEG 
DnDEG 
up_genes



# 7. Heatmap ##########################

library(pheatmap)

up_genes <- rownames(UpDEG)
dn_genes <- rownames(DnDEG)

# Combine both sets of genes for the heatmap
all_genes <- c(up_genes, dn_genes)


normalized_counts <- counts(dds, normalized = TRUE)

# Subset the normalized counts for the DEGs
deg_counts <- normalized_counts[all_genes, ]


pheatmap(deg_counts, 
         scale = "row",  # Scale the values across each row (gene) to z-scores
         cluster_rows = TRUE,  # Cluster genes
         cluster_cols = TRUE,  # Cluster samples
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE)  # Show sample names


normalized_counts[up_genes[1]  ,]


# 8. VolcanoPlot ##########################

library(ggplot2)


volcano_data <- as.data.frame(res_sham14vs_TBI14)
volcano_data$significant <- volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) >= 1



ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(aes(color=ifelse(significant, 
                              ifelse(log2FoldChange > 0, "red", "blue"), 
                              "grey")), 
             alpha=0.7, size=1.5) +  # 색상 및 점의 투명도와 크기 조정
  scale_color_identity() +  # 색상 매핑
  geom_vline(xintercept=c(-2, 2), linetype="dashed", color="black") +  # 수직선 추가
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black") +  # 수평선 추가
  theme_classic() +  # 기본 테마
  labs(x="Log2 Fold Change", 
       y="-Log10(adjusted p-value)", 
       title="Volcano Plot of Differential Expression") +  # 레이블 및 제목 설정
  theme(plot.title = element_text(hjust = 0.5))  # 제목 중앙 정렬




# 9. Enrichment analysis ##########################


# 1. Install and load the necessary packages
if (!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

library(enrichR)
library(writexl)

up_genes <- rownames(UpDEG)
dn_genes <- rownames(DnDEG)


enrichr_libraries <- listEnrichrDbs()

dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2021")




# Run enrichment for UpDEG
up_enrichment <- enrichr(up_genes, dbs)

# Run enrichment for DnDEG
dn_enrichment <- enrichr(dn_genes, dbs)

# 5. Convert enrichment results to a list of data frames
up_enrichment_list <- lapply(up_enrichment, as.data.frame)
dn_enrichment_list <- lapply(dn_enrichment, as.data.frame)

# 6. Save the enrichment results to an Excel file
write_xlsx(list("UpDEG_GO_BP" = up_enrichment_list$`GO_Biological_Process_2021`,
                "UpDEG_KEGG" = up_enrichment_list$`KEGG_2021_Human`,
                "UpDEG_Reactome" = up_enrichment_list$`Reactome_2021`,
                "DnDEG_GO_BP" = dn_enrichment_list$`GO_Biological_Process_2021`,
                "DnDEG_KEGG" = dn_enrichment_list$`KEGG_2021_Human`,
                "DnDEG_Reactome" = dn_enrichment_list$`Reactome_2021`),
           path = "Enrichr_Results.xlsx")


write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "BioCarta_up" = up_enrichment_list$BioCarta_2016,
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/TBI_sham14vsTBI14discoveryset.xlsx")



## Fgsea 

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(fgsea)
library(gridExtra)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")


# IL-6/JAK/STAT3 Signaling
# KRAS Signaling Up
# p53 Pathway
# PI3K/AKT/mTOR  Signaling
# TGF-beta Signaling
# 
# CDK Regulation of DNA Replication Homo sapiens h mcmPathway
# IL-2 Receptor Beta Chain in T cell Activation Homo sapiens h il2rbPathway
# Stress Induction of HSP Regulation Homo sapiens h hsp27Pathway
# TNFR1 Signaling Pathway Homo sapiens h tnfr1Pathway
# NF-kB Signaling Pathway Homo sapiens h nfkbPathway




gene_set_genelist <- list(
  "IL6_JAK_STAT3_Signaling" = c("PTPN1", "TGFB1", "CCR1", "PTPN2", "IL4R", "JUN", "IL1R1", "IL1R2", "ITGA4", "TNFRSF12A", "IL10RB", "PLA2G2A", "DNTT", 
                                "PTPN11", "TYK2", "OSMR", "TNFRSF1A", "TNFRSF1B", "CXCL10", "HAX1", "TLR2", "CXCL11", "IL1B", "IRF1", "TNFRSF21", 
                                "IL7", "CXCL13", "IL6", "CD9", "FAS", "IL3RA", "IL9R", "IRF9", "IL6ST", "PF4", "STAM2", "CRLF2", "IL18R1", "CNTFR", 
                                "ACVRL1", "CSF3R", "CXCL9", "CSF2", "CSF1", "EBI3", "PIK3R5", "TNF", "CXCL1", "IL2RG", "CXCL3", "SOCS3", "SOCS1", 
                                "CD38", "CD36", "IL12RB1", "MAP3K8", "IL17RA", "IL15RA", "IFNGR1", "STAT1", "IFNGR2", "STAT2", "STAT3", "INHBE", 
                                "IL17RB", "MYD88", "LTBR", "IL2RA", "LTB", "GRB2", "IFNAR1", "CD44", "REG1A", "ITGB3", "PDGFC", "CSF2RA", "ACVR1B", 
                                "CBL", "CSF2RB", "CCL7", "HMOX1", "LEPR", "PIM1", "BAK1", "A2M", "CD14", "IL13RA1"),
  
  "KRAS_Signaling_Up" = c("ATG10", "CIDEA", "EREG", "MALL", "NIN", "RBP4", "SLPI", "ITGBL1", "SPARCL1", "ZNF277", "ADAMDEC1", "SATB1", "C3AR1", 
                          "RETN", "GLRX", "IL2RG", "TFPI", "ADGRA2", "AKT2", "RABGAP1L", "ST6GAL1", "GADD45G", "PPBP", "AMMECR1", "NR0B2", "IGF2", 
                          "BTC", "TNNT2", "CCSER2", "STRN", "CFB", "SPON1", "CFH", "PLEK2", "IKZF1", "GNG11", "EPB41L3", "TSPAN7", "ANKH", "SOX9", 
                          "KCNN4", "CMKLR1", "TSPAN1", "APOD", "SCN1B", "GABRA3", "DUSP6", "FUCA1", "ACE", "JUP", "CCL20", "IL10RA", "EMP1", 
                          "PDCD1LG2", "SNAP91", "DCBLD2", "MMP11", "WDR33", "MMP10", "VWA5A", "IL1B", "LAT2", "GPRC5B", "CFHR2", "ID2", "TLR8", 
                          "SDCCAG8", "ADGRL4", "AVL9", "SERPINA3", "TPH1", "GYPC", "PTPRR", "CAB39L", "CROT", "TOR1AIP2", "CTSS", "MAP7", "CD37", 
                          "FCER1G", "MPZL2", "PRELID3B", "HSD11B1", "IL33", "MMD", "MAP4K1", "ARG1", "TMEM176A", "G0S2", "TMEM176B", "WNT7A", 
                          "LAPTM5", "NR1H4", "ANO1", "TMEM158", "ADAM17", "ANGPTL4", "F2RL1", "ENG", "HBEGF", "NRP1", "ERO1A", "MTMR10", "F13A1", 
                          "PIGR", "HOXD11", "PTGS2", "PSMB8", "AKAP12", "CXCR4", "GPNMB", "CBR4", "MAP3K1", "PTCD2", "GUCY1A1", "PRRX1", "H2BC3", 
                          "CXCL10", "MYCN", "PCP4", "CLEC4A", "ADAM8", "SCG5", "SCG3", "PLVAP", "CSF2", "PPP1R15A", "RGS16", "PRDM1", "RELN", 
                          "CCND2", "KIF5C", "SPP1", "IGFBP3", "INHBA", "TRAF1", "MAFB", "DNMBP", "BIRC3", "YRDC", "IL7R", "PLAU", "PLAT", 
                          "CSF2RA", "ETS1", "TMEM100", "FGF9", "BTBD3", "PCSK1N", "CDADC1", "CBX8", "GALNT3", "PEG3", "GFPT2", "ITGA2", "ETV1", 
                          "NAP1L2", "ETV4", "MMP9", "ETV5", "TNFRSF1B", "BMP2", "LCP1", "CPE", "ZNF639", "IRF8", "PECAM1", "DOCK2", "SNAP25", 
                          "ABCB1", "TNFAIP3", "LY96", "PTBP2", "CA2", "EVI5", "PRKG2", "ANXA10", "PLAUR", "BPGM", "KLF4", "NGF", "ALDH1A3", 
                          "ALDH1A2", "HKDC1", "SPRY2", "TRIB2", "TRIB1", "ITGB2", "FLT4", "SEMA3B", "USP12", "HDAC9", "TSPAN13", "CBL", "USH1C", 
                          "RBM4", "IL1RL2", "FBXO4", "EPHB2"),
  
  "p53_Pathway" = c("HEXIM1", "PRKAB1", "TM4SF1", "NINJ1", "CCP110", "BAX", "TNFSF9", "EPHA2", "TP53", "APP", "ATF3", "CYFIP2", "RB1", "GLS2", 
                    "CCNK", "EI24", "RETSAT", "ABAT", "CDKN2AIP", "TOB1", "NHLH2", "STOM", "S100A10", "RPL18", "JAG2", "TSPYL2", "PMM1", "ST14", 
                    "GADD45A", "PLK3", "PLK2", "DNTTIP2", "H1-2", "ALOX15B", "VAMP8", "TCN2", "SDC1", "HSPA4L", "POM121", "RRP8", "SLC35D1", 
                    "PDGFA", "XPC", "IFI30", "FDXR", "HMOX1", "RPL36", "SFN", "IER5", "RPS12", "IER3", "SLC19A2", "TP63", "ZNF365", "FUCA1", 
                    "BLCAP", "JUN", "TSC22D1", "SPHK1", "LIF", "VWA5A", "IL1A", "PTPRE", "RAD51C", "RAP2B", "H2AJ", "KRT17", "CCNG1", "COQ8A", 
                    "TAX1BP3", "ADA", "H2AW", "CSRNP2", "BTG2", "BTG1", "RRAD", "ABHD4", "RPS27L", "TGFA", "SLC7A11", "KLK8", "HINT1", "LDHB", 
                    "SOCS1", "CTSF", "CTSD", "SERPINB5", "CDKN2A", "CDKN2B", "PRMT2", "ELP1", "FOS", "WWP1", "AEN", "DDB2", "CDK5R1", "HBEGF", 
                    "ANKRA2", "CDKN1A", "SERTAD3", "KIF13B", "HRAS", "RACK1", "NUDT15", "TGFB1", "FBXW7", "ZBTB16", "IRAG2", "TAP1", "CGRRF1", 
                    "DEF6", "NUPR1", "FAS", "TRAFD1", "MAPKAPK3", "ISCU", "PLXNB2", "S100A4", "STEAP3", "PCNA", "CD82", "CD81", "DGKA", "PPP1R15A", 
                    "PVT1", "RGS16", "NOL8", "ZFP36L1", "FOXO3", "RCHY1", "NDRG1", "RXRA", "CCND3", "RAB40C", "CCND2", "RNF19B", "PERP", 
                    "PHLDA3", "EPHX1", "F2R", "INHBB", "MXD1", "TRAF4", "PROCR", "CDH13", "FGF13", "MXD4", "RAD9A", "MDM2", "NOTCH1", "AK1", 
                    "TRIAP1", "SAT1", "TM7SF3", "SESN1", "ZMAT3", "FAM162A", "CLCA2", "ABCC5", "VDR", "GPX2", "APAF1", "OSGIN1", "TPRKB", 
                    "PTPN14", "BAIAP2", "SP1", "BMP2", "SLC3A2", "PITPNC1", "DRAM1", "MKNK2", "IRAK1", "CASP1", "IP6K2", "RALGDS", "PIDD1", 
                    "RHBDF2", "DCXR", "KLF4", "TPD52L1", "WRAP73", "PPM1D", "ERCC5", "DDIT4", "TXNIP", "DDIT3", "TRIB3", "CEBPA", "HDAC3", 
                    "ITGB4", "TCHH", "ACVR1B", "GM2A", "SEC61A1", "EPS8L2", "BAK1", "UPP1", "POLH", "TNNI1"),
  
  "PI3K_AKT_mTOR_Signaling" = c("TSC2", "PPP1CA", "MAPK10", "SQSTM1", "PFN1", "RAF1", "UBE2N", "PRKAA2", "RALB", "ATF1", "SLA", "UBE2D3", "NOD1", 
                                "ADCY2", "IL2RG", "EGFR", "HSP90B1", "E2F1", "AKT1", "CFL1", "MAP3K7", "PAK4", "PITX2", "FGF22", "SMAD2", 
                                "MAP2K3", "PRKCB", "TRAF2", "FGF17", "MYD88", "ARPC3", "CSNK2B", "GRB2", "PIN1", "SLC2A1", "FASLG", "DAPP1", 
                                "FGF6", "TBK1", "THEM4", "PRKAR2A", "PLCG1", "SFN", "PDK1", "MAP2K6", "AP2M1", "ACTR3", "ACTR2", "DUSP3", 
                                "PTPN11", "IRAK4", "TNFRSF1A", "VAV3", "NFKBIB", "IL4", "ECSIT", "PIKFYVE", "AKT1S1", "CDK4", "CDK2", 
                                "CAMK4", "CDK1", "CALR", "GSK3B", "YWHAB", "ARF1", "CAB39L", "MAPKAP1", "CLTC", "ITPR2", "PIK3R3", "ACACA", 
                                "GRK2", "GNGT1", "MKNK1", "MKNK2", "RIPK1", "EIF4E", "PLA2G12A", "STAT2", "NGF", "TIAM1", "LCK", "RIT1", 
                                "DDIT3", "TRIB3", "PLCB1", "CDKN1B", "CDKN1A", "CAB39", "PTEN", "PRKAG1", "GNA14", "ARHGDIA", "RPTOR", 
                                "MAPK9", "MAPK8", "RPS6KA1", "PPP2R1B", "CXCR4", "RPS6KA3", "HRAS", "MAPK1", "RAC1", "NCK1"),
  
  "TGF_beta_Signaling" = c("TGFB1", "ACVR1", "RHOA", "SMURF2", "SMURF1", "HIPK2", "TGFBR1", "PPP1CA", "CDK9", "SKI", "RAB31", "APC", 
                           "BMP2", "NCOR2", "ID1", "ID3", "ID2", "BMPR1A", "BCAR3", "BMPR2", "SERPINE1", "UBE2D3", "PPP1R15A", "ARID4B", 
                           "CDH1", "MAP3K7", "SKIL", "SPTBN1", "SMAD3", "WWTR1", "TGIF1", "LEFTY2", "NOG", "IFNGR2", "SMAD1", "KLF10", 
                           "SMAD7", "SMAD6", "PPM1A", "FKBP1A", "TJP1", "XIAP", "ENG", "TRIM33", "CDKN1C", "HDAC1", "FURIN", "CTNNB1", 
                           "SLC20A1", "THBS1", "LTBP2", "FNTA", "PMEPA1", "JUNB"),
  
  "CDK_Regulation_of_DNA_Replication" = c("CDKN1B", "CDT1", "CDC6", "ORC6", "ORC5", "KITLG", "ORC2", "ORC1", "CCNE1", "ORC4", "ORC3", 
                                          "CDK2", "MCM4", "MCM5", "MCM6", "MCM7", "MCM2", "MCM3"),
  
  "IL_2_Receptor_Beta_Chain_in_T_cell_Activation" = c("BAD", "NMI", "IL2", "PPIA", "PIK3CA", "FAS", "BCL2", "PTPN6", "SOS1", "RAF1", 
                                                      "PCNA", "RB1", "SHC1", "CCNH", "IL2RG", "SOCS3", "CRKL", "SOCS1", "CCND3", "CCND2", 
                                                      "CCND1", "E2F1", "MYC", "AKT1", "STAT5B", "MAP2K2", "MAP2K1", "PDPK1", "FOS", 
                                                      "TFDP1", "CCNA1", "PIK3R1", "IL2RA", "RPS6KB1", "CCNE1", "IL2RB", "BCL2L1", "GRB2", 
                                                      "IRS1", "FASLG", "IKZF3", "CBL", "CCNB1", "MAPK3", "HRAS", "MAPK1", "JAK3", "JAK1", 
                                                      "STAT5A"),
  
  "Stress_Induction_of_HSP_Regulation" = c("DAXX", "APAF1", "TNF", "FASLG", "HSPB1", "IL1A", "ACTA1", "CASP9", "FAS", "MAPKAPK2", 
                                           "MAPKAPK3", "BCL2", "CASP3", "CYCS"),
  
  "TNFR1_Signaling_Pathway" = c("MAP2K4", "TRADD", "CRADD", "LMNB1", "TNF", "TRAF2", "TNFRSF1A", "LMNB2", "MADD", "LMNA", "CASP8", 
                                "BAG4", "CASP2", "CASP3", "FADD", "RIPK1", "BIRC3"),
  
  "NF_kB_Signaling_Pathway" = c("RELA", "MAP3K1", "IL1R1", "CHUK", "TRADD", "TNF", "NFKB1", "TNFRSF1A", "IL1A", "IKBKB", "NFKBIA", 
                                "MYD88", "TRAF6", "IRAK1", "FADD", "RIPK1", "IKBKG", "TAB1", "MAP3K7", "MAP4K4", "TLR4"), 
  
  "TGF_beta_Signaling_Pathway" = c("TGFB1", "ACVR1", "RHOA", "TGFB3", "TGFB2", "SMURF2", "SMURF1", "BMP8B", "PPP2R1A", "HFE2", "BMP8A", 
                                   "TGFBR1", "ACVR2B", "RBX1", "RGMA", "BMP7", "TGFBR2", "BMP6", "ACVR2A", "RGMB", "BMP5", "BMP4", 
                                   "RBL1", "SP1", "ZFYVE16", "BMP2", "ID1", "BAMBI", "ID3", "IFNG", "ID2", "ID4", "NBL1", "BMPR1A", 
                                   "BMPR1B", "BMPR2", "ROCK1", "ZFYVE9", "TNF", "MYC", "E2F4", "AMH", "E2F5", "PITX2", "SMAD3", 
                                   "TGIF2", "SMAD2", "4930544G11RIK", "CREBBP", "TGIF1", "SMAD5", "LEFTY2", "SMAD4", "LEFTY1", 
                                   "NOG", "CDKN2B", "SMAD1", "FST", "INHBA", "GDF5", "INHBC", "INHBB", "GDF7", "INHBE", "SMAD7", 
                                   "GDF6", "SMAD6", "SMAD9", "PPP2CA", "DCN", "TFDP1", "RPS6KB2", "RPS6KB1", "HAMP", "FMOD", 
                                   "NODAL", "HAMP2", "AMHR2", "CHRD", "CUL1", "THBS1", "ACVR1C", "LTBP1", "ACVR1B", "PPP2CB", 
                                   "SKP1A", "PPP2R1B", "EP300", "MAPK3", "MAPK1", "NEO1")
  
  
  
)

neutrophil_genelist <- list(
  "neutrophil_degranulation" = c(
    "RETN", "ASAH1", "CRACR2A", "PSMD7", "PSMD6", "PSMD3", "ALOX5", "PSMD1", "PSMD2", 
    "STOM", "NCSTN", "PRSS3", "S100A12", "ELANE", "S100A11", "ATP6AP2", "PRKCD", "NCKAP1L", 
    "HRNR", "VAMP7", "VAMP8", "BIN2", "TCN1", "VAMP2", "PADI2", "CFD", "GMFG", "A1BG", 
    "GNS", "ABCA13", "CFP", "C3", "CYB5R3", "PRDX4", "PRDX6", "PLAC8", "NAPRT", "HMOX2", 
    "LTA4H", "ATG7", "FUCA1", "GSDMD", "ANXA2", "ANXA3", "JUP", "FUCA2", "NFAM1", "AZU1", 
    "TMC6", "LYZ", "SIGLEC5", "SIGLEC9", "GPI", "PAFAH1B2", "SERPINA3", "SERPINA1", "FPR2", 
    "FPR1", "IQGAP1", "LRRC7", "IQGAP2", "TIMP2", "AP1M1", "DYNLL1", "SERPINB3", "AHSG", 
    "ARG1", "DGAT1", "PGAM1", "SERPINB1", "ADAM10", "STING1", "LRG1", "RNASE3", "RNASE2", 
    "TUBB4B", "SERPINB6", "ERP44", "FABP5", "PSMA5", "PSMA2", "SLCO4C1", "LAMTOR1", "LAMTOR3", 
    "PLEKHO2", "LAMTOR2", "SERPINB12", "GRN", "PIGR", "PTAFR", "PYGB", "GSTP1", "HEXB", 
    "SERPINB10", "PYGL", "LPCAT1", "PSMB7", "PSMB1", "METTL7A", "AGA", "PGLYRP1", "OLR1", 
    "GSN", "IST1", "AGL", "AOC1", "DEFA4", "DEFA1", "ATP11B", "ATP11A", "ILF2", "TICAM2", 
    "CLEC4C", "CLEC4D", "RAB31", "OLFM4", "RAB37", "S100A7", "ARL8A", "ADAM8", "DSC1", 
    "PTPN6", "GHDC", "PKP1", "S100A9", "FRK", "S100A8", "CD177", "SLC27A2", "ROCK1", 
    "HP", "EPX", "APEH", "NRAS", "PYCARD", "LAMP2", "LAMP1", "RAB24", "KCMF1", "RAB4B", 
    "TRAPPC1", "ATP8B4", "ARPC5", "KCNAB2", "EEF2", "C6ORF120", "NFKB1", "ACLY", "DDOST", 
    "HPSE", "LCN2", "FTL", "NIT2", "RAB5B", "RAB3D", "RAB3A", "ARHGAP9", "ATP8A1", "BRI3", 
    "HVCN1", "GPR84", "RAB44", "STK10", "GLIPR1", "MMP25", "CLEC5A", "TXNDC5", "TOM1", 
    "GYG1", "TMEM179B", "SLC15A4", "ACTR2", "HSP90AA1", "APAF1", "PA2G4", "NDUFC2", 
    "TNFRSF1B", "TMBIM1", "DNAJC3", "FCGR2A", "NPC2", "DNAJC5", "DSG1", "PECAM1", 
    "VPS35L", "SNAP25", "SNAP23", "UBR4", "FCGR3B", "ALAD", "PSAP", "ACP3", "CEP290", 
    "SNAP29", "CR1", "DSP", "RNASET2", "CPPED1", "PKM", "CHIT1", "KPNB1", "NFASC", 
    "PFKL", "RAB10", "TSPAN14", "VAPA", "IMPDH1", "IMPDH2", "P2RX1", "RAB14", "MNDA", 
    "BPI", "RAB18", "FOLR3", "ITGB2", "CANT1", "HSP90AB1", "TCIRG1", "ITGAM", "FCAR", 
    "ITGAL", "TTR", "GM2A", "ITGAX", "ITGAV", "PGM2", "PGM1", "KRT1", "PTPRN2", 
    "NBEAL2", "MAPK14", "CRISPLD2", "SLPI", "TMEM63A", "FGL2", "STXBP2", "STXBP3", 
    "CALML5", "PRCP", "PNP", "CPNE1", "CPNE3", "B2M", "ARSA", "ARSB", "STK11IP", 
    "SURF4", "XRCC5", "SLC11A1", "XRCC6", "COMMD3", "HLA-C", "MIF", "PPBP", "PRG3", 
    "RAB27A", "PRG2", "HLA-B", "VAT1", "VNN1", "CSNK2B", "CYSTM1", "MGST1", "TRPM2", 
    "COMMD9", "LGALS3", "ANPEP", "FTH1", "CD14", "GOLGA7", "ATP6V0A1", "CTSB", 
    "CTSA", "RHOA", "PTGES2", "SIGLEC14", "TUBB", "AMPD3", "LILRB2", "RHOF", "LILRB3", 
    "RHOG", "NME2", "DERA", "TLR2", "BST2", "CEACAM3", "BST1", "PTPRC", "RAP2B", 
    "CEACAM1", "RAP2C", "CEACAM8", "PTPRB", "CEACAM6", "ALDOC", "ACTR10", "ALDOA", 
    "RAB7A", "RAB5C", "SLC44A2", "HBB", "CTSZ", "PTPRJ", "CTSS", "DSN1", "ACTR1B", 
    "RAP1A", "RAP1B", "MAGT1", "NEU1", "ORMDL3", "MLEC", "QSOX1", "CTSH", "CTSG", 
    "CD36", "FCER1G", "CTSD", "CAMP", "CD33", "CTSC", "RAB6A", "CAP1", "CD53", 
    "CCT2", "DBNL", "MME", "GAA", "CREG1", "QPCT", "ANO6", "PGRMC1", "APRT", 
    "SVIP", "RAB9B", "CD47", "DNASE1L1", "DYNC1H1", "CD44", "CD63", "CSTB", 
    "ORM2", "ORM1", "GDI2", "FAF2", "CAB39", "AP2A2", "LILRA3", "CST3", "CXCR2", 
    "CXCR1", "STBD1", "RAC1", "CD59", "ENPP4", "CCT8", "CD58", "SPTAN1", "CD55", 
    "LAIR1", "TBC1D10C", "CHRNB4", "CD300A", "GCA", "IRAG2", "PDAP1", "HUWE1", 
    "TYROBP", "ATAD3B", "OSCAR", "EEF1A1", "DOK3", "PPIA", "CD68", "DEGS1", 
    "YPEL5", "PPIE", "FCN1", "MPO", "ARHGAP45", "COTL1", "CAPN1", "OSTF1", 
    "IDH1", "IGF2R", "MGAM", "ATP6V0C", "MANBA", "CD93", "DNAJC13", "CYBB", 
    "CYBA", "SDCBP", "CKAP4", "DIAPH1", "TARM1", "SELL", "PRTN3", "VCL", 
    "CRISP3", "ATP6V1D", "VCP", "C5AR1", "COPB1", "SLC2A3", "PLAU", "PSEN1", 
    "SLC2A5", "SRP14", "CNN2", "ADGRG3", "ALDH3B1", "SIRPA", "ACAA1", "HGSNAT", 
    "MMP8", "HSPA6", "CMTM6", "GGH", "HSPA8", "MMP9", "ADA2", "MAN2B1", 
    "MS4A3", "DOCK2", "CHI3L1", "PSMD11", "B4GALT1", "DDX3X", "PSMD13", 
    "PSMD12", "PSMD14", "TNFAIP6", "MOSPD2", "HMGB1", "CXCL1", "AGPAT2", 
    "PLD1", "HK3", "SIRPB1", "HEBP2", "ADGRE3", "SYNGR1", "DPP7", "ADGRE5", 
    "TMEM30A", "PLAUR", "FGR", "GLB1", "CAT", "S100P", "DYNLT1", "PTX3", 
    "LTF", "CLEC12A", "MVP", "CAND1", "MAPK1", "GUSB"),
  
  "neutrophil_activation_involved_in_immuneresponse"= c(
    "UNC13D", "ARMC8", "MCEMP1", "SCAMP1", "DYNC1LI1", "FLG2", "GALNS", "TOLLIP", "PSMC2", 
    "PSMC3", "HSPA1A", "GLA", "CDK13", "NHLRC3", "CDA", "HSPA1B", "CYFIP1", "PRSS2", 
    "PDXK", "C3AR1", "FRMPD3", "RETN", "ASAH1", "CRACR2A", "PSMD7", "PSMD6", "PSMD3", 
    "ALOX5", "PSMD1", "PSMD2", "STOM", "NCSTN", "PRSS3", "S100A12", "ELANE", "S100A11", 
    "ATP6AP2", "PRKCD", "NCKAP1L", "HRNR", "VAMP7", "VAMP8", "BIN2", "TCN1", "VAMP2", 
    "PADI2", "CFD", "GMFG", "A1BG", "GNS", "ABCA13", "CFP", "C3", "CYB5R3", "PRDX4", 
    "PRDX6", "PLAC8", "NAPRT", "HMOX2", "LTA4H", "ATG7", "FUCA1", "GSDMD", "ANXA2", 
    "ANXA3", "JUP", "FUCA2", "NFAM1", "AZU1", "TMC6", "LYZ", "SIGLEC5", "SIGLEC9", 
    "GPI", "PAFAH1B2", "SERPINA3", "SERPINA1", "FPR2", "FPR1", "IQGAP1", "LRRC7", 
    "IQGAP2", "TIMP2", "AP1M1", "DYNLL1", "SERPINB3", "AHSG", "ARG1", "DGAT1", "PGAM1", 
    "SERPINB1", "ADAM10", "STING1", "LRG1", "RNASE3", "RNASE2", "TUBB4B", "SERPINB6", 
    "ERP44", "FABP5", "PSMA5", "PSMA2", "SLCO4C1", "LAMTOR1", "LAMTOR3", "PLEKHO2", 
    "LAMTOR2", "SERPINB12", "GRN", "PIGR", "PTAFR", "PYGB", "GSTP1", "HEXB", "SERPINB10", 
    "PYGL", "LPCAT1", "PSMB7", "PSMB1", "METTL7A", "AGA", "PGLYRP1", "OLR1", "GSN", 
    "IST1", "AGL", "AOC1", "DEFA4", "DEFA1", "ATP11B", "ATP11A", "ILF2", "TICAM2", 
    "CLEC4C", "CLEC4D", "RAB31", "OLFM4", "RAB37", "S100A7", "ARL8A", "ADAM8", "DSC1", 
    "PTPN6", "GHDC", "PKP1", "S100A9", "FRK", "S100A8", "CD177", "SLC27A2", "ROCK1", 
    "HP", "EPX", "APEH", "NRAS", "PYCARD", "LAMP2", "LAMP1", "RAB24", "KCMF1", "RAB4B", 
    "TRAPPC1", "ATP8B4", "ARPC5", "KCNAB2", "EEF2", "C6ORF120", "NFKB1", "ACLY", "DDOST", 
    "HPSE", "LCN2", "FTL", "NIT2", "RAB5B", "RAB3D", "RAB3A", "ARHGAP9", "ATP8A1", "BRI3", 
    "HVCN1", "GPR84", "RAB44", "STK10", "GLIPR1", "MMP25", "CLEC5A", "TXNDC5", "TOM1", 
    "GYG1", "TMEM179B", "SLC15A4", "ACTR2", "HSP90AA1", "APAF1", "PA2G4", "NDUFC2", 
    "TNFRSF1B", "TMBIM1", "DNAJC3", "FCGR2A", "NPC2", "DNAJC5", "DSG1", "PECAM1", "VPS35L", 
    "SNAP25", "SNAP23", "UBR4", "FCGR3B", "ALAD", "PSAP", "ACP3", "CEP290", "SNAP29", 
    "CR1", "DSP", "RNASET2", "CPPED1", "PKM", "CHIT1", "KPNB1", "NFASC", "PFKL", "RAB10", 
    "TSPAN14", "VAPA", "IMPDH1", "IMPDH2", "P2RX1", "RAB14", "MNDA", "BPI", "RAB18", 
    "FOLR3", "ITGB2", "CANT1", "HSP90AB1", "TCIRG1", "ITGAM", "FCAR", "ITGAL", "TTR", 
    "GM2A", "ITGAX", "ITGAV", "PGM2", "PGM1", "KRT1", "PTPRN2", "NBEAL2", "MAPK14", 
    "CRISPLD2", "SLPI", "TMEM63A", "FGL2", "STXBP2", "STXBP3", "CALML5", "PRCP", "PNP", 
    "CPNE1", "CPNE3", "B2M", "ARSA", "ARSB", "STK11IP", "SURF4", "XRCC5", "SLC11A1", 
    "XRCC6", "COMMD3", "HLA-C", "MIF", "PPBP", "PRG3", "RAB27A", "PRG2", "HLA-B", 
    "VAT1", "VNN1", "CSNK2B", "CYSTM1", "MGST1", "TRPM2", "COMMD9", "LGALS3", "ANPEP", 
    "FTH1", "CD14", "GOLGA7", "ATP6V0A1", "CTSB", "CTSA", "RHOA", "PTGES2", "SIGLEC14", 
    "TUBB", "AMPD3", "LILRB2", "RHOF", "LILRB3", "RHOG", "NME2", "DERA", "TLR2", "BST2", 
    "CEACAM3", "BST1", "PTPRC", "RAP2B", "CEACAM1", "RAP2C", "CEACAM8", "PTPRB", "CEACAM6", 
    "ALDOC", "ACTR10", "ALDOA", "RAB7A", "RAB5C", "SLC44A2", "HBB", "CTSZ", "PTPRJ", 
    "CTSS", "DSN1", "ACTR1B", "RAP1A", "RAP1B", "MAGT1", "NEU1", "ORMDL3", "MLEC", 
    "QSOX1", "CTSH", "CTSG", "CD36", "FCER1G", "CTSD", "CAMP", "CD33", "CTSC", "RAB6A", 
    "CAP1", "CD53", "CCT2", "DBNL", "MME", "SYK", "GAA", "CREG1", "QPCT", "ANO6", 
    "PGRMC1", "APRT", "SVIP", "RAB9B", "CD47", "DNASE1L3", "DNASE1L1", "DYNC1H1", "CD44", 
    "CD63", "CSTB", "ORM2", "ORM1", "GDI2", "FAF2", "CAB39", "AP2A2", "LILRA2", "LILRA3", 
    "CST3", "CXCR2", "CXCR1", "STBD1", "RAC1", "CD59", "ENPP4", "CCT8", "CD58", "SPTAN1", 
    "CD55", "LAIR1", "TBC1D10C", "CHRNB4", "CD300A", "GCA", "IRAG2", "PDAP1", "HUWE1", 
    "TYROBP", "ATAD3B", "OSCAR", "EEF1A1", "DOK3", "PPIA", "CD68", "DEGS1", "YPEL5", 
    "PPIE", "FCN1", "MPO", "ARHGAP45", "COTL1", "CAPN1", "OSTF1", "IDH1", "IGF2R", 
    "MGAM", "ATP6V0C", "MANBA", "CD93", "DNAJC13", "CYBB", "CYBA", "SDCBP", "CKAP4", 
    "DIAPH1", "TARM1", "SELL", "PRTN3", "VCL", "CRISP3", "ATP6V1D", "VCP", "C5AR1", 
    "COPB1", "SLC2A3", "PLAU", "PSEN1", "SLC2A5", "SRP14", "CNN2", "ADGRG3", "ALDH3B1", 
    "SIRPA", "ACAA1", "HGSNAT", "MMP8", "HSPA6", "CMTM6", "GGH", "HSPA8", "MMP9", 
    "ADA2", "DNASE1", "MAN2B1", "MS4A3", "DOCK2", "CHI3L1", "PSMD11", "B4GALT1", 
    "DDX3X", "PSMD13", "PSMD12", "PSMD14", "TNFAIP6", "MOSPD2", "HMGB1", "CXCL1", 
    "AGPAT2", "PLD1", "HK3", "SIRPB1", "HEBP2", "ADGRE3", "SYNGR1", "DPP7", "ADGRE5", 
    "TMEM30A", "PLAUR", "FGR", "GLB1", "CAT", "S100P", "DYNLT1", "PTX3", "LTF", 
    "CLEC12A", "MVP", "CAND1", "MAPK1", "GUSB"),
  
  "neutrophil_mediated_immunity"=c(
    "UNC13D", "ARMC8", "MCEMP1", "SCAMP1", "DYNC1LI1", "FLG2", "GALNS", "TOLLIP", "PSMC2", "PSMC3", 
    "HSPA1A", "GLA", "CDK13", "NHLRC3", "CDA", "HSPA1B", "CYFIP1", "PRSS2", "PDXK", "C3AR1", 
    "FRMPD3", "RETN", "ASAH1", "CRACR2A", "PSMD7", "PSMD6", "PSMD3", "ALOX5", "PSMD1", "PSMD2", 
    "STOM", "NCSTN", "PRSS3", "S100A12", "ELANE", "S100A11", "ATP6AP2", "PRKCD", "NCKAP1L", "HRNR", 
    "VAMP7", "VAMP8", "BIN2", "TCN1", "VAMP2", "PADI2", "CFD", "GMFG", "A1BG", "GNS", 
    "ABCA13", "CFP", "C3", "CYB5R3", "PRDX4", "PRDX6", "PLAC8", "NAPRT", "HMOX2", "LTA4H", 
    "ATG7", "FUCA1", "GSDMD", "ACE", "ANXA2", "ANXA3", "JUP", "CARD9", "FUCA2", "NFAM1", 
    "AZU1", "TMC6", "LYZ", "SIGLEC5", "SIGLEC9", "GPI", "PAFAH1B2", "SERPINA3", "SERPINA1", "FPR2", 
    "FPR1", "IQGAP1", "LRRC7", "IQGAP2", "TIMP2", "AP1M1", "DYNLL1", "SERPINB3", "AHSG", "ARG1", 
    "DGAT1", "PGAM1", "SERPINB1", "ADAM10", "STING1", "LRG1", "RNASE3", "RNASE2", "TUBB4B", "SERPINB6", 
    "ERP44", "ADAM17", "FABP5", "PSMA5", "PSMA2", "SLCO4C1", "LAMTOR1", "LAMTOR3", "PLEKHO2", "LAMTOR2", 
    "SERPINB12", "GRN", "PIGR", "PTAFR", "PYGB", "GSTP1", "HEXB", "SERPINB10", "PYGL", "LPCAT1", 
    "PSMB7", "PSMB1", "METTL7A", "AGA", "PGLYRP1", "OLR1", "GSN", "IST1", "AGL", "AOC1", 
    "DEFA4", "DEFA1", "ATP11B", "ATP11A", "ILF2", "TICAM2", "CLEC4C", "CLEC4D", "RAB31", "OLFM4", 
    "RAB37", "S100A7", "ARL8A", "ADAM8", "DSC1", "PTPN6", "GHDC", "PKP1", "S100A9", "FRK", 
    "S100A8", "CD177", "SLC27A2", "ROCK1", "HP", "EPX", "APEH", "NRAS", "PYCARD", "LAMP2", 
    "LAMP1", "RAB24", "KCMF1", "RAB4B", "TRAPPC1", "ATP8B4", "ARPC5", "KCNAB2", "EEF2", "C6ORF120", 
    "NFKB1", "ACLY", "DDOST", "HPSE", "LCN2", "FTL", "NIT2", "RAB5B", "KMT2E", "RAB3D", 
    "RAB3A", "ARHGAP9", "ATP8A1", "BRI3", "HVCN1", "GPR84", "RAB44", "STK10", "GLIPR1", "MMP25", 
    "CLEC5A", "TXNDC5", "TOM1", "GYG1", "TMEM179B", "SLC15A4", "ACTR2", "HSP90AA1", "APAF1", "PA2G4", 
    "NDUFC2", "TNFRSF1B", "TMBIM1", "DNAJC3", "FCGR2A", "NPC2", "DNAJC5", "DSG1", "PECAM1", "VPS35L", 
    "SNAP25", "SNAP23", "UBR4", "FCGR3B", "ALAD", "PSAP", "ACP3", "CEP290", "SNAP29", "CR1", 
    "DSP", "RNASET2", "CPPED1", "PKM", "CHIT1", "KPNB1", "NFASC", "PFKL", "RAB10", "TSPAN14", 
    "VAPA", "IMPDH1", "IMPDH2", "P2RX1", "RAB14", "MNDA", "BPI", "RAB18", "FOLR3", "ITGB2", 
    "CANT1", "HSP90AB1", "TCIRG1", "ITGAM", "FCAR", "ITGAL", "TTR", "GM2A", "ITGAX", "ITGAV", 
    "PGM2", "PGM1", "KRT1", "PTPRN2", "NBEAL2", "MAPK14", "CRISPLD2", "SLPI", "TMEM63A", "FGL2", 
    "STXBP2", "STXBP3", "CALML5", "PRCP", "PNP", "CPNE1", "CPNE3", "B2M", "ARSA", "ARSB", 
    "STK11IP", "SURF4", "XRCC5", "SLC11A1", "XRCC6", "COMMD3", "HLA-C", "MIF", "PPBP", "PRG3", 
    "RAB27A", "PRG2", "HLA-B", "VAT1", "VNN1", "CSNK2B", "CYSTM1", "MGST1", "TRPM2", "COMMD9", 
    "LGALS3", "ANPEP", "FTH1", "CD14", "GOLGA7", "ATP6V0A1", "CTSB", "CTSA", "RHOA", "PTGES2", 
    "SIGLEC14", "TUBB", "AMPD3", "LILRB2", "RHOF", "LILRB3", "RHOG", "NME2", "DERA", "TLR2", 
    "BST2", "CEACAM3", "BST1", "PTPRC", "RAP2B", "CEACAM1", "RAP2C", "CEACAM8", "PTPRB", "CEACAM6", 
    "ALDOC", "ACTR10", "ALDOA", "RAB7A", "RAB5C", "SLC44A2", "HBB", "CTSZ", "PTPRJ", "CTSS", 
    "DSN1", "ACTR1B", "RAP1A", "RAP1B", "MAGT1", "NEU1", "ORMDL3", "MLEC", "QSOX1", "CTSH", 
    "CTSG", "CD36", "FCER1G", "CTSD", "CAMP", "CD33", "CTSC", "RAB6A", "CAP1", "CD53", 
    "CCT2", "DBNL", "MME", "GAA", "CREG1", "QPCT", "ANO6", "PGRMC1", "APRT", "SVIP", 
    "RAB9B", "CD47", "DNASE1L1", "DYNC1H1", "CD44", "CD63", "CSTB", "ORM2", "ORM1", "GDI2", 
    "FAF2", "CAB39", "AP2A2", "LILRA3", "CST3", "CXCR2", "CXCR1", "STBD1", "RAC1", "CD59", 
    "ENPP4", "CCT8", "CD58", "SPTAN1", "CD55", "LAIR1", "TBC1D10C", "CHRNB4", "CD300A", "GCA", 
    "IRAG2", "PDAP1", "HUWE1", "TYROBP", "ATAD3B", "OSCAR", "EEF1A1", "DOK3", "PPIA", "CD68", 
    "DEGS1", "YPEL5", "PPIE", "FCN1", "MPO", "ARHGAP45", "COTL1", "CAPN1", "OSTF1", "IDH1", 
    "IGF2R", "MGAM", "ATP6V0C", "MANBA", "CD93", "DNAJC13", "CYBB", "CYBA", "SDCBP", "CKAP4", 
    "DIAPH1", "TARM1", "SELL", "PRTN3", "VCL", "CRISP3", "ATP6V1D", "VCP", "C5AR1", "COPB1", 
    "SLC2A3", "PLAU", "PSEN1", "SLC2A5", "SRP14", "CNN2", "ADGRG3", "ALDH3B1", "SIRPA", "ACAA1", 
    "HGSNAT", "MMP8", "HSPA6", "CMTM6", "GGH", "HSPA8", "IRAK4", "MMP9", "ADA2", "IL6", 
    "MAN2B1", "MS4A3", "DOCK2", "CHI3L1", "PSMD11", "B4GALT1", "DDX3X", "PSMD13", "PSMD12", "PSMD14", 
    "TNFAIP6", "MOSPD2", "HMGB1", "CXCL1", "AGPAT2", "PLD1", "HK3", "SIRPB1", "HEBP2", "ADGRE3", 
    "SYNGR1", "DPP7", "ADGRE5", "TMEM30A", "PLAUR", "FGR", "GLB1", "CAT", "S100P", "DYNLT1", 
    "PTX3", "LTF", "PLA2G1B", "CLEC12A", "MVP", "CAND1", "MAPK1", "GUSB"
  )
)


res_sham1vs_TBI14






logfc_cutoff <- 0.1
pval_cutoff <- 0.05

filtered_df <- as.data.frame(res_sham1vs_TBI14[, 1:5])
filtered_df <- filtered_df[!is.na(filtered_df$padj), ]
filtered_df <- filtered_df[abs(filtered_df$log2FoldChange) >= logfc_cutoff & filtered_df$padj < pval_cutoff, ]


filtered_df$log2FoldChange <- filtered_df$log2FoldChange * (-1)

sorted_df <- filtered_df[order(filtered_df$log2FoldChange, decreasing = TRUE), ]

sorted_df
rownames(sorted_df)



# Step 1: Create a named vector with log2FC values sorted by avg_log2FC
# Assuming your sorted_DEG_results has rownames as gene names and avg_log2FC as a column
ranked_genes <- sorted_df$log2FoldChange
names(ranked_genes) <- toupper(rownames(sorted_df))
ranked_genes

# Step 4: Run fgsea
fgsea_results <- fgsea(pathways = neutrophil_genelist, 
                       stats = ranked_genes, 
                       minSize = 10, 
                       maxSize = 500)  # Number of permutations, you can adjust this

# Step 5: View results
fgsea_results


ranked_df <- data.frame(
  Gene = names(ranked_genes),
  Correlation = ranked_genes
)




# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  # geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq( floor(round(min(ranked_genes),1)), ceiling(round(max(ranked_genes),1)), by = 3.5)) +
  scale_x_continuous(breaks = seq(0, 4000, by = 1000)) +
  labs(x = "Rank", y = "LogFC") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))


term = "Neutrophil Degranulation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_genes) %in% neutrophil_genelist[[term]])

# Enrichment plot 생성
enrichment_plot_cd8tem_a <- plotEnrichment(neutrophil_genelist[[term]], ranked_genes) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = term ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results$ES[fgsea_results$pathway == term ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results$padj[fgsea_results$pathway == term ] ,3 ) ,
                                                                               round( fgsea_results$NES[fgsea_results$pathway == term ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



grid.arrange(enrichment_plot_cd8tem_a,  bar_plot, scatter_plot, ncol = 1, heights = c(2.5, 0.3, 1.5))




term = "Neutrophil Degranulation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_genes) %in% neutrophil_genelist[[term]])

# Enrichment plot 생성
enrichment_plot_cd8tem_a <- plotEnrichment(neutrophil_genelist[[term]], ranked_genes) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = term ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results$ES[fgsea_results$pathway == term ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results$padj[fgsea_results$pathway == term ] ,3 ) ,
                                                                               round( fgsea_results$NES[fgsea_results$pathway == term ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



grid.arrange(enrichment_plot_cd8tem_a,  bar_plot, scatter_plot, ncol = 1, heights = c(2.5, 0.3, 1.5))






#### GSEA PDF Plot  ###################

cairo_pdf("/Users/nabinwon" , width = 5 , height = 4)

grid.arrange(enrichment_plot_cd8tem_a,  bar_plot, scatter_plot, ncol = 1, heights = c(2.5, 0.3, 1.5))

dev.off()










# Enrichment analysis ##########################

# 
# # 1. Install and load the necessary packages
# if (!requireNamespace("enrichR", quietly = TRUE)) {
#   install.packages("enrichR")
# }
# 
# if (!requireNamespace("writexl", quietly = TRUE)) {
#   install.packages("writexl")
# }

library(enrichR)
library(writexl)




enrichr_libraries <- listEnrichrDbs()

dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "Elsevier_Pathway_Collection")


ShamDay1_DEG 
TBIDay1_DEG

ShamDay14_DEG
TBIDay14_DEG


# Run enrichment for UpDEG
up_enrichment <- enrichr(rownames(ShamDay1_DEG), dbs)

# Run enrichment for DnDEG
dn_enrichment <- enrichr(rownames(TBIDay1_DEG), dbs)

# 5. Convert enrichment results to a list of data frames
up_enrichment_list <- lapply(up_enrichment, as.data.frame)
dn_enrichment_list <- lapply(dn_enrichment, as.data.frame)

# 6. Save the enrichment results to an Excel file
write_xlsx(list("ShamDay1_GO_BP" = up_enrichment_list$`GO_Biological_Process_2021`,
                "ShamDay1_KEGG" = up_enrichment_list$`KEGG_2021_Human`,
                "ShamDay1_MsigDB" = up_enrichment_list$`MSigDB_Hallmark_2020`,
                "ShamDay1_Elsevier" = up_enrichment_list$Elsevier_Pathway_Collection,
                "TBIDay1_GO_BP" = dn_enrichment_list$`GO_Biological_Process_2021`,
                "TBIDay1_KEGG" = dn_enrichment_list$`KEGG_2021_Human`,
                "TBIDay1_MsigDB" = dn_enrichment_list$MSigDB_Hallmark_2020,
                "TBIDay1_Elsevier" = dn_enrichment_list$Elsevier_Pathway_Collection),
           path = "/home/ch2860/RNA-seqCourse/Day1_Enrichr_Results.xlsx")




ShamDay14_DEG
TBIDay14_DEG


# Run enrichment for UpDEG
up_enrichment <- enrichr(rownames(ShamDay14_DEG), dbs)

# Run enrichment for DnDEG
dn_enrichment <- enrichr(rownames(TBIDay14_DEG), dbs)

# 5. Convert enrichment results to a list of data frames
up_enrichment_list <- lapply(up_enrichment, as.data.frame)
dn_enrichment_list <- lapply(dn_enrichment, as.data.frame)

# 6. Save the enrichment results to an Excel file
write_xlsx(list("ShamDay14_GO_BP" = up_enrichment_list$`GO_Biological_Process_2021`,
                "ShamDay14_KEGG" = up_enrichment_list$`KEGG_2021_Human`,
                "ShamDay14_MsigDB" = up_enrichment_list$`MSigDB_Hallmark_2020`,
                "ShamDay14_Elsevier" = up_enrichment_list$`Elsevier_Pathway_Collection`,
                "TBIDay14_GO_BP" = dn_enrichment_list$`GO_Biological_Process_2021`,
                "TBIDay14_KEGG" = dn_enrichment_list$`KEGG_2021_Human`,
                "TBIDay14_MsigDB" = dn_enrichment_list$`MSigDB_Hallmark_2020` ,
                "TBIDay14_Elsevier" = dn_enrichment_list$`Elsevier_Pathway_Collection`),
           path = "/home/ch2860/RNA-seqCourse/Day14_Enrichr_Results.xlsx")


getwd()



#### Fig5 TCR GSEA Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape  ###################
install.packages("Cairo")
cairo_pdf( sprintf ("/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/GSEA Plot/GSEA_KRAS.pdf") , width = 5 , height = 4)

grid.arrange(enrichment_plot_cd8tem_a,  bar_plot, scatter_plot, ncol = 1, heights = c(2.5, 0.3, 1.5))

dev.off()




#########################################################





BiocManager::install("impute")
library("WGCNA")
install.packages("WGCNA")
BiocManager::install("preprocessCore")

# 10. WGCNA ########################
library(WGCNA)
library(ggplot2)
library(gridExtra)

vsd -> dds_norm
norm.counts <- assay(dds_norm) %>% 
  t()


# normalized_counts -> norm.counts

power <- c(c(1:10), seq(from = 12, to = 80, by = 2))
sft <- pickSoftThreshold(norm.counts, powerVector = power, networkType = "signed", verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)




# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 9 # Please note grid.arrange(a1, a2, nrow = 2) result 
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)





module_eigengenes <- bwnet$MEs


# Print out a preview
# 각 모듈의 eigengene 값으로, 특정 모듈의 발현 패턴을 나타내며, 
# 샘플 간 비교를 통해 특정 모듈의 활성화 여부나 차이를 확인할 수 있음.
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
# grey module = all genes that doesn't fall into other modules were assigned to the grey module
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)






library(igraph)


# 'black' 모듈 유전자 추출
genes_black <- module_eigengenes$black
module_counts <- norm.counts[, genes_black, drop = FALSE]

# 유전자 간의 상관행렬 생성
adjacency <- abs(WGCNA::cor(module_counts))
adjacency[adjacency < 0.9] <- 0  # 필터링하여 낮은 상관 제거 (옵션)

diag(adjacency) <- 0

# igraph 객체 생성
network_black <- graph.adjacency(adjacency, mode = "directed", weighted = TRUE)

# 네트워크 시각화
plot(network_black, vertex.label = genes_black, main = "Black Module Network")

# 네트워크 시각화 (노드와 화살표 크기 줄이기)
plot(network_black, 
     vertex.label = genes_black, 
     vertex.size = 8,              # 노드 크기 줄이기
     vertex.label.cex = 0.7,       # 레이블 크기 줄이기
     layout = layout_with_fr,      # Fruchterman-Reingold 레이아웃 사용
     main = "Black Module Network",
     edge.width = 0.5,             # 간선 너비 조정
     edge.arrow.size = 0.1         # 화살표 크기 줄이기
)

##Dotplot#############################

library(ggplot2)
library(viridis)
up_enrichment_list$GO_Biological_Process_2021 -> data



# Adjusted.P.value로 정렬하여 상위 10개 선택
# 여기에서 특정 원하는 term 들만 subset 을 이용하여 데이터프레임을 구성하면 편함
top_data <- data[order(data$Adjusted.P.value), ][1:10, ]

# Overlap 열에서 분자 부분만 추출하고 Gene ratio 계산
top_data$Gene_Count <- as.numeric(sapply(strsplit(top_data$Overlap, "/"), function(x) as.numeric(x[1])))
top_data$Gene_Ratio <- top_data$Gene_Count / as.numeric(sapply(strsplit(top_data$Overlap, "/"), function(x) as.numeric(x[2])))

# Term 순서를 Gene_Ratio의 내림차순으로 설정
top_data$Term <- factor(top_data$Term, levels = top_data$Term[order(top_data$Gene_Ratio, decreasing = TRUE)])

# DotPlot 생성
ggplot(top_data, aes(x = Gene_Ratio, y = Term)) +
  geom_point(aes(size = Gene_Count, color = Adjusted.P.value)) +
  scale_color_viridis_c(option = "viridis") +
  scale_size_continuous(range = c(3, 10)) +  # 원의 크기 최소 3, 최대 10으로 설정
  scale_y_discrete(limits = rev(levels(top_data$Term))) +  # Y축을 역순으로 설정
  labs(
    x = "Gene ratio",
    y = "Term",
    color = "Adjusted P-value",
    size = "Gene Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

























