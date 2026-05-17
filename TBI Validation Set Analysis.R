#Validation Set Analysis

source("/Users/nabinwon/Desktop/RNA-seq/renv/activate.R")
library(DESeq2)
file_path <- ("/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/GSE163415_29DPI_samples_counts_table.txt")
file_path_2 <-("/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/GSE163415_3DPI_samples_counts_table.txt")

count_data_1 <- read.table(file_path, sep = "\t", header = TRUE)
count_data_2 <- read.table(file_path_2, sep = "\t", header = TRUE)
total_data<-cbind(count_data_1,count_data_2)
head(total_data)
dim(total_data)


total_dataframe <- as.data.frame(total_data)
write.csv(total_dataframe, file = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/validation_total_data.txt")

rownames(total_data) <- total_data$Gene

head(total_data)
total_data$Gene <- NULL
total_data_filtered<- total_data[, !grepl("Drug|Hypo|Thal", names(total_data))]
colnames(total_data_filtered)


ordered_cols <- c(
  # 3DPI Hippocampus
  "X3DPI.Hipp.Veh.NoTBI.1", "X3DPI.Hipp.Veh.NoTBI.2", "X3DPI.Hipp.Veh.NoTBI.3", "X3DPI.Hipp.Veh.NoTBI.4",
  "X3DPI.Hipp.Veh.TBI.1", "X3DPI.Hipp.Veh.TBI.2", "X3DPI.Hipp.Veh.TBI.3", "X3DPI.Hipp.Veh.TBI.4",
  
  # 29DPI Hippocampus
  "X29DPI.Hipp.Veh.NoTBI.1", "X29DPI.Hipp.Veh.NoTBI.2", "X29DPI.Hipp.Veh.NoTBI.3", "X29DPI.Hipp.Veh.NoTBI.4",
  "X29DPI.Hipp.Veh.TBI.1", "X29DPI.Hipp.Veh.TBI.2", "X29DPI.Hipp.Veh.TBI.3", "X29DPI.Hipp.Veh.TBI.4"
)

total_data_ordered <- total_data[, c(ordered_cols)]
colnames(total_data_ordered) <- gsub("\\.", "", ordered_cols)
colnames(total_data_ordered)

condition <- c(
  # 3DPI Hippocampus
  "NoTBI", "NoTBI", "NoTBI", "NoTBI",  # Sham
  "TBI", "TBI", "TBI", "TBI",          # TBI
  
  # 29DPI Hippocampus
  "NoTBI", "NoTBI", "NoTBI", "NoTBI",  # Sham
  "TBI", "TBI", "TBI", "TBI"          # TBI
)

day <- c(
  # 3DPI Hippocampus
  "3DPI", "3DPI", "3DPI", "3DPI",
  "3DPI", "3DPI", "3DPI", "3DPI",
  
  # 29DPI Hippocampus
  "29DPI", "29DPI", "29DPI", "29DPI",
  "29DPI", "29DPI", "29DPI", "29DPI"
)
conditionday <- paste(condition, day, sep = "")

coldata <- data.frame(
  condition = factor(condition),
  day = factor(day),
  conditionday = factor(conditionday)
)


rownames(coldata) <- colnames(total_data_ordered)

coldata

dds <- DESeqDataSetFromMatrix(countData = total_data_ordered, colData = coldata, design = ~ conditionday)
dds

smallestGroupSize <- ncol(total_data_ordered) * 0.5
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

plotPCA(vsd, intgroup=c("condition", "day"))

dds$conditionday <- relevel(dds$conditionday, ref = "NoTBI29DPI")
dds <- DESeq(dds)

res_day29_shamVSTBI <- results(dds, contrast = c("conditionday", "TBI29DPI", "NoTBI29DPI"))
resultsNames(dds)

validation_TBI <- as.data.frame(res_day29_shamVSTBI)
write.csv(validation_TBI, file = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/validaiton_TBI.csv")

res_day29_shamVSTBI_shrunk <- lfcShrink(dds, 
                                        coef = "conditionday_TBI29DPI_vs_NoTBI29DPI", 
                                        type = "apeglm")

head(res_day29_shamVSTBI_shrunk)
res_day29_shamVSTBI_shrunk

UpDEG <- subset(res_day29_shamVSTBI_shrunk, log2FoldChange > 0.5 & padj < 0.05)
DnDEG <- subset(res_day29_shamVSTBI_shrunk, log2FoldChange < -0.5 & padj < 0.05)

UpDEG
DnDEG

up_genes <- rownames(UpDEG)
dn_genes <- rownames(DnDEG)

up_genes
dn_genes


up_enrichment <- enrichr(up_genes, dbs)
dn_enrichment <- enrichr(dn_genes, dbs)

up_enrichment_list <- lapply(up_enrichment, as.data.frame)
dn_enrichment_list <- lapply(dn_enrichment, as.data.frame)

install.packages("writexl")
library(writexl)

write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "BioCarta_up" = up_enrichment_list$BioCarta_2016,
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/TBI_sham29vsTBI29validation.xlsx")

BiocManager::install("impute")
library("WGCNA")
install.packages("WGCNA")
BiocManager::install("preprocessCore")

# WGCNA ########################
library(WGCNA)
library(ggplot2)
library(gridExtra)

vsd -> dds_norm
norm.counts <- assay(dds_norm) %>% 
  t()

norm.counts
View(norm.counts)

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

soft_power <- 8 # Please note grid.arrange(a1, a2, nrow = 2) result 
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
genes_black <- names(bwnet$colors[bwnet$colors == "black"])
module_counts <- norm.counts[, genes_black, drop = FALSE]

saveRDS(bwnet, "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/validationset_totalWGCNA_pwr8.rds")

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


volcano_data <- as.data.frame(res_sham1vs_TBI14)
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


##Dotplot####################################
# DotPlot 생성
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







