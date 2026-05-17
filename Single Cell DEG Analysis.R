Idents(seurat_object) <- seurat_object$harmony_clusters

# sum(table(combined.integrated@meta.data$SCT_snn_res.0.2))
clusters <- seurat_object$harmony_clusters

cluster_names <- rep("Unknown", length(clusters))

cluster_names[clusters %in% c(0,3 ,5 ,15, 23)] <- "Microglia1"
cluster_names[clusters %in% c(2,6, 17)] <- "Microglia2"
cluster_names[clusters %in% c(1, 19)] <- "Neutrophil"
cluster_names[clusters %in% c(10)] <- "Neural_Progenitor"

cluster_names[clusters %in% c(20)] <- "B"
cluster_names[clusters %in% c(9)] <- "T"

cluster_names[clusters %in% c(4, 11,24)] <- "Ependymal"
cluster_names[clusters %in% c(7, 13,18, 22 )] <- "Stromal"
cluster_names[clusters %in% c(8, 12 )] <- "Astrocytes"
cluster_names[clusters %in% c(16)] <- "Neuron"
cluster_names[clusters %in% c(14,21)] <- "Oligodendrocyte"



seurat_object <- AddMetaData(seurat_object, metadata = cluster_names, col.name = "TBI_Cell_Annotation")
DimPlot(seurat_object, reduction = "umap.harmony",label =T, label.size = 8, repel = T, group.by = "harmony_clusters")
DimPlot(seurat_object, reduction = "umap.harmony",label =T, label.size = 8, repel = T, group.by = "TBI_Cell_Annotation")

saveRDS(seurat_object, "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/TBI_cell_annotation.rds")

Idents(seurat_object)<-seurat_object$Condition
table(seurat_object$Condition)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_object, ident.1 ='Naïve' , ident.2 ='CCI')
seurat_object<-()
head(cluster5.markers, n = 5)

rm(TBI_seurat_obj)


TBI_microglia1<-subset(seurat_object,TBI_Cell_Annotation %in% 'Microglia1')

Idents(TBI_microglia1)<-TBI_microglia1$Timepoint
table(TBI_microglia1$Timepoint)
cluster6.markers <- FindAllMarkers(TBI_microglia1,logfc.threshold = 1,min.pct = 0.1,only.pos = TRUE)
#cluster5.markers <- FindMarkers(TBI_microglia1, ident.1 ='6 month' , ident.2 ='24 hours',logfc.threshold = 0.5)

cluster6.markers[c("Ccr5","Cd44","Cd47","C3ar1","Csf1r","Cd8b1","Plxnb2"),]
cluster6.markers[cluster6.markers$gene %in% c("Ccr5","Cd44","Cd47","C3ar1","Csf1r","Cd8b1","Plxnb2"), ]

# Step 1: Subset DEGs for "7 days" time point
deg_7day <- cluster6.markers[cluster6.markers$cluster == "7 day", ]

# Step 2: Extract gene names into a vector
genes_7day <- deg_7day$gene

# Step 3: Display the first few gene names to verify
head(genes_7day)

# Optional: Save the gene vector to a file if needed

write.table(cluster6.markers,"/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/microglia1_deg.txt",quote = FALSE,sep = "\t",row.names = FALSE)


#####Neural_Progenitor_DEG########################################

TBI_neural_progenitor<-subset(seurat_object,TBI_Cell_Annotation %in% 'Neural_Progenitor')

Idents(TBI_neural_progenitor)<-TBI_neural_progenitor$Timepoint
table(TBI_neural_progenitor$Timepoint)
cluster7.markers <- FindMarkers(TBI_neural_progenitor, ident.1 ='6 month' , ident.2 ='24 hours',logfc.threshold = 0.5)
cluster7.markers$gene <- rownames(cluster7.markers)
cluster7.markers


write.table(cluster7.markers,"/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/neural_progenitor_deg_24hvs6months.txt",quote = FALSE,sep = "\t",row.names = FALSE)

library(enrichR)
library(writexl)

up_genes <- cluster7.markers$gene
dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "Elsevier_Pathway_Collection")
up_enrichment <- enrichr(up_genes, dbs)
up_enrichment_list <- lapply(up_enrichment, as.data.frame)

write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "BioCarta_up" = up_enrichment_list$BioCarta_2016,
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Neural_progenitor_enrichment_result.xlsx")



# Check and handle NULL values
up_enrichment_list <- lapply(up_enrichment_list, function(x) {
  if (is.null(x)) {
    data.frame()  # Replace NULL with an empty data frame
  } else {
    x
  }
})

# Write to Excel
write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "KEGG_up" = up_enrichment_list$KEGG_2021_Human, # Corrected name
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Neural_progenitor_enrichment_result.xlsx")











# Load required libraries
library(Seurat)
library(dplyr)

# Get top markers for each timepoint (adjust the number of markers per group if needed)
top_markers <- cluster7.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) # Adjust 'n' as desired

# Generate the heatmap
DoHeatmap(TBI_neural_progenitor, 
          features = top_markers$gene, 
          size = 3) + 
  theme(axis.text.y = element_text(size = 5)) # Adjust size for better visibility



# Extract all genes from cluster7.markers
all_genes <- cluster7.markers$gene

# Generate the heatmap with all genes
DoHeatmap(TBI_neural_progenitor, 
          features = all_genes, 
          size = 3) + 
  theme(axis.text.y = element_text(size = 5)) # Adjust size for better readability


pdf("/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Plots:Graphs/neural_progenitor_DEG.pdf", width = 10, height = 12)
DoHeatmap(TBI_neural_progenitor, features = all_genes, size = 3) +
  theme(axis.text.y = element_text(size = 5))
dev.off()

##########Microglia2##############################

TBI_microglia2<-subset(seurat_object,TBI_Cell_Annotation %in% 'Microglia2')

Idents(TBI_microglia2)<-TBI_microglia2$Timepoint
table(TBI_microglia2$Timepoint)
cluster8.markers <- FindAllMarkers(TBI_microglia2,logfc.threshold = 0.5,min.pct = 0.1,only.pos = TRUE)
cluster8.markers

write.table(cluster8.markers,"/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/microglia2_deg.txt",quote = FALSE,sep = "\t",row.names = FALSE)

library(ggplot2)
library(pheatmap)

all_genes <- cluster8.markers$gene

# Generate the heatmap with all genes
DoHeatmap(TBI_microglia2, 
          features = all_genes, 
          size = 3) + 
  theme(axis.text.y = element_text(size = 5)) # Adjust size for better readability


pdf("/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Plots:Graphs/microglia2_DEG.pdf", width = 10, height = 12)
DoHeatmap(TBI_microglia2, features = all_genes, size = 3) +
  theme(axis.text.y = element_text(size = 5))
dev.off()


library(enrichR)
library(writexl)

up_genes <- cluster8.markers$gene
dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "Elsevier_Pathway_Collection")
up_enrichment <- enrichr(up_genes, dbs)
up_enrichment_list <- lapply(up_enrichment, as.data.frame)

write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "BioCarta_up" = up_enrichment_list$BioCarta_2016,
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Microglia2_enrichment_result.xlsx")



# Check and handle NULL values
up_enrichment_list <- lapply(up_enrichment_list, function(x) {
  if (is.null(x)) {
    data.frame()  # Replace NULL with an empty data frame
  } else {
    x
  }
})

# Write to Excel
write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "KEGG_up" = up_enrichment_list$KEGG_2021_Human, # Corrected name
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Microglia2_enrichment_result.xlsx")













########Neutrophil################################

TBI_neutrophil<-subset(seurat_object,TBI_Cell_Annotation %in% 'Neutrophil')



Idents(TBI_neutrophil)<-TBI_neutrophil$Timepoint
table(TBI_neutrophil$Timepoint)
cluster9.markers <- FindAllMarkers(TBI_neutrophil,logfc.threshold = 1, min.pct = 0.1, only.pos = TRUE)
cluster9.markers <- FindMarkers(TBI_neutrophil, ident.1 ='6 month' , ident.2 ='24 hours',logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)

cluster9.markers

write.table(cluster9.markers,"/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/neutrophil_deg.txt",quote = FALSE,sep = "\t",row.names = FALSE)

all_genes <- cluster9.markers$gene

# Generate the heatmap with all genes
DoHeatmap(TBI_neutrophil, 
          features = all_genes, 
          size = 3) + 
  theme(axis.text.y = element_text(size = 5)) # Adjust size for better readability





TBI_neutrophil<-subset(seurat_object,TBI_Cell_Annotation %in% 'Neutrophil')

Idents(TBI_neutrophil)<-TBI_neutrophil$Timepoint
table(TBI_neutrophil$Timepoint)
cluster7.markers <- FindMarkers(TBI_neutrophil, ident.1 ='6 month' , ident.2 ='24 hours',logfc.threshold = 0.5)
cluster7.markers$gene <- rownames(cluster7.markers)
cluster7.markers


write.table(cluster7.markers,"/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/neutrophil_deg_24hvs6months.txt",quote = FALSE,sep = "\t",row.names = FALSE)

library(enrichR)
library(writexl)
UpDEG <- subset(cluster7.markers, avg_log2FC > 0.5 & padj < 0.05)
DnDEG <- subset(cluster7.markers, avg_log2FC < -0.5 & padj < 0.05)
up_genes <- 
  dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "Elsevier_Pathway_Collection")
up_enrichment <- enrichr(up_genes, dbs)
up_enrichment_list <- lapply(up_enrichment, as.data.frame)

write_xlsx(list("GOBP_up" = up_enrichment_list$`GO_Biological_Process_2021`,
                "MSigDB_up" = up_enrichment_list$MSigDB_Hallmark_2020,
                "BioCarta_up" = up_enrichment_list$BioCarta_2016,
                "Elsevier_up" = up_enrichment_list$Elsevier_Pathway_Collection),
           path = "/Users/nabinwon/Desktop/RNA-seq/TBI Analysis/Neural_progenitor_enrichment_result.xlsx")



# Check and handle NULL values
up_enrichment_list <- lapply(up_enrichment_list, function(x) {
  if (is.null(x)) {
    data.frame()  # Replace NULL with an empty data frame
  } else {
    x
  }
})


