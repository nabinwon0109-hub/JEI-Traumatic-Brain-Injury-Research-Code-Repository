library(patchwork)
library(Seurat)

TBI_seurat_obj  <- readRDS("/Users/nabinwon/Downloads/scTBI_seurat_object.rds")

names(table(TBI_seurat_obj$TBI_Cell_Annotation))

cell_colors <- c(
  "Astrocytes" = "#B2DDF7",  
  "B" = "#FFB347",           
  "Ependymal" = "#AFE9D1",  
  "Microglia1" = "#77DD77",  
  "Microglia2" = "#FF6961",  
  "Neural Progenitor" = "#FDFD96", 
  "Neuron" = "#C7CEEA",      
  "Neutrophil" = "#FFD1DC",  
  "Oligodendrocyte" = "#CB99C9", 
  "Stromal" = "#F49AC2",     
  "T" = "#C5E384",          
  "Unknown" = "#B0B0B0"     
)


TBI_seurat_obj$TBI_Cell_Annotation[TBI_seurat_obj$TBI_Cell_Annotation == "Neural_Progenitor"] <- "Neural Progenitor"

Idents(TBI_seurat_obj) <- TBI_seurat_obj$TBI_Cell_Annotation
TBI_umap <- DimPlot(TBI_seurat_obj ,label=T , label.size = 6, raster = F,  reduction = "umap.harmony" ,  repel = T,
                    cols = cell_colors )+ ggtitle(NULL) + xlim(-16, 15) + ylim(-18, 12)  + NoAxes()  +theme(legend.text = element_text(margin = margin(t = 5, b = 5)) )



TBI_umap <- TBI_umap +
  # X?? ȭ??ǥ
  annotate("segment", x = -15, xend = -10, y = -17, yend = -17, 
           arrow = arrow(length = unit(0.2, "cm")), size = 0.8) +  # ȭ??ǥ
  annotate("text", x = -10, y = -17.5, label = "umap1", size = 5, vjust = 1, hjust=1) +  # X?? ???̺?
  
  # Y?? ȭ??ǥ
  annotate("segment", x = -15, xend = -15, y = -17, yend = -13, 
           arrow = arrow(length = unit(0.2, "cm")), size = 0.8) +  # ȭ??ǥ
  annotate("text", x = -15.5, y = -13, label = "umap2", size = 5, hjust = 1, vjust=0 , angle= 90)     # Y?? ???̺?

TBI_umap
