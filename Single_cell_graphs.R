# Neutrophil single Cell ##################################


library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)
library(tidyverse)

set.seed(1234)
options(future.globals.maxSize = 1e9)




# File Download and Preprocessing for Loading ##########################################
## File Extraction/Unzip and Renaming  ##########################################


# 압축 파일 경로와 작업할 디렉토리 설정
tar_file <- "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/raw/GSE269748_RAW.tar"  # tar 파일 경로
output_dir <- "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/raw/GSE269748_RAW"     # 압축을 해제할 디렉토리

# 1. tar 파일 압축 해제
untar(tar_file, exdir = output_dir)

# 2. 압축 해제된 파일 목록 확인
file_list <- list.files(output_dir, full.names = TRUE)

# 3. 각 파일을 Sample ID별 폴더로 이동
for (file_path in file_list) {
  # 파일 이름에서 Sample ID 추출 (예: GSM8326331)
  sample_id <- str_extract(basename(file_path), "^GSM[0-9]+")
  
  # Sample ID 폴더 생성 (폴더가 없으면 생성)
  sample_dir <- file.path(output_dir, sample_id)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir)
  }
  
  # 파일을 Sample ID 폴더로 이동
  file.rename(file_path, file.path(sample_dir, basename(file_path)))
}

# 4. 파일 이름 변경 함수 정의
rename_files <- function(folder_path) {
  # 폴더 내 파일 목록
  files <- list.files(folder_path, full.names = TRUE)
  
  # 각 파일에 대해 이름을 수정
  for (file_path in files) {
    if (grepl("features", basename(file_path))) {
      file.rename(file_path, file.path(folder_path, "features.tsv.gz"))
    } else if (grepl("matrix", basename(file_path))) {
      file.rename(file_path, file.path(folder_path, "matrix.mtx.gz"))
    } else if (grepl("barcodes", basename(file_path))) {
      file.rename(file_path, file.path(folder_path, "barcodes.tsv.gz"))
    }
  }
}

# 5. 각 Sample ID 폴더에서 파일 이름 변경 적용
sample_folders <- list.dirs(output_dir, recursive = FALSE)  # Sample ID 폴더 목록 가져오기
for (folder in sample_folders) {
  rename_files(folder)
}


# 아래의 형태로 파일을 만드려고 위의 코드들을 실행했어요
# GSE269748_RAW/
#   ├── GSM8326331/
#   │   ├── features.tsv.gz
#   │   ├── matrix.mtx.gz
#   │   └── barcodes.tsv.gz
#   ├── GSM8326332/
#   │   ├── features.tsv.gz
#   │   ├── matrix.mtx.gz
#   │   └── barcodes.tsv.gz



## Metadata ##########################################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269748 해당 사이트의 글을 참고하여 메타데이터를 만들었어요

sample_metadata <- data.frame(
  Sample_ID = c("GSM8326331", "GSM8326332", "GSM8326333", "GSM8326334", "GSM8326335", "GSM8326336", 
                "GSM8326337", "GSM8326338", "GSM8326339", "GSM8326340", "GSM8326341", "GSM8326342", 
                "GSM8326343", "GSM8326344", "GSM8326345", "GSM8326346", "GSM8326347", "GSM8326348", 
                "GSM8326349", "GSM8326350", "GSM8326351", "GSM8326352", "GSM8326353", "GSM8326354", 
                "GSM8326355", "GSM8326356", "GSM8326357", "GSM8326358", "GSM8326359", "GSM8326360", 
                "GSM8326361", "GSM8326362", "GSM8326363", "GSM8326364", "GSM8326365", "GSM8326366"),
  Timepoint = c(rep("24 hours", 28), rep("6 month", 4), rep("7 day", 4)),
  Sex = c(rep("male", 19), rep("female", 3), rep("male", 3),  rep("female", 3), rep("male", 8) ),
  Condition = c("HSCCI", "CCI", "Naïve", "CCI", "HSCCI", "Naïve", "rCHI", "Naïve", "rCHI", "rCHI", 
                "HSCCI", "CCI", "Naïve", "CCI", "CCI", "CCI", "CCI", "CCI", "CCI", "CCI", "CCI", 
                "CCI", "Naïve", "Naïve", "Naïve", "Naïve", "Naïve", "Naïve", "CCI", "CCI", "CCI", 
                "CCI", "CCI", "CCI", "CCI", "CCI"),
  Region = c("Adjacent", "ipsilateral", "Adjacent", "ipsilateral", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "ipsilateral", 
             "Adjacent", "ipsilateral", "contralateral", "ipsilateral", "contralateral", "ipsilateral", 
             "contralateral", "ipsilateral", "ipsilateral", "ipsilateral", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", 
             "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent", "Adjacent") ,
  Replicate = c("24h_male_HSCCI_1", "24h_male_CCI_ipsi_1", "24h_male_Naïve_1", "24h_male_CCI_ipsi_2", 
                "24h_male_HSCCI_2", "24h_male_Naïve_2", "24h_male_rCHI_1", "24h_male_Naïve_3", 
                "24h_male_rCHI_2", "24h_male_rCHI_3", "24h_male_HSCCI_3", "24h_male_CCI_ipsi_3", 
                "24h_male_Naïve_4", "24h_male_CCI_ipsi_4", "24h_male_CCI_contra_1", "24h_male_CCI_ipsi_5", 
                "24h_male_CCI_contra_2", "24h_male_CCI_ipsi_6", "24h_male_CCI_contra_3", "24h_female_CCI_ipsi_1", 
                "24h_female_CCI_ipsi_2", "24h_female_CCI_ipsi_3", "24h_male_Naïve_5", "24h_male_Naïve_6", 
                "24h_male_Naïve_7", "24h_female_Naïve_1", "24h_female_Naïve_2", "24h_female_Naïve_3", 
                "6m_male_CCI_1", "6m_male_CCI_2", "6m_male_CCI_3", "6m_male_CCI_4", "7d_male_CCI_1", 
                "7d_male_CCI_2", "7d_male_CCI_3", "7d_male_CCI_4")
)



## Select and Load Sample data  ##########################################

filtered_samples <- sample_metadata %>%
  filter(Sex == "male", Region != "contralateral", Condition != "rCHI") %>%
  filter(!Replicate %in% c("24h_male_Naïve_5","24h_male_Naïve_6","24h_male_Naïve_7","24h_male_CCI_ipsi_5","24h_male_CCI_ipsi_6", "6m_male_CCI_4" , "7d_male_CCI_4"))

nrow(filtered_samples)


# 파일 경로와 Sample ID에 맞게 Read10X 함수로 데이터 불러오기 및 Barcode 수정
for (i in 1:nrow(filtered_samples)) {
  sample_id <- filtered_samples$Sample_ID[i]
  file_path <- file.path("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/raw/GSE269748_RAW", sample_id)
  
  # 파일이 존재할 경우 Read10X로 불러오기
  if (dir.exists(file_path)) {
    sample_data <- Read10X(data.dir = file_path)
    
    # 기존 -1 접미사를 -i로 변경
    colnames(sample_data) <- gsub("-1$", paste0("-", i), colnames(sample_data))
    
    # 데이터 할당
    assign(sample_id, sample_data)
    cat(paste("Loaded and modified Barcode for:", sample_id, "nCell :", ncol(sample_data), "nGene :", nrow(sample_data),  "\n"))
  } else {
    cat(paste("File not found for:", sample_id, "\n"))
  }
}


## Metadata Processing  ##########################################

formatted_metadata <- data.frame()

# 각 Sample_ID에 대해 반복 수행
for (i in 1:nrow(filtered_samples)) {
  sample_id <- filtered_samples$Sample_ID[i]
  
  # 각 샘플의 데이터 가져오기
  sample_data <- get(sample_id)  # `assign`으로 생성한 오브젝트 가져오기
  
  # 세포 수에 맞게 sample_metadata 복제
  num_cells <- ncol(sample_data)
  sample_meta <- filtered_samples[i, ]
  
  # 각 샘플 메타데이터 복제하여 행 추가
  expanded_meta <- sample_meta[rep(1, num_cells), ]
  
  # 각 세포의 바코드 정보 추가
  expanded_meta$Barcode <- colnames(sample_data)
  
  # 전체 메타데이터에 결합
  formatted_metadata <- rbind(formatted_metadata, expanded_meta)
}

rownames(formatted_metadata) <- formatted_metadata$Barcode
formatted_metadata$Barcode <- NULL



## Count Matrix Combine ###############################

sample_list <- lapply(filtered_samples$Sample_ID, get)

# 모든 샘플 간 공통되는 rownames를 찾기
common_rownames <- Reduce(intersect, lapply(sample_list, rownames))
length(common_rownames)

# 각 샘플에서 공통되는 rownames만 선택하여 데이터 결합
sample_list_common <- lapply(sample_list, function(x) x[common_rownames, ])

# 공통 rownames만 남긴 상태에서 cbind로 결합하여 GSM_combined 선언
GSM_combined <- do.call(cbind, sample_list_common)





# Seurat Object Creating  ##########################################
TBI_seurat_obj <- CreateSeuratObject(counts = GSM_combined, project = "TBI", min.cells = 3, min.features = 200, meta.data = formatted_metadata)

# 각 샘플 별 세포의 수
table(TBI_seurat_obj$Sample_ID)


##object 안에 뭐가 있는지 보는 방법들 (용도에 맞게 사용)
# TBI_seurat_obj[["SampleID"]]
# TBI_seurat_obj@meta.data$SampleID
# TBI_seurat_obj$SampleID

#QC
TBI_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(TBI_seurat_obj, pattern = "^mt-")
VlnPlot(TBI_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# TBI_seurat_obj$nFeature_RNA

output_folder_name <- "StandardAnalysis"  # 원하는 폴더 이름으로 변경 가능
output_dir <- file.path("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result", output_folder_name)

# 지정한 폴더가 존재하지 않으면 생성
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# VlnPlot 생성
vln_plot <- VlnPlot(TBI_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3 , pt.size = 0 )

# ggsave를 사용하여 높은 해상도(예: 300dpi)로 저장
output_path <- file.path(output_dir, "1_VlnPlot_TBI_beforeQC.png")
ggsave(output_path, plot = vln_plot, dpi = 300, width = 12, height = 6)



#cutoff 정하기 어려우면 통계적 방식 접근해서 정함
median(TBI_seurat_obj$nFeature_RNA)
mean(TBI_seurat_obj$nFeature_RNA)
sd(TBI_seurat_obj$nFeature_RNA)
mean(TBI_seurat_obj$nFeature_RNA) + 3 * sd(TBI_seurat_obj$nFeature_RNA)

median(TBI_seurat_obj$percent.mt)
mean(TBI_seurat_obj$percent.mt)
sd(TBI_seurat_obj$percent.mt)
mean(TBI_seurat_obj$percent.mt) + 2 * sd(TBI_seurat_obj$percent.mt)

## Sample QC  ##########################################
TBI_seurat_obj <- subset(TBI_seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
table(TBI_seurat_obj$Sample_ID)

vln_plot <- VlnPlot(TBI_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3  , pt.size = 0 )

output_path <- file.path(output_dir, "2_VlnPlot_TBI_afterQC.png")
ggsave(output_path, plot = vln_plot, dpi = 300, width = 12, height = 6)



# Integration의 결과를 판단하기 위하여 Auto annotation tool인 Azimuth 실행 
# reference 비워둔 채로 실행 후 제시해주는 것 중 사용, 여기서는 그냥 pbmcref사용
TBI_seurat_obj <- RunAzimuth(TBI_seurat_obj, reference = "mousecortexref")

saveRDS(TBI_seurat_obj, "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/1_TBI_scData_notintegrated.rds")

TBI_seurat_obj <- readRDS("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/1_TBI_scData_notintegrated.rds")





# layer를 split함으로써 object는 하나로 유지
# 이때 normalization and variable feature identification is performed for each batch independently
TBI_seurat_obj[["RNA"]] <- split(TBI_seurat_obj[["RNA"]], f = TBI_seurat_obj$Sample_ID)
TBI_seurat_obj

# Normalization
TBI_seurat_obj <- NormalizeData(TBI_seurat_obj)
TBI_seurat_obj <- FindVariableFeatures(TBI_seurat_obj)
TBI_seurat_obj <- ScaleData(TBI_seurat_obj)
TBI_seurat_obj <- RunPCA(TBI_seurat_obj)
ElbowPlot(TBI_seurat_obj, ndims = 30)

dimensionNum <- 22

output_path <- file.path(output_dir, "3_ElbowPlot_TBI.png")
ggsave(output_path, plot = ElbowPlot(TBI_seurat_obj, ndims = 30), dpi = 300, width = 6, height = 6)



TBI_seurat_obj <- FindNeighbors(TBI_seurat_obj, dims = 1:dimensionNum, reduction = "pca")
TBI_seurat_obj <- FindClusters(TBI_seurat_obj, resolution = 0.5 , cluster.name = "unintegrated_clusters")


TBI_seurat_obj <- RunUMAP(TBI_seurat_obj, dims = 1:dimensionNum, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
# Dimplot이 dimensional reduction 함수
DimPlot(TBI_seurat_obj, reduction = "umap.unintegrated", group.by = c("Sample_ID", "predicted.subclass"))

output_path <- file.path(output_dir, "4_Unintegrated_UMAP_TBI.png")
ggsave(output_path, plot = DimPlot(TBI_seurat_obj, reduction = "umap.unintegrated", group.by = c("Sample_ID", "predicted.subclass")), dpi = 300, width = 15, height = 6)

saveRDS(TBI_seurat_obj, "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/1_TBI_scData_beforeintegrated.rds")
TBI_seurat_obj <- readRDS("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/1_TBI_scData_beforeintegrated.rds")





# TBI_seurat_objtmp <- subset (TBI_seurat_obj , Sample_ID %in%   c("GSM8326331", "GSM8326332", "GSM8326333", "GSM8326334") )
# 
# TBI_seurat_objtmp <- IntegrateLayers(
#   object = TBI_seurat_objtmp, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "harmony",
#   verbose = FALSE
# )
# TBI_seurat_objtmp <- IntegrateLayers(
#   object  = TBI_seurat_objtmp, method = CCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.cca",
#   verbose = FALSE
# )

# 여기까지는  visualize the results of a standard analysis without integration

## Sample Integration  ##########################################

# Integration
# 이제 여러 integration tool 사용 
TBI_seurat_obj <- IntegrateLayers(
  object = TBI_seurat_obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

# we can now visualize and cluster the datasets - harmony
TBI_seurat_obj <- FindNeighbors(TBI_seurat_obj, reduction = "harmony", dims = 1:dimensionNum)
TBI_seurat_obj <- FindClusters(TBI_seurat_obj, resolution = 0.5 , cluster.name = "harmony_clusters")
TBI_seurat_obj <- RunUMAP(TBI_seurat_obj, reduction = "harmony", dims = 1:dimensionNum, reduction.name = "umap.harmony")


TBI_seurat_obj <- IntegrateLayers(
  object  = TBI_seurat_obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# we can now visualize and cluster the datasets - cca
TBI_seurat_obj <- FindNeighbors(TBI_seurat_obj, reduction = "integrated.cca", dims = 1:dimensionNum)
TBI_seurat_obj <- FindClusters(TBI_seurat_obj, resolution = 0.5, cluster.name = "cca_clusters")
TBI_seurat_obj <- RunUMAP(TBI_seurat_obj, reduction = "integrated.cca", dims = 1:dimensionNum, reduction.name = "umap.cca")



a <- DimPlot(
  TBI_seurat_obj,
  reduction = "umap.unintegrated",
  group.by = c("Sample_ID"),
  combine = FALSE, label.size = 2
)

b <-  DimPlot(
  TBI_seurat_obj,
  reduction = "umap.harmony",
  group.by = c("Sample_ID"),
  combine = FALSE, label.size = 2
)


c <-  DimPlot(
  TBI_seurat_obj,
  reduction = "umap.cca",
  group.by = c("Sample_ID"),
  combine = FALSE, label.size = 2
)

a[[1]] | b[[1]] | c[[1]]


TBI_seurat_obj <- JoinLayers(TBI_seurat_obj)

saveRDS(TBI_seurat_obj, "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/2_TBI_scData_integrated.rds")

TBI_seurat_obj <- JoinLayers(TBI_seurat_obj)

output_path <- file.path(output_dir, "5_integrated_UMAP_TBI.png")
ggsave(output_path, plot = (a[[1]] | b[[1]] | c[[1]]) , dpi = 300, width = 18, height = 5)






# annotation ########################################
## Automated annotation tools #######################
library(celldex)
library(SingleR)
library(SingleCellExperiment)

as.SingleCellExperiment(TBI_seurat_obj) -> TBI_seurat_obj_SCE

# celldex::
celldex::ImmGenData() -> ImmGen
celldex::MouseRNAseqData() -> MouseRNAseq


TBI_seurat_obj_SCE_ImmGen <- SingleR(test = TBI_seurat_obj_SCE, ref = ImmGen, assay.type.test=1,
                                     labels = ImmGen$label.main)

TBI_seurat_obj_SCE_MouseRNAseq <- SingleR(test = TBI_seurat_obj_SCE, ref = MouseRNAseq, assay.type.test=1,
                                          labels = MouseRNAseq$label.main)

TBI_seurat_obj$ImmGen <-TBI_seurat_obj_SCE_ImmGen$labels 

DimPlot( TBI_seurat_obj, group.by = "ImmGen", reduction = "umap.harmony" , label=T)


TBI_seurat_obj$MouseRNAseq <-TBI_seurat_obj_SCE_MouseRNAseq$labels 

DimPlot( TBI_seurat_obj, group.by = "MouseRNAseq", reduction = "umap.harmony" , label=T)




library(HGNChelper)

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# get cell-type-specific gene sets from our in-built database (DB)
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(Tonly_CRC_filtered[["RNA"]]$scale.data) else as.matrix(Tonly_CRC_filtered[["RNA"]]@scale.data)
scRNAseqData_scaled <- as.matrix(TBI_seurat_obj[["RNA"]]$scale.data)


es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# View(es.max)


# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(TBI_seurat_obj@meta.data$cca_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(TBI_seurat_obj@meta.data[TBI_seurat_obj@meta.data$cca_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(TBI_seurat_obj@meta.data$cca_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])



TBI_seurat_obj@meta.data$sctype_classification_ImmuneCCA = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  TBI_seurat_obj@meta.data$sctype_classification_ImmuneCCA[TBI_seurat_obj@meta.data$cca_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(TBI_seurat_obj, reduction = "umap.cca", label = TRUE, repel = TRUE, group.by = 'sctype_classification_ImmuneCCA') 




cL_resutls <- do.call("rbind", lapply(unique(TBI_seurat_obj@meta.data$harmony_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(TBI_seurat_obj@meta.data[TBI_seurat_obj@meta.data$harmony_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(TBI_seurat_obj@meta.data$harmony_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])



TBI_seurat_obj@meta.data$sctype_classification_ImmuneHarmony = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  TBI_seurat_obj@meta.data$sctype_classification_ImmuneHarmony[TBI_seurat_obj@meta.data$harmony_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(TBI_seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = 'sctype_classification_ImmuneHarmony')       



gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Brain") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(Tonly_CRC_filtered[["RNA"]]$scale.data) else as.matrix(Tonly_CRC_filtered[["RNA"]]@scale.data)
scRNAseqData_scaled <- as.matrix(TBI_seurat_obj[["RNA"]]$scale.data)


es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# View(es.max)


# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(TBI_seurat_obj@meta.data$cca_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(TBI_seurat_obj@meta.data[TBI_seurat_obj@meta.data$cca_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(TBI_seurat_obj@meta.data$cca_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])



TBI_seurat_obj@meta.data$sctype_classification_BrainCCA = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  TBI_seurat_obj@meta.data$sctype_classification_BrainCCA[TBI_seurat_obj@meta.data$cca_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(TBI_seurat_obj, reduction = "umap.cca", label = TRUE, repel = TRUE, group.by = 'sctype_classification_BrainCCA') 




cL_resutls <- do.call("rbind", lapply(unique(TBI_seurat_obj@meta.data$harmony_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(TBI_seurat_obj@meta.data[TBI_seurat_obj@meta.data$harmony_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(TBI_seurat_obj@meta.data$harmony_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])



TBI_seurat_obj@meta.data$sctype_classification_BrainHarmony = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  TBI_seurat_obj@meta.data$sctype_classification_BrainHarmony[TBI_seurat_obj@meta.data$harmony_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(TBI_seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE, group.by = 'sctype_classification_BrainHarmony')       









## Manual annotation #######################
### Harmony ############################

DimPlot(TBI_seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE,  label.size = 6,group.by = 'sctype_classification_ImmuneHarmony') 
DimPlot(TBI_seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE, label.size = 6, group.by = 'sctype_classification_BrainHarmony') 
DimPlot(TBI_seurat_obj, reduction = "umap.harmony", label = TRUE, repel = TRUE, label.size = 6, group.by = 'Sample_ID') 


Microglia_marker <- c("Ccl4","Lyz2","Tmem119","Enpp2")
Neutrophil_marker <- c("Retnlg","Lyz2","Ccl4","Ifitm3")

B_marker <- c("Cd74","Ifitm3","H2-Eb1","Lyz2","Ccl4")
T_marker <- c("Cd2","Cd69","Cd74","H2-Eb1","Ccl4")
NK_marker <- c("Cd2","Cd69","Ccl4","Enpp2")

Ependymal_marker <- c("Calml4","Enpp2","Clu")
Endothelial_marker <- c("Spock2","Cldn5","Ifitm3")

Astrocytes_marker <- c("Clu","Slc1a2","Aldoc","Aqp4","Gpr37l1","Ptprz1")
Neuron_marker <- c("Stmn2","Tubb3","Sox11","Ptprz1")
Oligodendrocytes_marker <- c("Opalin","Plp1","Mog")

all_markers_union <- unique(c(
  Microglia_marker, Neutrophil_marker, B_marker, T_marker, NK_marker,
  Ependymal_marker, Endothelial_marker, Astrocytes_marker, Neuron_marker, Oligodendrocytes_marker
))



microglia <- c("Cx3cr1", "P2ry12",  "Tmem119")
neutrophil <- c("Retnlg","Ly6g", "Mpo", "S100a9", "Cxcr2")
b_cell <- c("Cd19", "Cd79a", "Ms4a1", "Cd22")
t_cell <- c("Cd3e", "Cd4", "Cd8a", "Trac")
nk_cell <- c("Ncr1", "Klrb1c", "Klra8", "Gzmb")
ependymal <- c("Foxj1", "Dcxr", "S100b", "Ttr","Pifo","Dynlrb2","Gfap","Foxj1" ,"Tuba1a","Rarres2")
endothelial <- c("Spock2","Cldn5","Pecam1", "Cdh5", "Vwf", "Nos3")
astrocytes <- c("Gfap", "Aldh1l1", "S100b", "Slc1a3","Slc1a2","Aqp4")
neuron <- c("Rbfox3", "Map2", "Tubb3", "Syp","Stmn2")
oligodendrocytes <- c("Mog", "Olig1", "Cnp", "Plp1","Opalin")

# Combine all markers into a single vector (union)
all_markers <- unique(c(
  microglia, neutrophil, b_cell, t_cell, nk_cell,
  ependymal, endothelial, astrocytes, neuron, oligodendrocytes
))


Idents(TBI_seurat_obj) <- factor(TBI_seurat_obj$harmony_clusters, levels = 0:25)
VlnPlot(TBI_seurat_obj ,features =  all_markers, flip = T, stack = T) + theme(legend.position = "none")
DimPlot(TBI_seurat_obj, reduction = "umap.harmony", group.by = "predicted.subclass", label =T, label.size = 8, repel = T)

VlnPlot(TBI_seurat_obj, features =  all_markers_union, flip = T, stack = T)+ theme(legend.position = "none")

DimPlot(TBI_seurat_obj, reduction = "umap.harmony",label =T, label.size = 8, repel = T)



Bulk_intersected_DEGs <- c("Aspg", "C1qa", "C1qb", "C1qc", "C3", "C3ar1", "C4b", "Capg", "Cd14", "Cd22", "Cd44", 
                           "Cd48", "Cd68", "Cd84", "Cmklr1", "Csf3r", "Ctsh", "Ctss", "Cyba", "Cybb", 
                           "Fcgr2b", "Flnc", "Gbp2", "Gfap", "Gpnmb", "Itgb2", "Lilrb4a", "Lyz2", "Mpeg1", 
                           "Myo1f", "Ptprc", "Slc11a1", "Tyrobp", "Unc93b1", "Vim")

# C4a -> C4b
# Lilrb4  ->  Lilrb4a
# Clec7a -> ??? not detected


VlnPlot(TBI_seurat_obj, features =  Bulk_intersected_DEGs , flip = T, stack = T)+ theme(legend.position = "none")


table(TBI_seurat_obj$Sample_ID, TBI_seurat_obj$harmony_clusters )
table(TBI_seurat_obj$Condition, TBI_seurat_obj$harmony_clusters )
table(TBI_seurat_obj$Replicate, TBI_seurat_obj$harmony_clusters )
prop.table(table(TBI_seurat_obj$Replicate, TBI_seurat_obj$harmony_clusters ), margin =1)



Idents(TBI_seurat_obj) <- factor(TBI_seurat_obj$harmony_clusters, levels = 0:25)
TBI.markers <- FindAllMarkers(TBI_seurat_obj, only.pos = TRUE)
TBI.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_marker

top10_marker

write.table(top10_marker, "/home/Data_Drive_8TB_3/ch2860/RNAcourse/TBI_Top10marker.txt", sep = "\t",quote = F)

FeaturePlot(TBI_seurat_obj , features =  c("Tmem119","Ifitm3","Cd74","H2-Eb1"), ncol = 2 , reduction = "umap.harmony")







Microglia1_Seunum <- c(0,3 ,5 ,15, 23 )
Microglia2_Seunum <- (2,6)
Neutrophil_Seunum <- c(1, 19 ) # 17
Neural_Progenitor_Seunum <- c( 10 )

B_Seunum <- c(20)
T_Seunum <- c(9)
NK_Seunum <- c()

Ependymal_Seunum <- c(4, 11,24)   # 24 Ciliated cell
Stromal_Seunum <-c(7, 13,18, 22 )  # Endothelial 7 13 18  # 22 Fibroblast
Astrocytes_Seunum <-c(8, 12 )
Neuron_Seunum <- c( 16)
Oligodendrocytes_Seunum <- c(14, 21)


Unknown <- c(25 )



Idents(TBI_seurat_obj) <- TBI_seurat_obj$harmony_clusters

# sum(table(combined.integrated@meta.data$SCT_snn_res.0.2))
clusters <- TBI_seurat_obj$harmony_clusters
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



TBI_seurat_obj <- AddMetaData(TBI_seurat_obj, metadata = cluster_names, col.name = "TBI_Cell_Annotation")

DimPlot(TBI_seurat_obj, reduction = "umap.harmony",label =T, label.size = 8, repel = T, group.by = "TBI_Cell_Annotation")

saveRDS(TBI_seurat_obj, "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_seurat_object.rds")




### CCA ############################
DimPlot(TBI_seurat_obj, reduction = "umap.cca", label = TRUE, repel = TRUE, group.by = 'sctype_classification_ImmuneCCA') 
DimPlot(TBI_seurat_obj, reduction = "umap.cca", label = TRUE, repel = TRUE, group.by = 'sctype_classification_BrainCCA') 


Idents(TBI_seurat_obj) <- factor(TBI_seurat_obj$cca_clusters, levels = 0:22)
DimPlot(TBI_seurat_obj, reduction = "umap.cca",label =T, label.size = 8, repel = T)

DimPlot(TBI_seurat_obj, reduction = "umap.cca", group.by = "predicted.subclass", label =T, label.size = 8, repel = T)


VlnPlot(TBI_seurat_obj ,features =  all_markers, flip = T, stack = T) + theme(legend.position = "none")

VlnPlot(TBI_seurat_obj, features =  all_markers_union, flip = T, stack = T)+ theme(legend.position = "none")


2     5     9  11  13  19 20  22

Microglia_Seunum <- c(0, 3 ,7 )
Neutrophil_Seunum <- c(1, 8)

B_Seunum <- c(18)
T_Seunum <- c(12)
NK_Seunum <- c()

Ependymal_Seunum <- c(6 ,10 , 14 ,21 )
Endothelial_Seunum <-c(4,15 , 17  )

Astrocytes_Seunum <-c()
Neuron_Seunum <- c()
Oligodendrocytes_Seunum <- c(16)


# CiberSortX signature Mat  #############################



DefaultAssay(TBI_seurat_obj) <- "RNA"
Idents(TBI_seurat_obj) <- TBI_seurat_obj$TBI_Cell_Annotation

TBI_only <- subset(TBI_seurat_obj, Condition %in% c("CCI","HSCCI")  )


# DEG 추출 (클러스터링된 그룹 간 비교)
deg_markers <- FindAllMarkers(TBI_only, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

# DEG 필터링 - 높은 p-value와 logFC로 필터링해 Signature Matrix 구성 유전자만 선별
deg_filtered <- deg_markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) %>%
  arrange(cluster, desc(avg_log2FC))


deg_filtered

write.table(deg_markers,  "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/Result/scTBI_celltype_DEG_logfc05.txt", row.names = T, col.names = F ,sep = "\t", quote = F, )


sc_TBI_only <- subset(TBI_only, TBI_Cell_Annotation %in% 
                        c("Ependymal" , "Microglia2" , "Microglia1" , "Neutrophil" , "T" , "Stromal" , "Neural_Progenitor" , "Astrocytes" , "Neuron", "B" , "Oligodendrocyte" )  )


phenotype_labels <- sc_TBI_only$TBI_Cell_Annotation

counts_subset <- as.data.frame(sc_TBI_only@assays$RNA$counts[unique(deg_filtered$gene), ])
singature_matrix <- rbind( phenotype_labels, counts_subset )
rownames(singature_matrix)  <- cbind( c("GeneSymbol", rownames(counts_subset))  )

head(singature_matrix)[1:5,1:5]

write.table(singature_matrix,  "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/Result/scTBI_neutrophi_singature_matrix.tsv", row.names = T, col.names = F ,sep = "\t", quote = F, )





# Cell chat #######################################

library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

## Whole Sample ###################################

Idents(TBI_seurat_obj) <- (TBI_seurat_obj$TBI_Cell_Annotation)

data.input <- TBI_seurat_obj[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(TBI_seurat_obj)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = TBI_seurat_obj, group.by = "ident", assay = "RNA")

# 추후에 Metadata 넣는 방법
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <-  CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# except for "Non-protein Signaling"
# CellChatDB.use <- subsetDB(CellChatDB)

# 최근 Non-protein Signaling ( metabolic and synaptic signaling)  업데이트 되었음
CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

# availableCores()  // parallel::detectCores().
# parallel::detectCores()
future::plan("multisession", workers = 32)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")




# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

saveRDS(cellchat, file = "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/Result/cellchat_scTBI.rds")


## naive Sample ###################################

TBI_seurat_obj_naive <- subset(TBI_seurat_obj ,  Condition %in% "Naïve" )

Idents(TBI_seurat_obj_naive) <- (TBI_seurat_obj_naive$TBI_Cell_Annotation)

data.input <- TBI_seurat_obj_naive[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(TBI_seurat_obj_naive)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = TBI_seurat_obj_naive, group.by = "ident", assay = "RNA")

# 추후에 Metadata 넣는 방법
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <-  CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# except for "Non-protein Signaling"
# CellChatDB.use <- subsetDB(CellChatDB)

# 최근 Non-protein Signaling ( metabolic and synaptic signaling)  업데이트 되었음
CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

# availableCores()  // parallel::detectCores().
# parallel::detectCores()
future::plan("multisession", workers = 32)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat_naive <- cellchat
saveRDS(cellchat_naive, file = "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/Result/cellchat_scTBI_Naive.rds")





## TBI Sample ###################################


TBI_seurat_obj_tbi <- subset(TBI_seurat_obj ,  Condition %in% c("CCI", "HSCCI") )

Idents(TBI_seurat_obj_tbi) <- (TBI_seurat_obj_tbi$TBI_Cell_Annotation)

data.input <- TBI_seurat_obj_tbi[["RNA"]]$data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(TBI_seurat_obj_tbi)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = TBI_seurat_obj_tbi, group.by = "ident", assay = "RNA")

# 추후에 Metadata 넣는 방법
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <-  CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

# except for "Non-protein Signaling"
# CellChatDB.use <- subsetDB(CellChatDB)

# 최근 Non-protein Signaling ( metabolic and synaptic signaling)  업데이트 되었음
CellChatDB.use <- CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

# availableCores()  // parallel::detectCores().
# parallel::detectCores()
future::plan("multisession", workers = 32)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)


cellchat_tbi <- cellchat

saveRDS(cellchat_tbi, file = "/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/Result/cellchat_scTBI_TBI.rds")














# Data analysis #######################################
library(CellChat)
library(patchwork)
library(Seurat)

setwd("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/")

TBI_seurat_obj  <- readRDS("scTBI_seurat_object.rds")
TBI_cellchat <- readRDS("cellchat_scTBI.rds")


# chatGPT에 이런식으로 물어볼 수 있습니다
# TBI_cellchat@DB$interaction 이라는 데이터프레임의 interaction_name  이라는 열에는
# 은 
# EBI3_IL27RA_IL6ST
# IL1F9_IL1RL2_IL1RAP
# SPP1_ITGAV_ITGB1
# POSTN_ITGAV_ITGB5
# SEMA3C_NRP1_PLXNA3
# LGALS9_CD44
# LAMB3_ITGA3_ITGB1
# COL6A6_ITGA9_ITGB1
# VTN_ITGAV_ITGB5
# COL1A2_CD44
# COL1A2_SDC4
# CADM1_CADM1
# CD226_NECTIN2
# 이런식으로 되어 있는데
# 
# 여기에서 
# 
# Bulk_intersected_DEGs <- c("Aspg", "C1qa", "C1qb", "C1qc", "C3", "C3ar1", "C4b", "Capg", "Cd14", "Cd22", "Cd44", 
#                            "Cd48", "Cd68", "Cd84", "Cmklr1", "Csf3r", "Ctsh", "Ctss", "Cyba", "Cybb", 
#                            "Fcgr2b", "Flnc", "Gbp2", "Gfap", "Gpnmb", "Itgb2", "Lilrb4a", "Lyz2", "Mpeg1", 
#                            "Myo1f", "Ptprc", "Slc11a1", "Tyrobp", "Unc93b1", "Vim")
# 에 해당하는 유전자들이 겹치는게 있는지 확인하고TBI_cellchat@DB$interaction  데이터프레임에서 pathway_name 열의 해당 행의 값이 뭔지 확인해줘 
# Bulk_intersected_DEGs 대소문자를 구분하지 않고 확인해주는 코드를 짜줘


Bulk_intersected_DEGs <- c("Aspg", "C1qa", "C1qb", "C1qc", "C3", "C3ar1", "C4b", "Capg", "Cd14", "Cd22", "Cd44", 
                           "Cd48", "Cd68", "Cd84", "Cmklr1", "Csf3r", "Ctsh", "Ctss", "Cyba", "Cybb", 
                           "Fcgr2b", "Flnc", "Gbp2", "Gfap", "Gpnmb", "Itgb2", "Lilrb4a", "Lyz2", "Mpeg1", 
                           "Myo1f", "Ptprc", "Slc11a1", "Tyrobp", "Unc93b1", "Vim")

TBI_cellchat@DB$interaction$interaction_name

interaction_data <- TBI_cellchat@DB$interaction

interaction_genes <- unique(unlist(strsplit(interaction_data$interaction_name, "_")))
interaction_genes_upper <- toupper(interaction_genes)

# 2. Bulk_intersected_DEGs를 대문자로 변환
Bulk_intersected_DEGs_upper <- toupper(Bulk_intersected_DEGs)

# 3. 겹치는 유전자 찾기
matching_genes <- intersect(interaction_genes_upper, Bulk_intersected_DEGs_upper)

# 4. interaction_name에서 겹치는 유전자를 포함한 행 필터링
matching_rows <- interaction_data[sapply(interaction_data$interaction_name, function(x) {
  any(toupper(unlist(strsplit(x, "_"))) %in% matching_genes)
}), ]

# 5. 결과 확인
matching_rows




library(openxlsx)

# 예제 데이터 (실제 데이터로 교체 필요)
LR_data <- TBI_cellchat@LR  # LR 데이터
DB_data <- TBI_cellchat@DB  # DB 데이터 (리스트 형태, 데이터프레임 이름 포함)

# 새 엑셀 워크북 생성
wb <- createWorkbook()

# matching_rows 데이터를 첫 번째 시트에 추가
addWorksheet(wb, "BulkDEG_Interaction")
writeData(wb, sheet = "BulkDEG_Interaction", x = matching_rows)

# LR 데이터를 두 번째 시트에 추가
addWorksheet(wb, "LR_Data")
writeData(wb, sheet = "LR_Data", x = LR_data)

# DB 데이터프레임들을 이름을 사용하여 각각의 시트에 추가
for (name in names(DB_data)) {
  addWorksheet(wb, name)  # 이름을 시트로 추가
  writeData(wb, sheet = name, x = DB_data[[name]])  # 해당 데이터 작성
}

# 엑셀 파일 저장
saveWorkbook(wb, "TBI_cellchat_Data.xlsx", overwrite = TRUE)








# TBI_cellchat@idents <- as.data.frame(TBI_cellchat@idents)
# TBI_cellchat@idents[TBI_cellchat@idents %in% "Neural_Progenitor"] <- "Neural Progenitor"
# TBI_cellchat@idents <- factor(TBI_cellchat@idents)

levels(TBI_cellchat@idents)
levels(TBI_cellchat@idents)[4] # --> Neutrophil
TBI_cellchat@idents <- factor( TBI_cellchat@idents , levels = c("Astrocytes", "B", "Ependymal", "Microglia1", "Microglia2", 
                                                                "Neural_Progenitor", "Neuron", "Neutrophil", "Oligodendrocyte", 
                                                                "Stromal", "T", "Unknown"))


netVisual_bubble(TBI_cellchat, sources.use = 2, targets.use = c(1:11), remove.isolate = FALSE)

netVisual_bubble(TBI_cellchat, sources.use = 4, targets.use = c(1:11), remove.isolate = FALSE) + scale_x_discrete( c(1:11) )


bubble_plot <- netVisual_bubble(TBI_cellchat, sources.use = "Neutrophil", targets.use = c(1:11), remove.isolate = FALSE, sort.by.target = T) 
TBI_Cellchat_bubble_plot <- bubble_plot + scale_x_discrete(labels = c("Neutrophil -> Astrocytes", 
                                                                      "Neutrophil -> B", 
                                                                      "Neutrophil -> Ependymal", 
                                                                      "Neutrophil -> Microglia1", 
                                                                      "Neutrophil -> Microglia2", 
                                                                      "Neutrophil -> Neural_Progenitor", 
                                                                      "Neutrophil -> Neuron", 
                                                                      "Neutrophil -> Neutrophil", 
                                                                      "Neutrophil -> Oligodendrocyte", 
                                                                      "Neutrophil -> Stromal", 
                                                                      "Neutrophil -> T", 
                                                                      "Neutrophil -> Unknown"))

TBI_Cellchat_bubble_plot

png("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_cellchat_Neutrophil_interaction.png", width = 4, height = 6, units = "in", res = 300)  # 300 DPI 설정
TBI_Cellchat_bubble_plot
dev.off()


bubble_plot <- netVisual_bubble(TBI_cellchat, sources.use = c(1:11), targets.use = "Neural_Progenitor" , remove.isolate = FALSE ,sort.by.source = T) 
TBI_Cellchat_bubble_plot <- bubble_plot + scale_x_discrete(labels = c("Astrocytes -> Neural Progenitor", 
                                                                      "B -> Neural Progenitor", 
                                                                      "Ependymal -> Neural Progenitor", 
                                                                      "Microglia1 -> Neural Progenitor", 
                                                                      "Microglia2 -> Neural Progenitor", 
                                                                      "Neural Progenitor -> Neural Progenitor", 
                                                                      "Neuron -> Neural Progenitor", 
                                                                      "Neutrophil -> Neural Progenitor", 
                                                                      "Oligodendrocyte -> Neural Progenitor", 
                                                                      "Stromal -> Neural Progenitor", 
                                                                      "T -> Neural Progenitor", 
                                                                      "Unknown -> Neural Progenitor") )


TBI_Cellchat_bubble_plot

png("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_cellchat_NP_interaction.png", width = 4, height = 6, units = "in", res = 300)  # 300 DPI 설정
TBI_Cellchat_bubble_plot
dev.off()






cellchat_naive
cellchat_naive@idents <- factor( cellchat_naive@idents , levels = c("Astrocytes", "B", "Ependymal", "Microglia1", "Microglia2", 
                                                                    "Neural_Progenitor", "Neuron", "Neutrophil", "Oligodendrocyte", 
                                                                    "Stromal", "T", "Unknown"))




bubble_plot <- netVisual_bubble(cellchat_naive, sources.use = "Neutrophil", targets.use = c(1:11), remove.isolate = FALSE, sort.by.target = T) 
cellchat_naive_bubble_plot <- bubble_plot + scale_x_discrete(labels = c("Neutrophil -> Astrocytes", 
                                                                        "Neutrophil -> B", 
                                                                        "Neutrophil -> Ependymal", 
                                                                        "Neutrophil -> Microglia1", 
                                                                        "Neutrophil -> Microglia2", 
                                                                        "Neutrophil -> Neural_Progenitor", 
                                                                        "Neutrophil -> Neuron", 
                                                                        "Neutrophil -> Neutrophil", 
                                                                        "Neutrophil -> Oligodendrocyte", 
                                                                        "Neutrophil -> Stromal", 
                                                                        "Neutrophil -> T", 
                                                                        "Neutrophil -> Unknown"))

cellchat_naive_bubble_plot

png("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_cellchat_naive_Neutrophil_interaction.png", width = 4, height = 6, units = "in", res = 300)  # 300 DPI 설정
cellchat_naive_bubble_plot
dev.off()


bubble_plot <- netVisual_bubble(cellchat_naive, sources.use = c(1:11), targets.use = "Neural_Progenitor" , remove.isolate = FALSE ,sort.by.source = T) 
cellchat_naive_bubble_plot <- bubble_plot + scale_x_discrete(labels = c("Astrocytes -> Neural Progenitor", 
                                                                        "B -> Neural Progenitor", 
                                                                        "Ependymal -> Neural Progenitor", 
                                                                        "Microglia1 -> Neural Progenitor", 
                                                                        "Microglia2 -> Neural Progenitor", 
                                                                        "Neural Progenitor -> Neural Progenitor", 
                                                                        "Neuron -> Neural Progenitor", 
                                                                        "Neutrophil -> Neural Progenitor", 
                                                                        "Oligodendrocyte -> Neural Progenitor", 
                                                                        "Stromal -> Neural Progenitor", 
                                                                        "T -> Neural Progenitor", 
                                                                        "Unknown -> Neural Progenitor") )


cellchat_naive_bubble_plot

png("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_cellchat_naive_NP_interaction.png", width = 4, height = 6, units = "in", res = 300)  # 300 DPI 설정
cellchat_naive_bubble_plot
dev.off()








cellchat_tbi
cellchat_tbi@idents <- factor( cellchat_tbi@idents , levels = c("Astrocytes", "B", "Ependymal", "Microglia1", "Microglia2", 
                                                                "Neural_Progenitor", "Neuron", "Neutrophil", "Oligodendrocyte", 
                                                                "Stromal", "T", "Unknown"))




bubble_plot <- netVisual_bubble(cellchat_tbi, sources.use = "Neutrophil", targets.use = c(1:11), remove.isolate = FALSE, sort.by.target = T) 
cellchat_tbi_bubble_plot <- bubble_plot + scale_x_discrete(labels = c("Neutrophil -> Astrocytes", 
                                                                      "Neutrophil -> B", 
                                                                      "Neutrophil -> Ependymal", 
                                                                      "Neutrophil -> Microglia1", 
                                                                      "Neutrophil -> Microglia2", 
                                                                      "Neutrophil -> Neural_Progenitor", 
                                                                      "Neutrophil -> Neuron", 
                                                                      "Neutrophil -> Neutrophil", 
                                                                      "Neutrophil -> Oligodendrocyte", 
                                                                      "Neutrophil -> Stromal", 
                                                                      "Neutrophil -> T", 
                                                                      "Neutrophil -> Unknown"))

cellchat_tbi_bubble_plot

png("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_cellchat_tbi_Neutrophil_interaction.png", width = 4, height = 6, units = "in", res = 300)  # 300 DPI 설정
cellchat_tbi_bubble_plot
dev.off()


bubble_plot <- netVisual_bubble(cellchat_tbi, sources.use = c(1:11), targets.use = "Neural_Progenitor" , remove.isolate = FALSE ,sort.by.source = T) 
cellchat_tbi_bubble_plot <- bubble_plot + scale_x_discrete(labels = c("Astrocytes -> Neural Progenitor", 
                                                                      "B -> Neural Progenitor", 
                                                                      "Ependymal -> Neural Progenitor", 
                                                                      "Microglia1 -> Neural Progenitor", 
                                                                      "Microglia2 -> Neural Progenitor", 
                                                                      "Neural Progenitor -> Neural Progenitor", 
                                                                      "Neuron -> Neural Progenitor", 
                                                                      "Neutrophil -> Neural Progenitor", 
                                                                      "Oligodendrocyte -> Neural Progenitor", 
                                                                      "Stromal -> Neural Progenitor", 
                                                                      "T -> Neural Progenitor", 
                                                                      "Unknown -> Neural Progenitor") )


cellchat_tbi_bubble_plot

png("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/scTBI_cellchat_tbi_NP_interaction.png", width = 4, height = 6, units = "in", res = 300)  # 300 DPI 설정
cellchat_tbi_bubble_plot
dev.off()









# scTBI analysis #########################################
library(patchwork)
library(Seurat)
library(reshape2)

setwd("/home/Data_Drive_8TB_3/ch2860/RNAcourse/scData/CellNeuron_data/result/")

TBI_seurat_obj  <- readRDS("scTBI_seurat_object.rds")

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

Idents(TBI_seurat_obj) <- factor(TBI_seurat_obj$TBI_Cell_Annotation, levels = c("Astrocytes", "B", "Ependymal", "Microglia1", "Microglia2", 
                                                                                "Neural Progenitor", "Neuron", "Neutrophil", "Oligodendrocyte", 
                                                                                "Stromal", "T", "Unknown"))

TBI_umap <- DimPlot(TBI_seurat_obj ,label=T , label.size = 6, raster = F,  reduction = "umap.harmony" ,  repel = T,
                    cols = cell_colors )+ ggtitle(NULL) + xlim(-16, 15) + ylim(-18, 12)  + NoAxes()  +theme(legend.text = element_text(margin = margin(t = 5, b = 5)) )



TBI_umap <- TBI_umap +
  # X축 화살표
  annotate("segment", x = -15, xend = -10, y = -17, yend = -17, 
           arrow = arrow(length = unit(0.2, "cm")), size = 0.8) +  # 화살표
  annotate("text", x = -10, y = -17.5, label = "umap1", size = 5, vjust = 1, hjust=1) +  # X축 레이블
  
  # Y축 화살표
  annotate("segment", x = -15, xend = -15, y = -17, yend = -13, 
           arrow = arrow(length = unit(0.2, "cm")), size = 0.8) +  # 화살표
  annotate("text", x = -15.5, y = -13, label = "umap2", size = 5, hjust = 1, vjust=0 , angle= 90)     # Y축 레이블

TBI_umap


ggsave(
  filename = "DimPlot_TBI_harmony_UMAP.png", 
  plot = TBI_umap, 
  width = 10,        # 가로 크기 (인치)
  height = 8,        # 세로 크기 (인치)
  dpi = 300          # 해상도 설정 (300 DPI)
)







TBI_seurat_obj$TBI_Cell_Annotation
TBI_seurat_obj$Condition
table(TBI_seurat_obj$Replicate)


# TBI_seurat_obj$Condition ###########################
prop_table_tmp <- prop.table(table(TBI_seurat_obj$TBI_Cell_Annotation, TBI_seurat_obj$Condition, TBI_seurat_obj$Sample_ID),margin = 3)

means_df_condition <- data.frame(
  Naïve = rowMeans( prop_table_tmp[, "Naïve", names(table(TBI_seurat_obj$Sample_ID[TBI_seurat_obj$Condition == "Naïve"])) ] ) ,
  CCI = rowMeans( prop_table_tmp[, "CCI", names(table(TBI_seurat_obj$Sample_ID[TBI_seurat_obj$Condition == "CCI"])) ] ),
  HSCCI = rowMeans( prop_table_tmp[, "HSCCI", names(table(TBI_seurat_obj$Sample_ID[TBI_seurat_obj$Condition == "HSCCI"])) ] )
)

means_df_condition


means_df_condition$Cell_Type <- rownames(means_df_condition)  # 세포 타입 열 추가
means_df_melted <- melt(means_df_condition, id.vars = "Cell_Type", variable.name = "Condition", value.name = "Proportion")

# Stacked barplot 생성
stacked_barplot_condition <- ggplot(means_df_melted, aes(x = Condition, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors) +  # cell_colors를 사용한 색상 매핑
  labs(title = "Proportions of Cell Types by Condition",
       x = "Condition",
       y = "Proportion",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))

ggsave(
  filename = "StackedBar_Cellpopulation_Condition.png", 
  plot = stacked_barplot_condition, 
  width = 4.5,        # 가로 크기 (인치)
  height = 4,        # 세로 크기 (인치)
  dpi = 300          # 해상도 설정 (300 DPI)
)



# TBI_seurat_obj$Timepoint ###########################

prop_table_tmp <- prop.table(table(TBI_seurat_obj$TBI_Cell_Annotation, TBI_seurat_obj$Timepoint, TBI_seurat_obj$Sample_ID),margin = 3)

means_df_timepoint <- data.frame(
  Timepoint_24h = rowMeans(prop_table_tmp[, "24 hours", 
                                          names(table(TBI_seurat_obj$Sample_ID[TBI_seurat_obj$Timepoint == "24 hours"]))]),
  Timepoint_6m = rowMeans(prop_table_tmp[, "6 month", 
                                         names(table(TBI_seurat_obj$Sample_ID[TBI_seurat_obj$Timepoint == "6 month"]))]),
  Timepoint_7d = rowMeans(prop_table_tmp[, "7 day", 
                                         names(table(TBI_seurat_obj$Sample_ID[TBI_seurat_obj$Timepoint == "7 day"]))])
)

# 세포 타입 추가
means_df_timepoint$Cell_Type <- rownames(means_df_timepoint)

# 데이터 변환: 긴 형식으로 변환
means_df_melted <- melt(means_df_timepoint, id.vars = "Cell_Type", 
                        variable.name = "Timepoint", value.name = "Proportion")

# 수동으로 이름 변경
means_df_melted$Timepoint <- factor(means_df_melted$Timepoint, 
                                    levels = c("Timepoint_24h", "Timepoint_7d",  "Timepoint_6m"),
                                    labels = c("24 hours", "7 days" ,"6 months"))

# ggplot 생성
stacked_barplot_timepoint <- ggplot(means_df_melted, aes(x = Timepoint, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors) +  # 세포 타입별 색상 지정
  labs(title = "Proportions of Cell Types by Timepoint",
       x = "Timepoint",
       y = "Proportion",
       fill = "Cell Type") +
  theme_minimal()

# 플롯 출력
print(stacked_barplot_timepoint)



ggsave(
  filename = "StackedBar_Cellpopulation_timepoint.png", 
  plot = stacked_barplot_timepoint, 
  width = 4.5,        # 가로 크기 (인치)
  height = 4,        # 세로 크기 (인치)
  dpi = 300          # 해상도 설정 (300 DPI)
)

# 2883                4149                4602               10377                8089                4342 


# TBI_seurat_obj$Replicate ###########################
names(table(TBI_seurat_obj$Replicate) ) ->  replicates

CCI_24h <- replicates[grep("24h_male_CCI", replicates)]
HSCCI_24h <- replicates[grep("24h_male_HSCCI", replicates)]
Naïve_24h <- replicates[grep("24h_male_Naïve", replicates)]
CCI_6m <- replicates[grep("6m_male_CCI", replicates)]
CCI_7d <- replicates[grep("7d_male_CCI", replicates)]


prop_table_tmp <- prop.table(table(TBI_seurat_obj$TBI_Cell_Annotation, TBI_seurat_obj$Replicate, TBI_seurat_obj$Sample_ID),margin = 3)

prop_table_tmp <- prop.table(table(TBI_seurat_obj$TBI_Cell_Annotation, TBI_seurat_obj$Replicate),margin = 2)


means_df_group <- data.frame(
  Naïve_24h = rowMeans(prop_table_tmp[, Naïve_24h  ]  ),
  CCI_24h = rowMeans(prop_table_tmp[, CCI_24h]),
  HSCCI_24h = rowMeans(prop_table_tmp[, HSCCI_24h]),
  
  CCI_7d = rowMeans(prop_table_tmp[, CCI_7d]),
  CCI_6m = rowMeans(prop_table_tmp[,CCI_6m ])
)


# 세포 타입 추가
means_df_group$Cell_Type <- rownames(means_df_group)

# 데이터 변환: 긴 형식으로 변환
means_df_melted <- melt(means_df_group, id.vars = "Cell_Type", 
                        variable.name = "group", value.name = "Proportion")

# 수동으로 이름 변경
means_df_melted$group <- factor(means_df_melted$group, 
                                levels = c("Naïve_24h", "CCI_24h", "HSCCI_24h",    "CCI_7d", "CCI_6m"),
                                labels = c("Naïve\n24 hours", "CCI\n24 hours" ,"HSCCI\n24 hours","CCI\n7 days","CCI\n6 months"))

# ggplot 생성
stacked_barplot_group<- ggplot(means_df_melted, aes(x = group, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors) +  # 세포 타입별 색상 지정
  labs(title = "Proportions of Cell Types by group",
       x = "Group",
       y = "Proportion",
       fill = "Cell Type") +
  theme_minimal()

# 플롯 출력
print(stacked_barplot_group)



ggsave(
  filename = "StackedBar_Cellpopulation_group.png", 
  plot = stacked_barplot_group, 
  width = 5,        # 가로 크기 (인치) 
  height = 4,        # 세로 크기 (인치) 
  dpi = 300          # 해상도 설정 (300 DPI)
)






