## script as part of paraCell Theileria-cow atlas. 

## See paper for the accession number.

library(Seurat)
library(stringr)

# Run 1 ----------------------------------------------------------------------------------------------------------------------------------------

run1.data <- Read10X(data.dir = "Theileria_Data/82H_UTR/")

run1 <- CreateSeuratObject(counts = run1.data, project = "run1", min.cells = 3, min.features = 200)

#reformat gene names

geneNames = rownames(run1)

head(geneNames)
#[1] "Theileria-UTR-Tap370b08.q2ca38.01"  "Theileria-UTR-Tap370b08.q2ca38.02c" "Theileria-UTR-Tap370b08.q2ca38.03c"
#[4] "Theileria-UTR-TA16055"  

tail(geneNames)
#"Cow-----------ARSD" "Cow-----------GYG2" "Cow-----------XG" "Cow-----------CD99"              
#[5] "Cow-----------ENSBTAG00000049411" "Cow-----------ENSBTAG00000044391"

newGeneNames = vector()

for (x in geneNames){
  if (startsWith(x,"Cow")){ #if host
    genelist = as.list(strsplit(x, "-----------")[[1]])
    newname = genelist[2]
    newname = paste0("H-",newname)
    newGeneNames = c(newGeneNames,newname)
  }else{ #if parasite
    genelist = as.list(strsplit(x, "-UTR-")[[1]])
    newname = genelist[2]
    newname = paste0("P-",newname)
    newGeneNames = c(newGeneNames,newname)
  }
}

newGeneNames = as.character(newGeneNames)

head(newGeneNames)
#[1] "P-Tap370b08.q2ca38.01"  "P-Tap370b08.q2ca38.02c" "P-Tap370b08.q2ca38.03c" "P-TA16055" "P-TA16050"             
#[6] "P-TA16045" 

tail(newGeneNames)
#[1] "H-ARSD" "H-GYG2" "H-XG" "H-CD99" "H-ENSBTAG00000049411" "H-ENSBTAG00000044391"

counts <- GetAssayData(run1,assay = "RNA",slot = "counts")
rownames(counts) = newGeneNames

run1_new = CreateSeuratObject(counts = counts, meta.data = run1@meta.data)

#saveRDS(run1_new, "run1_new.rds")

#rm(run1)

#add host-parasite metadata

run1_new[["percent.parasite"]] <- PercentageFeatureSet(run1_new, pattern = "P-")
run1_new[["parasiteUMI"]] <- PercentageFeatureSet(run1_new, pattern = "P-") * run1_new$nCount_RNA/100
run1_new[["percent.host"]] <- PercentageFeatureSet(run1_new, pattern = "H-")
run1_new[["hostUMI"]] <- PercentageFeatureSet(run1_new, pattern = "H-") * run1_new$nCount_RNA/100

FeatureScatter(run1_new, feature1 = "percent.parasite", feature2 = "percent.host") #perfect correlation ,as expected

#quality control

VlnPlot(run1_new, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) 

VlnPlot(run1_new, features = c("percent.parasite", "parasiteUMI"), ncol = 2) 

VlnPlot(run1_new, features = c("percent.host", "hostUMI"), ncol = 2) 

FeatureScatter(run1_new, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

run1_new <- subset(run1_new, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500)

run1_new <- subset(run1_new, subset = percent.parasite < 50 & parasiteUMI < 5000)

run1_new <- subset(run1_new, subset = percent.host > 75)


run1_new <- NormalizeData(run1_new)

run1_new <- FindVariableFeatures(run1_new, selection.method = "vst", nfeatures = 2000)

head(run1_new@assays$RNA@var.features) # "H-CCL5"    "H-TNFSF15" "H-IFIT2"   "H-FLT1"    "H-CXCL2"   "H-VIM"
tail(run1_new@assays$RNA@var.features) # "H-TET3" "H-ENSBTAG00000004415" "H-ENSBTAG00000043570" "H-DGKK"  "H-SLC12A6" "P-TA14525"
#plenty of variation in host genes - a good sign.

# Run 2 ----------------------------------------------------------------------------------------------------------------------------------------

run2.data <- Read10X(data.dir = "Theileria_Data/12886_UTR/filtered_feature_bc_matrix")

run2 <- CreateSeuratObject(counts = run2.data, project = "run2", min.cells = 3, min.features = 200)

geneNames = rownames(run2)

head(geneNames)
#[1] "Theileria-UTR-Tap370b08.q2ca38.01"  "Theileria-UTR-Tap370b08.q2ca38.02c" "Theileria-UTR-Tap370b08.q2ca38.03c"
#[4] "Theileria-UTR-TA16055"              "Theileria-UTR-TA16050"              "Theileria-UTR-TA16045" 

tail(geneNames)
#"Cow-----------ARSD"               "Cow-----------GYG2"               "Cow-----------XG"                 "Cow-----------CD99"              
#[5] "Cow-----------ENSBTAG00000049411" "Cow-----------ENSBTAG00000044391"

newGeneNames = vector()

for (x in geneNames){
  if (startsWith(x,"Cow")){ #if host
    genelist = as.list(strsplit(x, "-----------")[[1]])
    newname = genelist[2]
    newname = paste0("H-",newname)
    newGeneNames = c(newGeneNames,newname)
  }else{ #if parasite
    genelist = as.list(strsplit(x, "-UTR-")[[1]])
    newname = genelist[2]
    newname = paste0("P-",newname)
    newGeneNames = c(newGeneNames,newname)
  }
}

newGeneNames = as.character(newGeneNames)

head(newGeneNames)
#[1] "P-Tap370b08.q2ca38.01"  "P-Tap370b08.q2ca38.02c" "P-Tap370b08.q2ca38.03c" "P-TA16055"              "P-TA16050"             
#[6] "P-TA16045" 

tail(newGeneNames)
#[1] "H-ARSD" "H-GYG2" "H-XG" "H-CD99" "H-ENSBTAG00000049411" "H-ENSBTAG00000044391"

counts <- GetAssayData(run2,assay = "RNA",slot = "counts")

rownames(counts) = newGeneNames

run2_new = CreateSeuratObject(counts = counts, meta.data = run2@meta.data)

#saveRDS(run2_new, "run2_new.rds")

#add host-parasite metadata

run2_new[["percent.parasite"]] <- PercentageFeatureSet(run2_new, pattern = "P-")
run2_new[["parasiteUMI"]] <- PercentageFeatureSet(run2_new, pattern = "P-") * run2_new$nCount_RNA/100
run2_new[["percent.host"]] <- PercentageFeatureSet(run2_new, pattern = "H-")
run2_new[["hostUMI"]] <- PercentageFeatureSet(run2_new, pattern = "H-") * run2_new$nCount_RNA/100

FeatureScatter(run2_new, feature1 = "percent.parasite", feature2 = "percent.host") 

#quality control

VlnPlot(run2_new, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) 

VlnPlot(run2_new, features = c("percent.parasite", "parasiteUMI"), ncol = 2) 

VlnPlot(run2_new, features = c("percent.host", "hostUMI"), ncol = 2) 

FeatureScatter(run2_new, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

run2_new <- subset(run2_new, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500)

run2_new <- subset(run2_new, subset = percent.parasite < 25 & parasiteUMI < 2500)

run2_new <- subset(run2_new, subset = percent.host > 75)

saveRDS(run2_new, "run2_new.rds")

#Normalize and Find Variable Features

run2_new <- NormalizeData(run2_new)

run2_new <- FindVariableFeatures(run2_new, selection.method = "vst", nfeatures = 2000)

head(run2_new@assays$RNA@var.features) # "H-HS6ST3"  "H-HMOX1"   "H-CD52"    "H-CCL5"    "H-VIM"     "H-TNFSF15"
tail(run2_new@assays$RNA@var.features) #  "H-BATF"    "H-SLC49A4" "P-TA07590" "H-NUP210L" "H-NCAM2"   "H-TCF7L1"

#Integrate Data -------------------------------------------------------------------------------------------------------------------------------

object_list = c(run1_new,run2_new)

features <- SelectIntegrationFeatures(object.list = object_list)

hp.anchors <- FindIntegrationAnchors(object.list = object_list, anchor.features = features)

hp.combined <- IntegrateData(anchorset = hp.anchors)

DefaultAssay(hp.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
hp.combined <- ScaleData(hp.combined, verbose = FALSE)
hp.combined <- RunPCA(hp.combined, npcs = 50, verbose = FALSE)

ElbowPlot(hp.combined, ndims = 50)

hp.combined <- RunUMAP(hp.combined, reduction = "pca", dims = 1:30)
hp.combined <- FindNeighbors(hp.combined, reduction = "pca", dims = 1:30)
hp.combined <- FindClusters(hp.combined, resolution = 0.5)

p1 <- DimPlot(hp.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(hp.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


DimPlot(hp.combined, reduction = "umap", split.by = "orig.ident")

#saveRDS(hp.combined, "Theileria_Data/theileria.combined.rds") - without filtering by host parasite metadate

saveRDS(hp.combined, "Theileria_Data/theileria.combined2.rds") # - with filtering by host parasite metadata





