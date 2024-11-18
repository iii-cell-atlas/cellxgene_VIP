## script as part of paraCell Toxoplasma-mouse atlas. 

## See paper for the accession number.

library(Seurat)

#Run1---------------------------------------------------------------------------------------------------------------------------------------------

m1.data <- Read10X(data.dir = "Mouse_Toxo/BMDMplusRH-IFNY_UTR/")

run1 <- CreateSeuratObject(counts = m1.data, project = "run1", min.cells = 3, min.features = 200)

geneNames = rownames(run1)

head(geneNames) 
#[1] "TgondiiGT1-2.5UTR-TGGT1-308090" "TgondiiGT1-2.5UTR-TGGT1-408710" "TgondiiGT1-2.5UTR-TGGT1-360830" "TgondiiGT1-2.5UTR-TGGT1-409590"
#[5] "TgondiiGT1-2.5UTR-TGGT1-409600" "TgondiiGT1-2.5UTR-TGGT1-363030"

tail(geneNames)
#[1] "Mmus--------------Csprs"          "Mmus--------------AC125149.3"     "Mmus--------------AC168977.1"     "Mmus--------------AC149090.1"    
#[5] "Mmus--------------CAAA01118383.1" "Mmus--------------CAAA01147332.1"

newGeneNames = vector()

for (x in geneNames){
  if (startsWith(x,"Mmus")){ #if host
    genelist = as.list(strsplit(x, "--------------")[[1]])
    newname = genelist[2]
    newname = paste0("H-",newname)
    newGeneNames = c(newGeneNames,newname)
  }else{ #if parasite
    genelist = as.list(strsplit(x, ".5UTR-")[[1]])
    newname = genelist[2]
    newname = paste0("P-",newname)
    newGeneNames = c(newGeneNames,newname)
  }
}

newGeneNames = as.character(newGeneNames)

head(newGeneNames)
#"P-TGGT1-308090" "P-TGGT1-408710" "P-TGGT1-360830" "P-TGGT1-409590" "P-TGGT1-409600" "P-TGGT1-363030"

tail(newGeneNames)
#"H-Csprs"  "H-AC125149.3" "H-AC168977.1" "H-AC149090.1"  "H-CAAA01118383.1" "H-CAAA01147332.1"

counts <- GetAssayData(run1,assay = "RNA",slot = "counts")

rownames(counts) = newGeneNames

run1_new = CreateSeuratObject(counts = counts, meta.data = run1@meta.data)

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

run1_new <- subset(run1_new, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)

run1_new <- NormalizeData(run1_new)

run1_new <- FindVariableFeatures(run1_new, selection.method = "vst", nfeatures = 2000)

head(run1_new@assays$RNA@var.features) # "H-S100a9" "H-Ngp"    "H-S100a8" "H-Ltf"    "H-Lcn2"   "H-Chil3" 
tail(run1_new@assays$RNA@var.features) # "H-Sirpb1a" "P-TGGT1-285250" "P-TGGT1-246550" "H-Timeless"  "P-TGGT1-259260" "P-TGGT1-229360"

#Run2-------------------------------------------------------------------------------------------------------------------------------------

m2.data <- Read10X(data.dir = "Mouse_Toxo/BMDMplusRHplusIFNY_UTR/")

run2 <- CreateSeuratObject(counts = m2.data, project = "run2", min.cells = 3, min.features = 200)

geneNames = rownames(run2)

head(geneNames) 
#[1] "TgondiiGT1-2.5UTR-TGGT1-408710" "TgondiiGT1-2.5UTR-TGGT1-360830" "TgondiiGT1-2.5UTR-TGGT1-409590" "TgondiiGT1-2.5UTR-TGGT1-409600"
#[5] "TgondiiGT1-2.5UTR-TGGT1-363030" "TgondiiGT1-2.5UTR-TGGT1-409990"

tail(geneNames)
#[1] "Mmus--------------AC125149.3"     "Mmus--------------AC125149.2"     "Mmus--------------AC168977.1"     "Mmus--------------AC149090.1"    
#[5] "Mmus--------------CAAA01118383.1" "Mmus--------------CAAA01147332.1"

newGeneNames = vector()

for (x in geneNames){
  if (startsWith(x,"Mmus")){ #if host
    genelist = as.list(strsplit(x, "--------------")[[1]])
    newname = genelist[2]
    newname = paste0("H-",newname)
    newGeneNames = c(newGeneNames,newname)
  }else{ #if parasite
    genelist = as.list(strsplit(x, ".5UTR-")[[1]])
    newname = genelist[2]
    newname = paste0("P-",newname)
    newGeneNames = c(newGeneNames,newname)
  }
}

newGeneNames = as.character(newGeneNames)

head(newGeneNames)
# "P-TGGT1-408710" "P-TGGT1-360830" "P-TGGT1-409590" "P-TGGT1-409600" "P-TGGT1-363030" "P-TGGT1-409990"

tail(newGeneNames)
#"H-AC125149.3" "H-AC125149.2" "H-AC168977.1" "H-AC149090.1" "H-CAAA01118383.1" "H-CAAA01147332.1"

counts <- GetAssayData(run2,assay = "RNA",slot = "counts")

rownames(counts) = newGeneNames

run2_new = CreateSeuratObject(counts = counts, meta.data = run2@meta.data)

#add host-parasite metadata

run2_new[["percent.parasite"]] <- PercentageFeatureSet(run2_new, pattern = "P-")
run2_new[["parasiteUMI"]] <- PercentageFeatureSet(run2_new, pattern = "P-") * run2_new$nCount_RNA/100
run2_new[["percent.host"]] <- PercentageFeatureSet(run2_new, pattern = "H-")
run2_new[["hostUMI"]] <- PercentageFeatureSet(run2_new, pattern = "H-") * run2_new$nCount_RNA/100

#quality control

VlnPlot(run2_new, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) 

VlnPlot(run2_new, features = c("percent.parasite", "parasiteUMI"), ncol = 2) 

VlnPlot(run2_new, features = c("percent.host", "hostUMI"), ncol = 2) 

run2_new <- subset(run2_new, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)

#normalize and find variable featuers

run2_new <- NormalizeData(run2_new)

run2_new <- FindVariableFeatures(run2_new, selection.method = "vst", nfeatures = 2000)

head(run2_new@assays$RNA@var.features) # "H-S100a9" "H-Ngp"    "H-S100a8" "H-Ltf"    "H-Lcn2"   "H-Chil3" 
tail(run2_new@assays$RNA@var.features) # "H-Sirpb1a" "P-TGGT1-285250" "P-TGGT1-246550" "H-Timeless"  "P-TGGT1-259260" "P-TGGT1-229360"

# Integration ---------------------------------------------------------------------------------------------------------------------------------

object_list = c(run1_new,run2_new)

features <- SelectIntegrationFeatures(object.list = object_list)

mm.anchors <- FindIntegrationAnchors(object.list = object_list, anchor.features = features)

mm.combined <- IntegrateData(anchorset = mm.anchors)

DefaultAssay(mm.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
mm.combined <- ScaleData(mm.combined, verbose = FALSE)
mm.combined <- RunPCA(mm.combined, npcs = 50, verbose = FALSE)

ElbowPlot(mm.combined, ndims = 50)

mm.combined <- RunUMAP(mm.combined, reduction = "pca", dims = 1:30)
mm.combined <- FindNeighbors(mm.combined, reduction = "pca", dims = 1:30)
mm.combined <- FindClusters(mm.combined, resolution = 0.5)

p1 <- DimPlot(mm.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(mm.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#save RDS object 

saveRDS(mm.combined, "Mouse_Toxo/mouse_toxo.combined.rds")


