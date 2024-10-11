library(Seurat)
library(ggplot2)


# Save & Load Main object
# saveRDS(Main_obj, 'GSE176171_Main_obj.rds')
Main_obj <- readRDS('scRNAseq/UMAP_1st/GSE176171_Main_obj.rds')


# Load raw data
list.files("/Data_1/Raw_Data/adipocyte/Other_papers/Emont_etal.2022.Nature/GSE176171_Mm10X/")
Main_obj <- Read10X(data.dir = "/Data_1/Raw_Data/adipocyte/Other_papers/Emont_etal.2022.Nature/GSE176171_Mm10X/")
Main_obj <- CreateSeuratObject(counts = Main_obj, project = "mouseWAT", min.cells = 3, min.features = 200)

# Basic processing -------------------------------------------------------------
Main_obj[["percent.mt"]] <- PercentageFeatureSet(Main_obj, pattern = "^mt-")

## Violin plots
# plot <- VlnPlot(Main_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=TRUE, pt.size = 0)
# options(repr.plot.width = 6, repr.plot.height = 6)
# pdf(file="Rplots.pdf")
# plot
# dev.off()

## Scatter plots
# plot1 <- FeatureScatter(Main_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Main_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# options(repr.plot.width = 12, repr.plot.height = 6)
# plot1 + plot2

Main_obj <- subset(Main_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 2)

Main_obj <- NormalizeData(Main_obj, normalization.method = "LogNormalize", scale.factor = 10000)

Main_obj <- FindVariableFeatures(Main_obj, selection.method = "vst", nfeatures = 4000)
# plot1 <- VariableFeaturePlot(Main_obj)
# pdf(file="VariableFeaturePlot.pdf")
# plot1
# dev.off()

all.genes <- rownames(Main_obj)
Main_obj <- ScaleData(Main_obj, features = all.genes)


# PCA & UMAP -------------------------------------------------------------------
# 1st
Main_obj <- RunPCA(Main_obj, features = VariableFeatures(object = Main_obj))

pdf(file="ElbowPlot.pdf")
ElbowPlot(Main_obj, ndims = 100)
dev.off()

Slct_dims <- c(1:9)
Main_obj <- FindNeighbors(Main_obj, dims = Slct_dims)
Main_obj <- FindClusters(Main_obj, resolution = 0.5)
Main_obj$cluster <- Idents(Main_obj)
#
Main_obj <- RunUMAP(Main_obj, dims = Slct_dims,
                 n.neighbors = 30, #default 30
                 min.dist = 0.3) #default 0.3


# UMAP plots -------------------------------------------------------------------
Idents(Main_obj) <- Main_obj$cluster
UMAP <- DimPlot(Main_obj, reduction = "umap", pt.size = 1)
UMAP_label <- LabelClusters(plot = UMAP, id = 'ident')

pdf(file="UMAP.pdf")
UMAP
dev.off()
pdf(file="UMAP_label.pdf")
UMAP_label
dev.off()


# cell type annotation -----------------------------------------------------------
pdf(file="Major_markers_on_paper.pdf", width = 16, height = 20)
FeaturePlot(object = Main_obj,
            features = c('Adipoq','Pdgfra','Msln','Jam2','Prox1','Steap4','Myocd','Mafb','Cybb','Flt3','Cpa3','Csf3r','Ms4a1','Klrd1','Il7r','Dcdc2a','Erbb4'),
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

pdf(file="Major_markers_on_other_paper.pdf", width = 12, height = 30) # https://doi.org/10.1016/j.celrep.2022.111046
FeaturePlot(object = Main_obj,
            features = c('Adipoq','Lipe','Plin1',
                        'Fkbp10','Col6a1','Col6a2',
                        'Upk3b','Msln','Krt19',
                        'Mmrn2','Esam','Cdh5',
                        'Kcnmb1','Cnn1','Myh11',
                        'Cd68','C1qc','Fcer1g',
                        'Csf3r','Fcgr3','Cxcr2',
                        'Cpa3','Tpsb2','Tpsab1',
                        'Trbc2','Cd6','Cd3e',
                        'Igkc','Jchain','Mzb1'),
            reduction = "umap", min.cutoff = 0,
            ncol=3,
            pt.size = 1)
dev.off()

pdf(file="Adipocyte_markers.pdf", width = 8, height = 12)
FeaturePlot(object = Main_obj,
            features = c('Adipoq','Pnpla2','Plin1','Lipe','Pparg','Fabp4'),
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

pdf(file="Adipocyte_subclustering_markers.pdf", width = 8, height = 12)
FeaturePlot(object = Main_obj,
            features = c('Ces1f','Btc','Apoe','Cacna1a','Prune2','Mt2'),
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

pdf(file="Thermogenic_Ad_markers.pdf", width = 8, height = 8)
FeaturePlot(object = Main_obj,
            features = c('Ucp1','Cidea','Hspb7','Zic1'),
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

pdf(file="ASPC_markers.pdf", width = 8, height = 4)
FeaturePlot(object = Main_obj,
            features = c('Pdgfra','Cspg4'),
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

# cell type labeling -----------------------------------------------------------
Main_obj$cell_type_1 <- Main_obj$cluster
Idents(Main_obj) <- Main_obj$cell_type_1
Main_obj <- RenameIdents(Main_obj,
                        '0'='ASPC', '1'='ASPC', '13'='ASPC', '22'='ASPC',
                        '2'='Adipocyte', '7'='Adipocyte', '10'='Adipocyte', '16'='Adipocyte',
                        '3'='Immune cells', '4'='Immune cells', '5'='Immune cells', '8'='Immune cells', '9'='Immune cells', '14'='Immune cells',
                        '6'='Mesothelial', 
                        '11'='Endothelial',
                        '12'='Epithelial', '15'='Epithelial', '19'='Others', '20'='Epithelial', '21'='Epithelial',
                        '17'='Smooth muscle', '18'='Smooth muscle')
Main_obj$cell_type_1 <- Idents(Main_obj)

Idents(Main_obj) <- Main_obj$cell_type_1
UMAP_cell_type1 <- DimPlot(Main_obj, reduction = "umap",pt.size = 1, shuffle = T)
UMAP_cell_type1 <- LabelClusters(plot = UMAP_cell_type1, id = 'ident')
pdf(file="UMAP_cell_type1.pdf")
UMAP_cell_type1
dev.off()

Main_obj$cell_type_2 <- Main_obj$cluster
Idents(Main_obj) <- Main_obj$cell_type_2
Main_obj <- RenameIdents(Main_obj,
                        '0'='ASPC', '1'='ASPC', '13'='ASPC', '22'='ASPC',
                        '2'='Adipocyte', '7'='Adipocyte', '10'='Adipocyte', '16'='Adipocyte',
                        '3'='Immune cells', '4'='Immune cells', '5'='Immune cells', '8'='Immune cells', '9'='Immune cells', '14'='Immune cells',
                        '6'='Mesothelial', 
                        '11'='Endothelial',
                        '12'='Epithelial', '15'='Epithelial', '19'='Epithelial', '20'='Epithelial', '21'='Epithelial',
                        '17'='Smooth muscle', '18'='Smooth muscle')
Main_obj$cell_type_2 <- Idents(Main_obj)

Idents(Main_obj) <- Main_obj$cell_type_2
UMAP_cell_type1 <- DimPlot(Main_obj, reduction = "umap",pt.size = 1, shuffle = T)
UMAP_cell_type1 <- LabelClusters(plot = UMAP_cell_type1, id = 'ident')
pdf(file="UMAP_cell_type2.pdf")
UMAP_cell_type1
dev.off()

# Gls Glul check -----------------------------------------------------------
pdf(file="Gls_Glul_scExp.pdf", width = 8, height = 8)
FeaturePlot(object = Main_obj,
            features = c('Gls','Glul','Glud1','Gpt2'),
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

pdf(file="Gls_Glul_violin.pdf", width = 12, height = 8)
VlnPlot(Main_obj, features = c('Gls','Glul','Glud1','Gpt2'),
        pt.size = 0)
dev.off()

# Metadata -----------------------------------------------------------
tmp_data <- data.frame(cell_ID = Main_obj@assays[["RNA"]]@data@Dimnames[[2]],
                       cell_type = as.character(Main_obj@meta.data[["cell_type_2"]]),
                       Adipoq = as.numeric(Main_obj@assays[['RNA']]@counts['Adipoq',]),
                       Lipe = as.numeric(Main_obj@assays[['RNA']]@counts['Lipe',]),
                       Plin1 = as.numeric(Main_obj@assays[['RNA']]@counts['Plin1',]),
                       Pdgfra = as.numeric(Main_obj@assays[['RNA']]@counts['Pdgfra',]),
                       Col6a1 = as.numeric(Main_obj@assays[['RNA']]@counts['Col6a1',]),
                       Col6a2 = as.numeric(Main_obj@assays[['RNA']]@counts['Col6a2',]),
                       Gls = as.numeric(Main_obj@assays[['RNA']]@counts['Gls',]),
                       Glul = as.numeric(Main_obj@assays[['RNA']]@counts['Glul',])
                       )
tmp_data <- cbind(Main_obj@reductions[["umap"]]@cell.embeddings, tmp_data)

write.table(tmp_data, 'Emont_UMAP_metadata1.tsv', sep = '\t', col.names = T, row.names = F)

# scExp data -----------------------------------------------------------
tmp <- data.frame(t(Main_obj@assays[["RNA"]]@scale.data))

tmp_data <- data.frame(cell_ID = Main_obj@assays[["RNA"]]@data@Dimnames[[2]],
                       cell_type = as.character(Main_obj@meta.data[["cell_type_2"]])
                       )
tmp_data <- cbind(Main_obj@reductions[["umap"]]@cell.embeddings, tmp_data)
tmp_data <- cbind(tmp_data, tmp)

write.table(tmp_data, 'scRNAseq/UMAP_1st/scExp_data.tsv', sep='\t', row.names = F, col.names = T)

# scExp data / selected 1 -----------------------------------------------------------
tmp_data <- data.frame(cell_ID = Main_obj@assays[["RNA"]]@data@Dimnames[[2]],
                       cell_type = as.character(Main_obj@meta.data[["cell_type_2"]]),
                       Adipoq = as.numeric(Main_obj@assays[['RNA']]@scale.data['Adipoq',]),
                       Pnpla2 = as.numeric(Main_obj@assays[['RNA']]@scale.data['Pnpla2',]),
                       Plin1 = as.numeric(Main_obj@assays[['RNA']]@scale.data['Plin1',]),
                       Lipe = as.numeric(Main_obj@assays[['RNA']]@scale.data['Lipe',]),
                       Pdgfra = as.numeric(Main_obj@assays[['RNA']]@scale.data['Pdgfra',]),
                       Gls = as.numeric(Main_obj@assays[['RNA']]@scale.data['Gls',]),
                       Glul = as.numeric(Main_obj@assays[['RNA']]@scale.data['Glul',]),
                       Glud1 = as.numeric(Main_obj@assays[['RNA']]@scale.data['Glud1',]),
                       Gpt2 = as.numeric(Main_obj@assays[['RNA']]@scale.data['Gpt2',])
                       )
tmp_data <- cbind(Main_obj@reductions[["umap"]]@cell.embeddings, tmp_data)

write.table(tmp_data, 'scRNAseq/UMAP_1st/scExp_data_slct1.tsv', sep='\t', row.names = F, col.names = T)

# ModuleScore data -----------------------------------------------------------

tmp_list <- readLines('scRNAseq/UMAP_1st/ModuleScore/All_MitoCarta.txt')
multiple <- Reduce(intersect,  list(tmp_list, Main_obj@assays[["RNA"]]@counts@Dimnames[[1]]))
multiple <- data.frame(multiple)
Main_obj <- AddModuleScore(object=Main_obj, features=multiple, name='All_MitoCarta')

pdf(file="scRNAseq/UMAP_1st/ModuleScore/All_MitoCarta.pdf", width = 4, height = 4)
FeaturePlot(object = Main_obj,
            features = 'All_MitoCarta1',
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

tmp_list <- readLines('scRNAseq/UMAP_1st/ModuleScore/TCA_cycle.txt')
multiple <- Reduce(intersect,  list(tmp_list, Main_obj@assays[["RNA"]]@counts@Dimnames[[1]]))
multiple <- data.frame(multiple)
Main_obj <- AddModuleScore(object=Main_obj, features=multiple, name='TCA_cycle')

pdf(file="scRNAseq/UMAP_1st/ModuleScore/TCA_cycle.pdf", width = 4, height = 4)
FeaturePlot(object = Main_obj,
            features = 'TCA_cycle1',
            reduction = "umap", min.cutoff = 0,
            pt.size = 1)
dev.off()

tmp_list <- readLines('scRNAseq/UMAP_1st/ModuleScore/OXPHOS_subunits.txt')
multiple <- Reduce(intersect,  list(tmp_list, Main_obj@assays[["RNA"]]@counts@Dimnames[[1]]))
multiple <- data.frame(multiple)
Main_obj <- AddModuleScore(object=Main_obj, features=multiple, name='OXPHOS_subunits')

pdf(file="scRNAseq/UMAP_1st/ModuleScore/OXPHOS_subunits.pdf", width = 4, height = 4)
FeaturePlot(object = Main_obj,
            features = 'OXPHOS_subunits1',
            reduction = "umap", min.cutoff = 0, max.cutoff = 0.2,
            pt.size = 1)
dev.off()