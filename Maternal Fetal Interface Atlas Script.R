library(dplyr)
library(Seurat)
library(patchwork)
library(tibble)

#Loading datasets
Wang <- readRDS(file = "D:/Calgary placenta collaboration/Full integrated/Wangs_SN.rds")
garcia <- readRDS(file = "D:/Calgary placenta collaboration/Full integrated/garcia_rivas.rds")
keenen <- readRDS(file = "D:/Calgary placenta collaboration/Full integrated/Keenen_integrated.rds")
parchem <- readRDS(file = "D:/Calgary placenta collaboration/Full integrated/Parchem.rds")

#Adding Metadata Columns for Wang
Wang$Delivery_type <- Wang$type
Wang$fetal_sex <- Wang$orig.ident
Wang$fetal_sex[Wang$orig.ident == 'Late_1'] <- 'F'
Wang$fetal_sex[Wang$orig.ident == 'Late_2'] <- 'F'
Wang$fetal_sex[Wang$orig.ident == 'Late_3'] <- 'F'
Wang$fetal_sex[Wang$orig.ident == 'Late_4'] <- 'M'
Wang$fetal_sex[Wang$orig.ident == 'Late_5'] <- 'M'
Wang$fetal_sex[Wang$orig.ident == 'Late_6'] <- 'M'

#Adding Metadata Columns for keenen
keenen$Delivery_type <- keenen$type
keenen$fetal_sex <- keenen$orig.ident
keenen$fetal_sex[keenen$orig.ident == 'AGPL2_SC'] <- 'M'
keenen$fetal_sex[keenen$orig.ident == 'AGPL3_SC'] <- 'F'
keenen$fetal_sex[keenen$orig.ident == 'AGPL3_SN'] <- 'F'
keenen$fetal_sex[keenen$orig.ident == 'AGPL4_SC'] <- 'M'
keenen$fetal_sex[keenen$orig.ident == 'AGPL4_SN'] <- 'M'
keenen$fetal_sex[keenen$orig.ident == 'AGPL5_SN'] <- 'M'

#Adding metadata columns for parchem
parchem$Delivery_type <- parchem$type
parchem$fetal_sex <- 'F'

#Checking mt.percent
VlnPlot(parchem, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), group.by = 'orig.ident')
Wang <- subset(Wang, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)

#Merging all datasets together
big_boss <- merge(x = garcia, y = c(keenen, parchem, Wang))

#Pre-processing
big_boss <- NormalizeData(big_boss, normalization.method = 'LogNormalize', scale.factor = 10000)
big_boss <- FindVariableFeatures(big_boss, selection.method = 'vst', nfeatures = 2000)
top10 <- head(VariableFeatures(big_boss), 10)
LabelPoints(plot = VariableFeaturePlot(big_boss), points = top10, repel = T)

#Scaling the data
all.genes <- rownames(big_boss)
big_boss <- ScaleData(big_boss, features = all.genes)

#PCA
big_boss <- RunPCA(big_boss, features = VariableFeatures(object = big_boss))
ElbowPlot(big_boss)
DimHeatmap(big_boss, dims = 1:20, cells = 500) 

#Making the UMAP
big_boss <- FindNeighbors(big_boss, dims = 1:15)
big_boss <- FindClusters(big_boss, resolution = 1, cluster.name = 'unintegrated clusters')
big_boss <- RunUMAP(big_boss, dims = 1:15)


#Integrating the layers
big_boss <- IntegrateLayers(big_boss, method = CCAIntegration, orig.reduction = 'pca', new.reduction = 'integrated.cca')
big_boss[['RNA']] <- JoinLayers(big_boss[['RNA']])

#Running clustering analysis
big_boss <- ScaleData(big_boss)
big_boss <- FindNeighbors(big_boss, reduction = "integrated.cca", dims = 1:15)
big_boss <- FindClusters(big_boss, resolution = 1, cluster.name = "cca_clusters")

#Making the UMAp for integrated layers
big_boss <- RunUMAP(big_boss, reduction = 'integrated.cca', dims = 1:15, reduction.name = 'umap.cca')

#Dotplots
p1 <- DimPlot(big_boss, reduction = 'umap.cca', label = T, raster = F) + NoLegend() + ggtitle('Integrated UMAP')
p2 <- DimPlot(big_boss, reduction = 'umap.cca', group.by = 'orig.ident', shuffle = T, raster = F) + ggtitle('Integrated UMAP by sample')
p3 <- DimPlot(big_boss, reduction = 'umap', label = T, raster = F) + NoLegend() + ggtitle('Unintegrated UMAP')
p4 <- DimPlot(big_boss, reduction = 'umap', group.by = 'orig.ident', shuffle = T, raster = F) + ggtitle('Unintegrated UMAP by sample')
p1 + p2 + p3 + p4

#Identifying clusters
p5 <- DimPlot(big_boss, reduction = 'umap.cca', label = T)
p5 <- FeaturePlot(big_boss, features = c('PTPRC', 'CD3D', 'CD3G', 'CD86', 'CD74', 'HLA-DRA'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('PTPRC', 'EOMES', 'NCR1', 'MS4A1', 'CD19'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('PECAM1', 'PROX1', 'LYVE1', 'GP9', 'ITGA2B', 'GP1BA'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('PAPPA', 'SDC1', 'CSH2', 'ERVW-1', 'ERVFRD-1', 'GATA3'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('GATA3', 'TEAD4', 'ITGA6', 'MKI67', 'HLA-G', 'NOTCH2'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('ACTA2', 'CNN1', 'CD34', 'GJA1', 'KCNN3'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('PRL', 'IGFBP1', 'DCN', 'IGF1', 'PGR', 'HAND2'), reduction = 'umap.cca')
p5 <- FeaturePlot(big_boss, features = c('DCN', 'KRT19', 'KRT7', 'CDH1', 'HLA-G'), reduction = 'umap.cca')


#Renaming the Clusters
big_boss <- RenameIdents(object = big_boss,
                         '0' = 'ST',
                         '1' = 'ST',
                         '2' = 'ST',
                         '3' = 'ST',
                         '4' = 'Stromal 1',
                         '5' = 'pCTs',
                         '6' = 'ST',
                         '7' = 'Macrophages',
                         '8' = 'EVTs/Chorionic Trophoblasts',
                         '9' = 'Hofbauer Cells',
                         '10' = 'NK Cells',
                         '11' = 'Early ST',
                         '12' = 'Decidualizing Cells', 
                         '13' = 'T Cells',
                         '14' = 'EVTs/Chorionic Trophoblasts',
                         '15' = 'ST', 
                         '16' = 'Macrophages',
                         '17' = 'Decidualized Cells', 
                         '18' = 'Lymphatic ECs',
                         '19' = 'B Cells',
                         '20' = 'Stromal 2',
                         '21' = 'Epithelial Cells',
                         '22' = 'ST',
                         '23' = 'EVTs/Chorionic Trophoblasts', 
                         '24' = 'Fus. Comp. CTs', 
                         '25' = 'ECs',
                         '26' = 'Myocytes',
                         '27' = 'ICCs',
                         '28' = 'ECs',
                         '29' = 'pCTs',
                         '30' = 'Activated T Cells',
                         '31' = 'Hofbauer Cells',
                         '32' = 'Macrophages')
big_boss$cluster.ident <- Idents(big_boss)

big_boss$tissue[big_boss$tissue == 'Pla_SN' & big_boss$origin == 'This Manuscript'] <- 'Basal_Plate_SN'


