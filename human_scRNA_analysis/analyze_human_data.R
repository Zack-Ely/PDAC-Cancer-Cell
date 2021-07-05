#This script can be used as a guide to recapitulate results from our analysis
#of human PDAC data from Peng et al (2019) (see ms for citation). Please note
#that the same version of Seurat (v3.2.2) may be necessary to reproduce
#the exact same plots,etc from our manuscript. We have included a R workspace
#called human_TIL_workspace.RData. The most relevant and final Seurat object in this
#workspace is called, patient.integrated
#This workspace can be loaded into R to examine our analysis
#in its final form. Finally, note that the raw count data from peng et al (2019)
#is publicly available.

#load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(reticulate)
library(Matrix)
library(hdf5r)
library(biomaRt)
library(magrittr)
library(dplyr)

# Load human ensembl attributes
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Load mouse ensembl attributes
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#load raw count data from peng et al PDAC scRNA-Seq study
countData <- read.table('/Users/zack/human_scrna_yang/final_yang_matrix.txt', header = TRUE, sep="\t")

#remove normal samples from the dataset
g2 <- countData[, -grep("N", colnames(countData))]
nmes <- g2[,1]
count_data <- countData[,2:41987]
rownames(count_data) <- nmes

#load all tumor-sample-derived cells across all 24 patients, into 1 Seurat object
pdac_hum <- CreateSeuratObject(counts=count_data, min.cells=10)

#run the standard scRNA-Seq workflow
pdac_hum <- NormalizeData(pdac_hum)
pdac_hum <- FindVariableFeatures(pdac_hum)
pdac_hum <- ScaleData(pdac_hum)
pdac_hum <- RunPCA(pdac_hum, verbose = FALSE)
ElbowPlot(pdac_hum, ndims=50)
pdac_hum <- FindNeighbors(pdac_hum, dims = 1:15)
pdac_hum <- FindClusters(pdac_hum, resolution = 1)
pdac_hum <- RunUMAP(pdac_hum, dims = 1:15, verbose = FALSE, metric = "euclidean")

#output some plots for the whole dataset
setEPS()
postscript("/Users/zack/NEW_AUG_2020/whole_umap.eps", height=6.8, width=6.8)
DimPlot(pdac_hum, label=TRUE)+ NoLegend()+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zack/cd8a_whole.eps", height=6.8, width=6.8)
FeaturePlot(pdac_hum, features = c("CD8A"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE)+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zack/cd3e_whole.eps", height=6.8, width=6.8)
FeaturePlot(pdac_hum, features = c("CD3E"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE)+ NoLegend()+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()

#subset cd4s and cd8s based on expression of major CD3 genes and CD8A and CD4
cd8es <- subset(x = pdac_hum, subset = CD3E > 0.1 & CD8A > 0.1)
cd8ds <- subset(x = pdac_hum, subset = CD3D > 0.1 & CD8A > 0.1)
cd8gs <- subset(x = pdac_hum, subset = CD3G > 0.1 & CD8A > 0.1)
allcd8_d <- setdiff(colnames(cd8ds), colnames(cd8es))
allcd8_g <- setdiff(colnames(cd8gs), c(colnames(cd8es),colnames(cd8ds)))
all_cd8 <- c(colnames(cd8es), allcd8_d, allcd8_g)


cd4es <- subset(x = pdac_hum, subset = CD3E > 0.1 & CD4 > 0.1)
cd4ds <- subset(x = pdac_hum, subset = CD3D > 0.1 & CD4 > 0.1)
cd4gs <- subset(x = pdac_hum, subset = CD3G > 0.1 & CD4 > 0.1)
allcd4_d <- setdiff(colnames(cd4ds), colnames(cd4es))
allcd4_g <- setdiff(colnames(cd4gs), c(colnames(cd4es),colnames(cd4ds)))
all_cd4 <- c(colnames(cd4es), allcd4_d, allcd4_g)

full4 <- setdiff(all_cd4, all_cd8)
full8 <- setdiff(all_cd8, all_cd4)

cts <- GetAssayData(object = pdac_hum, slot = "counts")
cd4_2 <- cts[, full4]
cd8_2 <- cts[, full8]
cd4_3 <- CreateSeuratObject(counts=cd4_2, min.cells=5)
cd8_3 <- CreateSeuratObject(counts=cd8_2, min.cells=5)

#merge all T cells
tcells <- merge(x=cd8_3, y=cd4_3)

#examine number of candidate T cells across all patients (note each "TX" refers to a unique tumor sample).
counts <- table(tcells[[]]$orig.ident)
count_test <- counts[order(-counts)]
postscript("/Users/zack/barplot_cd3_pos_cells_per_patient.eps", height=6.8, width=6.8)
barplot(count_test, main="CD3+ Cells Per Patient", xlab="Patient Number", cex.names=0.45, ylim=c(0,400), font.axis=2, cex.lab=1.4, cex.axis=1.1)
dev.off()


#subset on patients with > 50 T cells (i.e., exclude pts w < 50)
tcells_20 <- subset(x = tcells, subset = orig.ident == "T20", invert=TRUE)
tcells_20 <- subset(x = tcells_20, subset = orig.ident == "T18", invert=TRUE)
tcells_20 <- subset(x = tcells_20, subset = orig.ident == "T17", invert=TRUE)
tcells_20 <- subset(x = tcells_20, subset = orig.ident == "T14", invert=TRUE)
tcells_2 <- subset(x = tcells_20, subset = orig.ident == "T9", invert=TRUE)


#Run the Seurat integration workflow to resolve batch effects among remaining samples' T cells
#See the corresponding Seurat vignette for more information

patient.list <- SplitObject(tcells_2, split.by = "orig.ident")
#run SCTransform for each patient-specific Seurat object
for (i in 1:length(patient.list)) {
  patient.list[[i]] <- SCTransform(patient.list[[i]], verbose = FALSE)
}

#identify features to use in the Integration function and prepare the objects for integration
patient.features <- SelectIntegrationFeatures(object.list = patient.list, nfeatures = 3000)
options(future.globals.maxSize= 6891289600)
patient.list <- PrepSCTIntegration(object.list = patient.list, anchor.features = patient.features, verbose = FALSE)
#save.image(file="~/tcell_analysis.RData")

patient.anchors <- FindIntegrationAnchors(object.list = patient.list, normalization.method = "SCT", anchor.features = patient.features, verbose = FALSE, k.filter=50)
patient.integrated <- IntegrateData(anchorset = patient.anchors, normalization.method = "SCT", verbose = FALSE)

#Run standard scRNA-Seq workflow starting with PCA
patient.integrated <- RunPCA(patient.integrated, verbose = FALSE)
ElbowPlot(patient.integrated)

patient.integrated <- FindNeighbors(patient.integrated, dims = 1:14)
patient.integrated <- FindClusters(patient.integrated, resolution = 0.8)
patient.integrated <- RunUMAP(patient.integrated, dims = 1:14)

#output UMAPs displaying the batch-corrected CD3+ CD4+/CD8+ cells
setEPS()
postscript("/Users/zack/tcell_integrated.eps", height=6.8, width=6.8)
DimPlot(patient.integrated, label = TRUE, reduction="umap") + NoLegend()+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zack/tcell_integrated_by_patient.eps", height=6.8, width=6.8)
DimPlot(patient.integrated, reduction="umap",group.by = c("orig.ident"))+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()

#change to RNA assay and prep data to generate gene expression plots
DefaultAssay(patient.integrated) <- "RNA"
patient.integrated <- NormalizeData(patient.integrated)
patient.integrated <- FindVariableFeatures(patient.integrated)
patient.integrated <- ScaleData(patient.integrated)
#patient.integrated <- RunPCA(patient.integrated, verbose = FALSE)
#ElbowPlot(patient.integrated, ndims=50)

patient.integrated.markers <- FindAllMarkers(patient.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

#output gene expression plots (2 examples)
postscript("/Users/zack/integrated_cd4.eps", height=6.8, width=6.8)
FeaturePlot(patient.integrated, features = c("CD4"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE, reduction="umap")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zack/integrated_cd8a.eps", height=6.8, width=6.8)
FeaturePlot(patient.integrated, features = c("CD8A"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE, reduction="umap")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()


#proportion of patients CD3+ cells that are in each cluster
propo_cd3 <- prop.table(table(Idents(patient.integrated), patient.integrated$orig.ident), margin = 2)
write.csv(propo_cd3, file = "/Users/zack/proportion_of_patient_cd3_cells_in_each_cluster.csv",row.names = TRUE, col.names = NA)



#Establish gene lists for PAGODA modules derived from mouse TIL analyses
full_pg_45 <- c("Padi2,Trps1,Cd38,Lrrk1,Ccdc50,Itgb1,Tnfrsf1b,Runx2,Nedd9,Pglyrp1,Efhd2,Unc119,Fuca2,Nmb,Endod1,Zbtb38,Gpr68,Acadl,Serpinb6a,Plscr1,Gabarapl1,Cyfip1,Casp4,Nsmaf,Itga4,Ptpn11,Ric1,Ppp1r16b,Ttc39c,Sytl2,Neat1,Nabp1,Evi2a,Gna15,Il18r1,Il18rap,Myo1f,Gzma,Plek,Zeb2,Cx3cr1,Bicral,Gxylt1,Coch,Maf,Pear1,L1cam,Rpa2,Cers4,Klrg1,Ccr2,Ccr5,Gzmb,Eomes,Gzmk,Hip1r,Kif13b,Ciapin1,Mxi1,Stk32c,Azin1,Ncald,Lims1,Lilrb4a,Tm6sf1,Dok2,Lxn,Dapk2,Vmp1,Ikzf3,Fgfr1op2,Sh3glb1,Ppp3ca,Mapre2,Rab8b,Mif4gd,Lmnb1,1700017B05Rik,Cenpa,Baz1a,N4bp1,Mfsd14a,Prdm1,Rasgef1b,Entpd1,S100a11,Serpina3g,Tmsb4x,Anxa2,Txn1,Myl6,Adam19,AW112010,Glrx,Gng2,Prr13,Ctsd,Cox8a,Cox17,H2-D1,H2-K1,B2m,H2-Q7,Cd3g,Cd3e,Cd8a,Osbpl3,Fasl,Cdk6,Acot7,Chsy1,Camk2n1,Lilr4b,Litaf,Cxcr6,Ybx3,Atxn1,Ccl4,Ccl5,Ctla2a,Nkg7,Il2rb,S100a4,S100a6,Lgals3,Ahnak,Lgals1,Ifng,Rgs1,Id2")
full_pg_45 <- as.list(strsplit(full_pg_45,","))

full_pg_30 <- c("Tmsb10,Tcf7,Rgs10,Emb,Evl,Pik3ip1,Id3,Rapgef6,Foxp1,Cdkn2d,Il7r,Wdr89,Klk8,Satb1,Sell,Lef1,Ccr7,Ighm,Igfbp4,Dapl1,Actn1,Nsg2,Il6ra,Cmah,Slco3a1,Arl4c,Hsd11b1,Sidt1,Ly6c2,Ms4a4c,Ripor2,Dph5,Gm2682,Rasgrp2,Dtx1,Lcn4,Treml2,Gramd4,Tdrp,Rflnb,Atp1b1,Dennd2d,Cd2ap,Dgka,Kbtbd11,S1pr1,Klf2,Ssh2,Socs3,Tubb5,Hmgn1")
full_pg_30 <- as.list(strsplit(full_pg_30,","))

full_pg_36 <- c("Slc17a6,Zbtb32,Mt1,Map2k3,Irf4,Tnip3,Sdcbp2,Twsg1,Slc22a15,Havcr2,Ccl3,Elk3,Spry1,Eea1,Nrp1,Ubash3b,Mir155hg,Alcam,Irf8,Tnfrsf9,Tnfrsf4,Mdfic,Capg,Cish,Csf1,Itgav,Rbpj,Samsn1,Hip1,Stk39,Ptms,Lag3,Pdcd1,Casp3,Ier5l,Chst2,Tox,Ikzf2,2010111I01Rik,Ctsb,Bcl2a1b,Bcl2a1d,Gapdh,Nr4a2,Rgs16,Bhlhe40,Ctla4")
full_pg_36 <- as.list(strsplit(full_pg_36,","))


#convert PAGODA gene sets to sets comprised of human orthologs
for (i in full_pg_45)
{
  genes.list = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = i, mart = mouse, attributesL = c("hgnc_symbol"), martL = human)
}
pg_45_full <- genes.list
pg_45_full_module <- list(pg_45_full$HGNC.symbol)

for (i in full_pg_30)
{
  genes.list = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = i, mart = mouse, attributesL = c("hgnc_symbol"), martL = human)
}
pg_30_full <- genes.list
pg_30_full_module <- list(pg_30_full$HGNC.symbol)

for (i in full_pg_36)
{
  genes.list = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = i, mart = mouse, attributesL = c("hgnc_symbol"), martL = human)
}
pg_36_full <- genes.list
pg_36_full_module <- list(pg_36_full$HGNC.symbol)

#Add PAGODA module scores to batch-corrected Seurat object
patient.integrated <- AddModuleScore(object = patient.integrated, features = pg_30_full_module, ctrl = 8, name = "Pagoda_30_Full")
patient.integrated <- AddModuleScore(object = patient.integrated, features = pg_36_full_module, ctrl = 8, name = "Pagoda_36_Full")
patient.integrated <- AddModuleScore(object = patient.integrated, features = pg_45_full_module, ctrl = 8, name = "Pagoda_45_Full")

#output plots of PAGODA module expression
postscript("/Users/zack/integrated_pagoda36.eps", height=6.8, width=6.8)
FeaturePlot(patient.integrated, features = c("Pagoda_36_Full1"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE, reduction="umap")+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zack/integrated_pagoda45.eps", height=6.8, width=6.8)
FeaturePlot(patient.integrated, features = c("Pagoda_45_Full1"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE, reduction="umap")+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zack/integrated_pagoda30.eps", height=6.8, width=6.8)
FeaturePlot(patient.integrated, features = c("Pagoda_30_Full1"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1, order=TRUE, reduction="umap")+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
