#note Seurat version 3.2.2 was used in our analysis. This
#may be	required to recapitulate the exact same plots and results. 
#One can refer to the script below to reproduce our filtering
#steps, etc. Alternatively, one	could simply load the
#R workspace included in our GitHub repository to explore
#our scRNA-Seq data in its final form.


#load necessary libaries

library(biomaRt)
library(Seurat)
library(reticulate)
library(Matrix)
library(ggplot2)


#load biomart data for mapping mouse and human orthologs

# Load human ensembl attributes
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Load mouse ensembl attributes
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#read filtered sequencing data output by Cell Ranger (mtx file, etc)

data <- Read10X(data.dir = "/Users/zely/Downloads/scdata/")

#create Seurat object and filter cells and genes (545 cells remain after)

#remove genes not expressed in any cell
seurat_object = CreateSeuratObject(counts = data, min.cells=1)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

#remove cells with percentage of reads derived from mitochondrial genes > 10%
seurat_2 <- subset(seurat_object, subset = percent.mt < 10)

#remove genes not expressed in remaining cells and removed cells with < 100 features detected
temp <- GetAssayData(object = seurat_2, slot = "counts")
mod_ant = CreateSeuratObject(counts = dog, min.features=100, min.cells=1)

#clean out unnecessary objects
rm(temp, seurat_object, seurat_2)


#Run standard normalization steps, PCA, etc
mod_ant <- NormalizeData(mod_ant)
mod_ant <- FindVariableFeatures(mod_ant)
mod_ant <- ScaleData(mod_ant, features=rownames(mod_ant))
mod_ant <- RunPCA(mod_ant, verbose = FALSE)

#Run standard nearest neighbor finding, cluster dectection, and UMAP construction

mod_ant <- FindNeighbors(mod_ant, dims = 1:20)
mod_ant <- FindClusters(mod_ant, resolution = 0.26)
mod_ant <- RunUMAP(mod_ant, dims = 1:20, verbose = FALSE, umap.method = "umap-learn", metric = "correlation")

#examine CD8/CD3 expression patterns (reveals a subset of CD3-/CD8- cells - likely contaminant cell type)
FeaturePlot(mod_ant, features=c("Cd8a")

#remove cells lacking expression of both Cd8a and Cd3e
mod_ant2 <- subset(x=mod_ant, subset = Cd8a > 0 & Cd3e > 0)

#remove cells exceeding 97th percentile for genes detected
mod_ant2 <- subset(mod_ant2, subset = nFeature_RNA < 4065 & nFeature_RNA > 100)
mod_ant3 <- mod_ant2

#run standard scRNA-Seq analysis steps on remaining post-QC data
mod_ant2 <- NormalizeData(mod_ant2)
mod_ant2 <- FindVariableFeatures(mod_ant2, selection.method = "vst", nfeatures = 2000)
mod_ant2 <- ScaleData(mod_ant2, features=rownames(mod_ant2))
mod_ant2 <- RunPCA(mod_ant2, verbose = FALSE, npcs = 130)

mod_ant2 <- FindNeighbors(mod_ant2, dims = 1:30)

#Note resolution = 0.69 generated 4 clusters in our environment/version of Seurat.
mod_ant2 <- FindClusters(mod_ant2, resolution = 0.69)

#Generate final UMAP for expression plots, etc
mod_ant2 <- RunUMAP(mod_ant2, dims = 1:30, verbose = FALSE, umap.method = "umap-learn", metric = "correlation")


#calculate DEGs for each clsuter
mod_ant2.markers <- FindAllMarkers(mod_ant2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)


#output raw count matrix for the final 447 post-QC cells, as input to PAGODA
seurat_object <- CreateSeuratObject(counts = data, min.features=100)
to_out <- subset(seurat_object, cells=colnames(mod_ant2))
to_out <- GetAssayData(object = to_out, slot = "counts")
write.csv(to_out, file = "/Users/zely/Downloads/resolution/prep_scde.csv",row.names = TRUE, col.names = NA)


#load published gene modules into R and the seurat object
wherry <- read.table("/Users/zely/Downloads/wherry_comp.csv", sep = ",", header=TRUE)
regev <- read.table("/Users/zely/Downloads/mmc4.csv", sep = ",", header=TRUE)
#prepare Wherry modules
AM2 <- list(wherry$Ô..AM2[2:207])
AM12 <- list(wherry$AM12[2:89])
AM8 <- list(wherry$AM8[2:17])
AM5 <- list(wherry$AM5[2:181])
AM3 <- list(wherry$AM3[2:25])
AM13 <- list(wherry$AM13[2:41])
AM7 <- c("BC050254,Btg2,Cpne3,Dusp1,Ier2,Pea15a,Ptpn4,Pvt1,Rabgap1l,Rbl2,Rbm38,Rhob,Rnf138,Rybp,Spry2,St8sia4,Syt11,Taf1d,Zbtb1")
AM7 <- as.list(strsplit(AM7,","))
AM9 <- list(wherry$AM9[2:266])
AM4 <- list(wherry$AM4[2:30])
AM10 <- list(wherry$AM10[2:282])
AM1 <- list(wherry$AM1[2:49])
AM11 <- list(wherry$AM11[2:77])
AM6 <- list(wherry$AM6[2:52])
CM2 <- list(wherry$CM2[2:258])
CM5 <- list(wherry$CM5[2:681])
CM3 <- list(wherry$CM3[2:70])
CM1 <- list(wherry$CM1[2:245])
CM6 <- list(wherry$CM6[2:249])
CM4 <- list(wherry$CM4[2:16])

#Prepare Davis modules
davis <- c("Samd3,Plek,Ccl5,Ctla2a,Klrb1c,Klrk1,Lrrk1,Cxcr3,Serpina3g,Klra7,Gpr15,Cd44,Ctla2b,Osbpl3,1700025G04Rik,Il18rap,Ly6c2,Klrc1,Ppm$
davis <- as.list(strsplit(davis, ","))
davis_updted <- c("Klra3,Klra7,Klre1,Ccr5,Klra23,Klra6,Klrc2,Klra1,S1pr5,Klrc1,Klrb1c,Klrk1,Samd3,Ccl5,Plek,Dmrta1,Olfm1,Lrrk1,Ccr2,Ctla2b$
davis_updted <- as.list(strsplit(davis_updted, ","))

#Prepare Regev modules
#note these modules are derived from human genes. We convert these to mouse
#orthologs in the loops below.

act_module <- list(regev$Ô..Activation_module_100[2:100])
for (i in act_module)
{
genes.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = i, mart = human, attributesL = c("mgi_symbol"), martL = mouse)
}
act_module <- genes.list
act_module <- list(act_module$MGI.symbol)

dys_module <- list(regev$Dysfunction_module_100[2:100])
for (i in dys_module)
{
genes2.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = i, mart = human, attributesL = c("mgi_symbol"), martL = mouse)
}
dys_module <- genes2.list
dys_module <- list(dys_module$MGI.symbol)

act_dys_module <- list(regev$Activation_Dysfunction_module_100[2:100])
for (i in act_dys_module)
{
genes3.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = i, mart = human, attributesL = c("mgi_symbol"), martL = mouse)
}
act_dys_module <- genes3.list
act_dys_module <- list(act_dys_module$MGI.symbol)

naive_mem_module <- list(regev$Naive_Memory_like_module_100[2:100])
for (i in naive_mem_module)
{
genes4.list = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = i, mart = human, attributesL = c("mgi_symbol"), martL = mouse)
}
naive_mem_module <- genes4.list
naive_mem_module <- list(naive_mem_module$MGI.symbol)


#Prepare PAGODA modules derived from our de novo module analysis (see related folder in GitHub repository)
pagoda_36 <- c("Slc17a6,Zbtb32,Mt1,Map2k3,Irf4,Tnip3,Sdcbp2,Twsg1,Slc22a15,Havcr2,Ccl3,Elk3,Spry1,Eea1,Nrp1,Ubash3b,Mir155hg,Alcam,Irf8,Tnfrsf9,Tnfrsf4,Mdfic,Capg,Cish,Csf1,Itgav,Rbpj,Samsn1,Hip1,Stk39,Ptms,Lag3,Pdcd1,Casp3,Ier5l,Chst2,Tox,Ikzf2,2010111I01Rik,Ctsb,Bcl2a1b,Bcl2a1d,Gapdh,Nr4a2,Rgs16,Bhlhe40,Ctla4")
pagoda_36 <- as.list(strsplit(pagoda_36,","))
pagoda_30 <- c("Tmsb10,Tcf7,Rgs10,Emb,Evl,Pik3ip1,Id3,Rapgef6,Foxp1,Cdkn2d,Il7r,Wdr89,Klk8,Satb1,Sell,Lef1,Ccr7,Ighm,Igfbp4,Dapl1,Actn1,Nsg2,Il6ra,Cmah,Slco3a1,Arl4c,Hsd11b1,Sidt1,Ly6c2,Ms4a4c,Ripor2,Dph5,Gm2682,Rasgrp2,Dtx1,Lcn4,Treml2,Gramd4,Tdrp,Rflnb,Atp1b1,Dennd2d,Cd2ap,Dgka,Kbtbd11,S1pr1,Klf2,Ssh2,Socs3,Tubb5,Hmgn1")
pagoda_30 <- as.list(strsplit(pagoda_30,","))
pagoda_45 <- c("Padi2,Trps1,Cd38,Lrrk1,Ccdc50,Itgb1,Tnfrsf1b,Runx2,Nedd9,Pglyrp1,Efhd2,Unc119,Fuca2,Nmb,Endod1,Zbtb38,Gpr68,Acadl,Serpinb6a,Plscr1,Gabarapl1,Cyfip1,Casp4,Nsmaf,Itga4,Ptpn11,Ric1,Ppp1r16b,Ttc39c,Sytl2,Neat1,Nabp1,Evi2a,Gna15,Il18r1,Il18rap,Myo1f,Gzma,Plek,Zeb2,Cx3cr1,Bicral,Gxylt1,Coch,Maf,Pear1,L1cam,Rpa2,Cers4,Klrg1,Ccr2,Ccr5,Gzmb,Eomes,Gzmk,Hip1r,Kif13b,Ciapin1,Mxi1,Stk32c,Azin1,Ncald,Lims1,Lilrb4a,Tm6sf1,Dok2,Lxn,Dapk2,Vmp1,Ikzf3,Fgfr1op2,Sh3glb1,Ppp3ca,Mapre2,Rab8b,Mif4gd,Lmnb1,1700017B05Rik,Cenpa,Baz1a,N4bp1,Mfsd14a,Prdm1,Rasgef1b,Entpd1,S100a11,Serpina3g,Tmsb4x,Anxa2,Txn1,Myl6,Adam19,AW112010,Glrx,Gng2,Prr13,Ctsd,Cox8a,Cox17,H2-D1,H2-K1,B2m,H2-Q7,Cd3g,Cd3e,Cd8a,Osbpl3,Fasl,Cdk6,Acot7,Chsy1,Camk2n1,Lilr4b,Litaf,Cxcr6,Ybx3,Atxn1,Ccl4,Ccl5,Ctla2a,Nkg7,Il2rb,S100a4,S100a6,Lgals3,Ahnak,Lgals1,Ifng,Rgs1,Id2")
pagoda_45 <- as.list(strsplit(pagoda_45,","))


#Add modules to Seurat object
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM1, ctrl = 8, name = "AM1")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM2, ctrl = 8, name = "AM2")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM3, ctrl = 8, name = "AM3")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM4, ctrl = 8, name = "AM4")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM5, ctrl = 8, name = "AM5")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM6, ctrl = 8, name = "AM6")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM7, ctrl = 8, name = "AM7")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM8, ctrl = 8, name = "AM8")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM9, ctrl = 8, name = "AM9")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM10, ctrl = 8, name = "AM10")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM11, ctrl = 8, name = "AM11")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM12, ctrl = 8, name = "AM12")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = AM13, ctrl = 8, name = "AM13")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = CM1, ctrl = 8, name = "CM1")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = CM2, ctrl = 8, name = "CM2")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = CM3, ctrl = 8, name = "CM3")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = CM4, ctrl = 8, name = "CM4")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = CM5, ctrl = 8, name = "CM5")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = CM6, ctrl = 8, name = "CM6")

mod_ant2 <- AddModuleScore(object = mod_ant2, features = act_module, ctrl = 8, name = "Activation_module")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = dys_module, ctrl = 8, name = "Dysfunction_module")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = act_dys_module, ctrl = 8, name = "Active_dysfunct_module")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = naive_mem_module, ctrl = 8, name = "Naive_mem_module")

mod_ant2 <- AddModuleScore(object = mod_ant2, features = pagoda_36, ctrl = 8, name = "PAGODA_36")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = pagoda_30, ctrl = 8, name = "PAGODA_30")
mod_ant2 <- AddModuleScore(object = mod_ant2, features = pagoda_45, ctrl = 8, name = "PAGODA_45")

#Output heatmap with top genes
setEPS()
postscript('/Users/zely/Downloads/heatmap.eps', height=6.8, width=6.8)
DoHeatmap(mod_ant2, features=c("S100a4","S100a6","Pdcd1","Lag3","Ccl4","Bhlhe40","Rgs16","Lgals3","Gzmb","Fasl","Klrc1","Ccl5","Bcl2a1d","Nr4a2","Gzmk","Litaf","Ccl3","Ptms","Klrk1","Ahnak","Camk2n1","Ctla4","Lilr4b","Hip1","Anxa2","Socs3","Dapl1","Igfbp4","Sbno2","Fam241a","Tdrp","Ccr7","Sell","Fos","Satb1","Maff","Zfp36","Socs1","Ppp1r15a","Klf2","Gm20186","Csrnp1","Nsg2","Lef1","Ssh2","Myc","Gm43698","Gm14085","Gadd45b","Irs2","Pdk1","Igfbp4","Sell","Ramp1","Dapl1","Lef1","Dtx1","Ighm","Actn1","Art2b","Treml2","Il6ra","Ccr9","Tmem108","Cd79b","Klra7","Klra6","1700025G04Rik","Ly6c2","Samd3","Klra9","Rbpms","Slc11a2","F2rl1","Yes1","Ms4a4c","Gas7","Gzmm","Clec2i","Xcl1","Rps20"))+theme(axis.text.y=element_text(size=6,face="bold"))+theme(axis.ticks.y.left=element_line(size=0.5))
dev.off()



#output EPS plots depicting gene module expression (a few examples)
#Note that (in our version) Seurat added "1" to the names of added modules. 
#So, the feature plot reference (e.g.,) "AM81" plots the expression of 
#module, AM8. 
setEPS()
postscript("/Users/zely/Downloads/resolution/Active_dysfunct_module.eps", height=6.8, width=6.8)
FeaturePlot(mod_ant2, features = c("Active_dysfunct_module1"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1,pt.size=1.6, order=TRUE)+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()

postscript("/Users/zely/Downloads/resolution/am8.eps", height=6.8, width=6.8)
FeaturePlot(mod_ant2, features = c("AM81"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1,pt.size=1.6, order=TRUE)+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
setEPS()

postscript("/Users/zely/Downloads/resolution/pagoda_36.eps", height=6.8, width=6.8)
FeaturePlot(mod_ant2, features = c("PAGODA_361"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1,pt.size=1.6, order=TRUE)+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()

postscript("/Users/zely/Downloads/resolution/pagoda_45.eps", height=6.8, width=6.8)
FeaturePlot(mod_ant2, features = c("PAGODA_451"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1,pt.size=1.6, order=TRUE)+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()
postscript("/Users/zely/Downloads/resolution/pagoda_30.eps", height=6.8, width=6.8)
FeaturePlot(mod_ant2, features = c("PAGODA_301"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1,pt.size=1.6, order=TRUE)+scale_colour_gradient(low = "cadetblue2", high = "darkorange")+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()


#output umap plot without labels
setEPS()
postscript("/Users/zely/Downloads/resolution/umap_nolabel.eps", height=6.8, width=6.8)
DimPlot(mod_ant2, label = FALSE, pt.size=1.6) + NoLegend()+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)
dev.off()

#output gene expression plots (one example)
postscript('/Users/zely/Downloads/resolution/cd3_pattern.eps', height=6.8, width=6.8)
FeaturePlot(mod_ant2, features = c("Cd3e"), min.cutoff = "q01", max.cutoff = "q99", ncol = 1,pt.size=1.6, order=TRUE)+FontSize(x.title=34, y.title=34, main=40, x.text=32, y.text=32)+ theme(legend.text=element_text(size=24))
dev.off()

#output raw count matrix of final cells for GEO
proc_count <- GetAssayData(object = mod_ant2, slot = "counts")
write.csv(proc_count, file = "/Users/zely/Downloads/resolution/processed_data_count_matrix.csv",row.names = TRUE, col.names = NA)

#output normalized data for GEO
proc_norm <- GetAssayData(object = mod_ant2, slot = "data")
write.csv(proc_norm, file = "/Users/zely/Downloads/resolution/processed_data_normalized_matrix.csv",row.names = TRUE, col.names = NA)
