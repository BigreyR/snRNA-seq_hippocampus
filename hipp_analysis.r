#! /usr/bin/env Rscript
# library -----------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(scater)
library(Matrix)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(destiny)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(monocle)
library(circlize)
library(ggplot2)
library(Seurat)
library(scater)
library(Matrix)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(destiny)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(yyplot)
library(devtools)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(pheatmap)
library(hrbrthemes)
library(ggrepel)
library(ggVennDiagram)
library(RVenn)
library(SCENIC)
library(FateID)
library(Seurat)
library(devtools)
library(velocyto.R)
library(SeuratDisk)
library(SeuratObject)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(loomR)
library(hdf5r)
library(dyno)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(scater)
library(Matrix)
library(cowplot)
library(SingleCellExperiment)
library(scater)
library(destiny)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(yyplot)
library(devtools)
library(ggplotify)
library(PCAtools)
library(corrplot)
library(DoubletFinder)
library(reshape2)
library(ComplexHeatmap)
library(UpSetR)
library(future)
library(paletteer)
library(data.table)
library(tidydr)
library(clusterProfiler)
library(SingleR)
library(ggpubr)
library(Augur)
library(gg.gap)
library(SCopeLoomR)

# — 1 Fig1: Single-nucleus transcriptomic atlas of the human hippocampus across different ages and after stroke injury -----------------------------------------------------------
# —— 1.1 Fig1B Atlas umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))
# celltype
cell<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

(p1<-DimPlot(tmp_Sobj, reduction = "umap",cols = color16,label.size = 6,raster = T,
             label = T, pt.size = 1,group.by = "Major",repel = T) + 
    # theme_bw()+
    theme(text = element_text(size = 15),legend.position = "none",
          plot.title = element_blank(),title = element_blank(),
          axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_blank(),panel.grid = element_blank(),
          strip.switch.pad.grid = element_blank())+NoLegend())
ggsave("figure/Fig1B.Atlas_umap.pdf",p1,width = 9,height = 9)

# —— 1.2 Fig1C Cell marker bubble plot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))
# celltype
cell<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

DefaultAssay(tmp_Sobj) <- "RNA"
table(tmp_Sobj$Celltype)
gene<-c("ALDH1L1","GFAP","HOPX","VIM","PAX6","SOX2","CCND2","SOX4","STMN2","SOX11","PROX1","NEUROD2","SYT1","SV2B","SYN1","GAD1","CCK","VIP","OLIG1","OLIG2","SOX10","MOG","MBP","PLP1","CSF1R","CTSS","VWF","PNN","MAP3K15","RELN","TBX18")
celltype<-c("M-AS","AS/qNSC","pNSC","aNSC","NB","GC","IN","OPC","OLG","MG","EC","Pyr","CR","Per")
sub<-subset(tmp_Sobj,subset = Celltype %in% celltype) %>% 
  {.$Celltype<-factor(.$Celltype,levels=celltype);.}

pp<-DotPlot(object = sub, features = gene, group.by  = "Celltype",scale = T,scale.by = "size",col.min=-1,col.max = 3,
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>% 
  {.$Cluster<-factor(.$Cluster,levels = rev(celltype));.}
(p1<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(0,10))+
    scale_color_gradient2(low="grey90",mid="#FFE162",high ="#EA2027",
                          midpoint = mean(c(min(result$`Average Expression`),max(result$`Average Expression`))))+
    theme_classic()+
    theme(axis.text=element_text(size=15, color="black"),legend.position = "top",axis.title = element_blank()) + RotatedAxis())
ggsave("figure/Fig1C.Cell_marker.bubble_plot.pdf",p1,width = 12,height = 6.5)

# —— 1.3 Fig1D Cell marker feature plot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))
# celltype
cell<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

DefaultAssay(tmp_Sobj) <- "RNA"
gene<-c("ALDH1L1","GFAP","PAX6","VIM","STMN2","SV2B","OLIG1","CSF1R","GAD1","RELN","MAP3K15","VWF")
(p<-FeaturePlot(tmp_Sobj, features = gene,repel = T,pt.size = 0.5,ncol = 6,raster = T,
                reduction = "umap",label = F) &
    scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")))&
    scale_y_continuous(breaks=NULL)&
    scale_x_continuous(breaks=NULL)&
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
          axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank()))

ggsave("figure/Fig1D.Cell_marker.FeaturePlot.pdf",p,width = 24,height = 8)

# —— 1.4 Fig1E Atlas umap split by Stage -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))
# celltype
cell<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

(p1<-DimPlot(tmp_Sobj, reduction = "umap",cols = color16,label.size = 6,raster = T,
             label = F, pt.size = 1,group.by = "Major",repel = T) + 
    # theme_bw()+
    theme(text = element_text(size = 15),legend.position = "none",
          plot.title = element_blank(),title = element_blank(),
          axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_blank(),panel.grid = element_blank(),
          strip.switch.pad.grid = element_blank())+
    facet_wrap(~tmp_Sobj$Stage,nrow = 1)+NoLegend())
ggsave("figure/Fig1E.Atlas_umap.split_by_Stage.pdf",p1,width = 16,height = 6)

# — 2 FigS1 Related to Fig1: Cell atlas of human hippocampus across different ages and post stoke-induced injury -----------------------------------------------------------
# —— 2.1 FigS1AB QC Violin -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))

load("../result/Before_filtered.alldata.RData")
load("../result/seurat_filtered.v1.RData")

alldata$orig.ident<-factor(alldata$orig.ident,levels = levels(tmp_Sobj$orig.ident))
p1<-VlnPlot(alldata, features = c("nFeature_RNA","nCount_RNA","percent.mito"),
            group.by = "orig.ident",pt.size = 0)&
  NoLegend()&
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14,colour = "black"))
data.filt$orig.ident<-factor(data.filt$orig.ident,levels = levels(tmp_Sobj$orig.ident))
p2<-VlnPlot(data.filt, features = c("nFeature_RNA","nCount_RNA","percent.mito"),
            group.by = "orig.ident",pt.size = 0)&
  NoLegend()&
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14,colour = "black"))

(p<-p1/p2)
ggsave("figure/Fig1E.Atlas_umap.split_by_Stage.pdf",p,width = 17,height = 8)


# —— 2.3 FigS1C Atlas umap(3D) -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))
# celltype
cell<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

library(plotly)

# Construct a dataframe using data from your pre-clustered Seurat v3.1.1 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.1.1 (works for v3.0.0 as well)
yourseuratobject <- tmp_Sobj
Idents(tmp_Sobj)
# Re-run UMAPs that you have accurate calculations for all UMAP(s)
yourseuratobject <- RunUMAP(yourseuratobject,seed.use = 42,
                            dims = 1:20,
                            n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- yourseuratobject[["umap"]]@cell.embeddings[,1]
umap_2 <- yourseuratobject[["umap"]]@cell.embeddings[,2]
umap_3 <- yourseuratobject[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Celltype"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~Celltype, 
        colors = color16,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 2, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text")

# —— 2.4 FigS1D Atlas umap color by samples -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

color10<-c("#4DBBD5B2","#00A087B2","#FF82AB","#9196f7","#6ab04c","#DC0000B2",
           "#FFA500","#40407a","#FF7256","#7E6148B2")

(p1<-DimPlot(tmp_Sobj, reduction = "umap",cols = color10,label.size = 6,raster = T,
             label = F, pt.size = 0.5,group.by = "orig.ident",repel = T) + 
    # theme_bw()+
    theme(text = element_text(size = 15),
          plot.title = element_blank(),title = element_blank(),
          axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_blank(),panel.grid = element_blank(),
          strip.switch.pad.grid = element_blank()))
ggsave("figure/FigS1D.Atlas_umap.color_by_samples.pdf",p1,width = 10,height = 9)

# —— 2.5 FigS1E Top marker heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

top=50
mark_gene <- read.delim("data/FigS1E.all.markers.celltype.v1.xls", stringsAsFactors=FALSE) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC) %>% .$gene
cell_type<-c("Granule cells","Interneuron","Neuroblast","aNSC","Intermediate NSC",
             "Astrocytes-qNSC","Mature astrocytes","OPC","Oligodendrocytes","Microglia",
             "Endothelial cells","Pyramidal neuron","CR","Pericytes","Unknown1","Unknown2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

DefaultAssay(tmp_Sobj) <- "RNA"
mat <- tmp_Sobj@assays[["RNA"]]@data

cluster_info <- sort(tmp_Sobj$seurat_cluster_cell_typer)
cluster_info <- factor(cluster_info, levels = cell_type)
mat <- as.matrix(mat[mark_gene, names(cluster_info)])
mat[mat>4] =4

df <- as.data.frame(sort(tmp_Sobj$seurat_cluster_cell_typer),stringsAsFactors=F) %>% 
  {colnames(.)<-"cluster";.} %>% {.$cluster<-factor(.$cluster,levels=cell_type);.}
top_ann = HeatmapAnnotation(type = df$cluster,show_annotation_name = F,show_legend = F,
                            col = list(type = c("Granule cells"="#4DBBD5B2","Interneuron"="#6ab04c",
                                                "Neuroblast"="#9196f7","aNSC"="#DC0000B2",
                                                "Intermediate NSC"="#00A087B2","Astrocytes-qNSC"="#40407a",
                                                "Mature astrocytes"="#FF82AB","OPC"="#FFA500",
                                                "Oligodendrocytes"="#FF7256","Microglia"="#4682B4",
                                                "Endothelial cells"="#8B658B","Pyramidal neuron"="#FFFF00",
                                                "CR"="#00F5FF","Pericytes"="#000000",
                                                "Unknown1"="#7E6148B2","Unknown2"="#7FFF00")))

library(circlize)
f2 = colorRamp2(seq(min(mat), max(mat), length = 5),c("#0B5FA5","#3F8FD2","#FFE4C4","#FF7D73","red"),space = "RGB")
ht<-Heatmap(mat,
            col = f2,
            cluster_rows = F,
            cluster_columns = F,
            show_column_names = FALSE,
            show_row_names = T,
            column_split = cluster_info,
            column_title = NULL,
            # column_gap=0.1,
            top_annotation = top_ann,
            use_raster = F,
            # right_annotation = right_ann,
            row_names_gp = gpar(fontsize = 2),
            column_names_gp = gpar(fontsize = 2),
            heatmap_legend_param = list(
              title = "LogNormalize",
              title_position = "leftcenter-rot"),
            width = unit(40, "cm"), height = unit(20, "cm")) 
filename<-paste("figure/FigS1E.Top_marker.heatmap.png",sep = "")
png(file = filename,width = 5400,height = 2700,res=300)
ht
dev.off()

# —— 2.6 FigS1F Mean gene summary by celltype -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))
# celltype
cell<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")

(p1<-tmp_Sobj@meta.data %>% .[,c("nCount_RNA","nFeature_RNA","Major","Stage")] %>% 
    group_by(Stage,Major) %>% dplyr::summarise(Mean_genes = mean(nFeature_RNA),
                                               Mean_transcripts=mean(nCount_RNA)) %>% 
    as.data.frame() %>% 
    ggplot(data = ., mapping = aes(x = Stage, y = Mean_genes,fill=Stage)) + 
    geom_bar(stat= 'identity', position = 'dodge') +
    # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
    #               position=position_dodge(.9))+
    # scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")) +
    scale_fill_manual(values=c("#f5d662","#EB5230","#800026","#2A4E9D")) +
    ylab(label = "Mean genes") + xlab(label = "")+
    scale_x_discrete(expand = c(0, 0.6)) + 
    scale_y_continuous(expand = c(0, 0.002),limits=c(0,6200)) +
    # scale_fill_npg()+
    #ylim(0,0.35) + 
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.x = element_text(angle = 30,hjust = 1),
          legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")) +
    facet_wrap(~Major,nrow = 2)+NoLegend())
ggsave("figure/FigS1F.Mean_gene_summary.by_celltype.pdf",p1,width = 12,height = 4)

# — 3 Fig2: Confirmation of neurogenic lineage and dissecting of NSC molecular heterogeneity in the postnatal human hippocampus -----------------------------------------------------------
# —— 3.1 Fig2A Cross-species comparison umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

load(file = paste("../Interation1234/figure/MMU_Pig_Mouse_Human.5types.rawObj.seurat.combined.RData",sep = ""))

subObj<-subset(seurat.combined,subset=Group%in%c("Human","Mouse_C","Pig","Rhesus"))
subObj$Group<-factor(subObj$Group,levels = c("Human","Mouse_C","Pig","Rhesus"))
levels(subObj$Group) %<>% {.[.%in%"Mouse_C"]<-"Mouse";.}
cell<-c("Astro","Astro-adult","Immature-Astro","Astro-juv","Astrocytes-qNSC",
        "RGL","RGL_young","Primed NSC",
        "nIPC","nIPC-perin",
        "NB","Neuroblast",
        "Granule cells","GC","Immature-GC","GC-juv","GC-adult",
        "aNSC")
subObj$Celltype<-factor(subObj$Celltype,levels = cell)
color18<-c("#FF8585","#FF4646","#8C0000","#BD2000","#DE4463",
           "#B8B5FF","#7868E6","#6155A6",
           "#A7D129","#C6E377",
           "#F1C550","#F9ED69",
           "#39A2DB","#00AFC1","#E0FCFF","#0B409C","#00FFF0",
           "black")

(pd<-DimPlot(subObj, reduction = "umap", group.by = "Celltype",
             split.by = "Group",raster=T,pt.size = 1,
             label = T,label.size = 4,repel = T,
             cols = color18)+NoLegend())

ggsave(filename = "figure/Fig2A.Cross-species_comparison.umap.pdf",pd,width = 18,height = 5)



# —— 3.2 Fig2B Hochgerner_H.et_al.2018 Marker check -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

DefaultAssay(tmp_Sobj) <- "RNA"
Cell<-c("Astrocytes-qNSC","Primed NSC","aNSC","Neuroblast","Granule cells")
color6<-c("#40407a","#00A087B2","#DC0000B2","#9196f7","#4DBBD5B2")
Subobj<-subset(tmp_Sobj,seurat_cluster_cell_typer%in%Cell) %>% 
  {.$seurat_cluster_cell_typer<-factor(.$seurat_cluster_cell_typer,levels = Cell);.}
DefaultAssay(Subobj) <- "RNA"
table(Subobj$Celltype)

gene1<-c("HOPX","HES1","SOX9","GFAP","SOX2","EMX2","NOTCH2","SOX4","VIM","ID3","NES","SMC4","CDK4","MXD3","EOMES","TBR1","BHLHE22 ","NEUROD2","FXYD7","SOX11","STMN2","NEUROD6","NNAT")

gene<-intersect(gene1,rownames(Subobj))
(pp<-DotPlot(object = Subobj, features = gene, group.by  = "Celltype",scale = T,scale.by = "size",
             cols = c("#21aba5", "#e84a5f")) + RotatedAxis())
Cell<-c("AS/qNSC","pNSC","aNSC","NB","GC","IN")
result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>% 
  {.$Cluster<-factor(.$Cluster,levels = rev(Cell));.}
(p1<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(2,10))+
    # scale_color_gradientn(colours = c("grey90","black","#7FB5FF","#FFEA85","#F73F52"))+
    scale_color_gradient2(low="black",mid="#007965",high ="#ffcc29",
                          midpoint = 0)+
    # scale_color_manual(values=color16)+
    theme_classic()+
    theme(axis.text=element_text(size=15, color="black"),legend.position = "top",axis.title = element_blank()) + RotatedAxis())
ggsave("figure/Fig2B.Hochgerner_H.et_al.2018.Marker_check.pdf",p1,width = 10,height = 3.5)

# —— 3.3 Fig2C AS_qNSC subtype umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

data.filt<- subset(tmp_Sobj,subset = seurat_cluster_cell_typer %in% "Astrocytes-qNSC") %>%
  {.$seurat_cluster_cell_typer<-factor(.$seurat_cluster_cell_typer,levels = "Astrocytes-qNSC");.}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.filt.list <- SplitObject(data.filt, split.by = "orig.ident")
data.filt.list <- lapply(X = data.filt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

astrocytes.anchors <- FindIntegrationAnchors(object.list = data.filt.list , dims = 1:20)
astrocytes.combined <- IntegrateData(anchorset = astrocytes.anchors, dims = 1:20)
DefaultAssay(astrocytes.combined) <- "RNA"
astrocytes.combined <- ScaleData(astrocytes.combined,verbose = T)
astrocytes.combined <- CellCycleScoring(astrocytes.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(astrocytes.combined) <- "integrated"
astrocytes.combined <- ScaleData(astrocytes.combined,verbose = T)
# astrocytes.combined <- ScaleData(astrocytes.combined, vars.to.regress = "CC.Difference")
astrocytes.combined <- RunPCA(astrocytes.combined, npcs = 30, verbose = T)
# t-SNE and Clustering
astrocytes.combined <- RunUMAP(astrocytes.combined, seed.use = 42, reduction = "pca", dims = 1:20)
astrocytes.combined <- RunTSNE(astrocytes.combined, seed.use = 2, reduction = "pca", dims = 1:20)
astrocytes.combined <- FindNeighbors(astrocytes.combined, reduction = "pca", dims = 1:20)
astrocytes.combined <- FindClusters(astrocytes.combined, resolution = seq(0.1,1,0.1))

astrocytes.combined$Celltype<-astrocytes.combined$integrated_snn_res.0.1
levels(astrocytes.combined$Celltype) %<>% {.[.%in%"2"]<-"qNSC2";.}
levels(astrocytes.combined$Celltype) %<>% {.[.%in%"1"]<-"AS";.}
levels(astrocytes.combined$Celltype) %<>% {.[.%in%"0"]<-"qNSC1";.}
save(astrocytes.combined,file = "../result/astrocytes.combined.Celltype.RData")

load("../result/astrocytes.combined.Celltype.RData")

(p1<-DimPlot(astrocytes.combined,group.by = "Celltype",raster = T,pt.size = 1,
             cols = c("#795548","#40407a","#DC0000B2"),label = T,label.size = 6)&
  tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                   arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
  theme(panel.grid = element_blank(),plot.title = element_blank(),
        # legend.position = "none",
        axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
        axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold"))&
  NoLegend())

ggsave("figure/Fig2C.AS_qNSC_subtype.Umap.pdf",p1,width = 4.5,height = 4.5)


# —— 3.4 Fig2D AS_qNSC subtype top marker heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

load("../result/astrocytes.combined.Celltype.RData")

ident<-"Celltype"
Idents(astrocytes.combined)<-ident
DefaultAssay(astrocytes.combined) <- "RNA"
all.markers <- FindAllMarkers(astrocytes.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top3 <- all.markers %>% group_by(cluster) %>% do(head(.,n=10))
(p4<-DoHeatmap(astrocytes.combined, features = unique(top3$gene),angle = 0,hjust = 0.5,
               group.colors =c("#4682B4","#40407a","#DC0000B2"),raster = F)+
    scale_fill_gradientn(colors = colorRampPalette(c("#1a2a6c", "white", "#c21e20"))(10),na.value = "white"))

ggsave("figure/Fig2D.AS_qNSC_subtype.top_marker.heatmap.pdf",p4,width = 12,height = 8)

# —— 3.5 Fig2EF AS_qNSC subtype Addmodulescore -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))
load("../result/astrocytes.combined.Celltype.RData")

tmp_Sobj<-astrocytes.combined
tmp_Idents<-"Celltype"
Idents(tmp_Sobj)<-tmp_Idents
DefaultAssay(tmp_Sobj) <- "RNA"

all.markers <- read.delim(paste("data/Fig2EF.AS_qNSC.subtype.subtype.res0.1.findallmarker.xls",sep = ""), stringsAsFactors=FALSE)

ortholog <- read.delim("data/Homo_MFA_MMU_Pig_Mouse.Gid_Gname.xls", header=T,stringsAsFactors=FALSE,na.strings = "") %>% 
  .[!is.na(.$Homo_Gname),]
diff_gene<-all.markers 
Marker<-read.delim("data/Fig2EF.Gene_set.sheet3.xls") %>%
  .[.$Celltype%in%c("RGL_progenitor","astrocytes"),] %>% 
  left_join(.,ortholog[,c("Homo_Gname","Mouse_Gname")],by=c("Marker"="Mouse_Gname")) %>% 
  {.[is.na(.$Homo_Gname),]$Homo_Gname<-.[is.na(.$Homo_Gname),]$Marker;.} %>% 
  unique() %>% subset(.,select=-Marker) %>% {colnames(.)<-c("Celltype","Marker");.} %>% 
  .[.$Marker%in%diff_gene$gene,] %>% unique() %>% .[.$Marker!="GFAP",]

setdiff(Marker$Marker,rownames(tmp_Sobj))
# "Glast"  "CD49f"  "GLT1"   "Cx43"   "SPOT14"

lt<-split(Marker$Marker, Marker$Celltype)

tmp_Plist<-list()
for (i in 1:length(lt)) {
  # i<-1
  gene <- lt[i]
  tmp_Sobj2 <- AddModuleScore(
    object = tmp_Sobj,
    features = gene,
    assay = "RNA",
    name = names(lt[i]),
    ctrl = 100)
  Idents(tmp_Sobj2)<-tmp_Idents
  pf<-FeaturePlot(tmp_Sobj2, features = paste(names(lt[i]),"1",sep = ""),repel = T,pt.size = 0.1,
                  reduction = "umap",label =T,label.size = 6)+
    scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold"))
  
  
  pp<-VlnPlot(tmp_Sobj2, features = paste(names(lt[i]),"1",sep = ""),split.plot = F,
              cols = paletteer_d("khroma::soil"),
              pt.size = 0, combine = T)&
    coord_flip()&
    theme(plot.title = element_blank(),axis.title = element_blank(),
          legend.position = "none")&
    geom_boxplot(alpha=.2,outlier.shape = NA)
  my_comparisons=list(c("qNSC2","AS"),c("AS","qNSC1"))
  pv<-pp$data %>% {colnames(.)<-c("Cell","ident");.} %>% 
    ggviolin(., x='ident', y='Cell', fill = 'ident')+
    geom_boxplot(alpha=0,outlier.shape = NA) +
    # stat_compare_means(comparisons = my_comparisons,method = "t.test") + 
    ggpubr::stat_compare_means(comparisons = my_comparisons,
                               symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                symbols = c("****", "***", "**", "*", "ns")),
                               method = "wilcox.test",
                               label = "p.signif",
                               size = 5)+
    scale_fill_manual(values = c("#795548","#40407a","#DC0000B2"))+
    theme(axis.title = element_blank())+
    coord_flip()+NoLegend()
  
  p<-((pf|pv)+plot_layout(widths = c(2, 1))) %>% as.ggplot()
  p<-p+labs(title = names(lt[i]))
  tmp_Plist[[i]]<-p
}
p<-(tmp_Plist[[2]]+tmp_Plist[[1]])+plot_layout(ncol = 2)

ggsave("figure/Fig2EF.AS_qNSC_subtype.Addmodulescore.png",p,width = 16,height = 5)
ggsave("figure/Fig2EF.AS_qNSC_subtype.Addmodulescore.pdf",p,width = 16,height = 5)


# —— 3.6 Fig2G AS_qNSC subtype Featureplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))
load("../result/astrocytes.combined.Celltype.RData")

DefaultAssay(astrocytes.combined) <- "RNA"

Gene<-c("S100B","GFAP","LPAR1","HOPX","STMN1","PROX1","SIRT2","ST18")
(p<-FeaturePlot(astrocytes.combined, features = Gene,repel = T,pt.size = 1,ncol = 4,raster = T,
                reduction = "umap",label = F)&
    scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")))&
    scale_y_continuous(breaks=NULL)&
    scale_x_continuous(breaks=NULL)&
    theme(legend.key.size = unit(0.8,'cm'),legend.text = element_blank(),
          # title = element_text(size = 0),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
          axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())&
    NoLegend())

ggsave("figure/Fig2G.AS_qNSC_subtype.Featureplot.pdf",p,width = 14,height = 8)
# —— 3.7 Fig2HI aNSC_pNSC top1000 GOBP -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_df <- read.delim("data/Fig2HI.aNSC.select.GO.xls") %>% {.$Group<-"aNSC";.}
p1<-tmp_df %>% 
  .[rev(order(.$pvalue)),] %>% 
  {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
  ggplot(., aes(x=Description,y=-log10(pvalue))) + 
  geom_bar(stat = "identity",fill="grey80") +
  geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))+
  # labs(title = "UP")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=13,colour = "black"),
        axis.text.y = element_blank(),axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),axis.line.y = element_blank(),
        legend.text= element_text(size=13),
        legend.title = element_text(size=13),
        strip.background = element_rect(color="#CF5051",fill = "#CF5051"),
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.wrap = unit(0.1, "mm"),
        plot.title = element_blank(),
        strip.text = element_text(size =16,face = "bold",colour = "white"))+
  facet_wrap(Group~.,strip.position = "left",dir="v")

tmp_df <- read.delim("data/Fig2HI.qNSC.select.GO.xls") %>% {.$Group<-"qNSC";.}
p2<-tmp_df %>% 
  .[rev(order(.$pvalue)),] %>% 
  {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
  ggplot(., aes(x=Description,y=-log10(pvalue))) + 
  geom_bar(stat = "identity",fill="grey80") +
  geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))+
  # labs(title = "UP")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=13,colour = "black"),
        axis.text.y = element_blank(),axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),axis.line.y = element_blank(),
        legend.text= element_text(size=13),
        legend.title = element_text(size=13),
        strip.background = element_rect(color="#CF5051",fill = "#CF5051"),
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.wrap = unit(0.1, "mm"),
        plot.title = element_blank(),
        strip.text = element_text(size =16,face = "bold",colour = "white"))+
  facet_wrap(Group~.,strip.position = "left",dir="v")

(p<-p2|p1)

ggsave("figure/Fig2HI.aNSC_pNSC.top1000_GOBP.barplot.pdf",p,width = 16,height = 5)

# —— 3.8 Fig2J Cell-cycle phases umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

color3<-c("#4682B4","#ffd66b","#ec4646")

cell<-c("qNSC1","qNSC2","pNSC","aNSC","NB")
tmp_obj <- subset(tmp_Sobj,subset=Minor2%in%cell)

(p1<-DimPlot(tmp_obj, reduction = "umap",cols = color3,label.size = 6,raster = T,
             label = F, pt.size = 1,group.by = "Phase",repel = T) + 
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

ggsave("figure/Fig2J.Cell-cycle_phases.umap.pdf",p1,width = 7,height = 6)

# — 4 FigS2 Related to Fig2: Distinguish qNSCs and astrocytes molecular heterogeneity in the postnatal human hippocampus. -----------------------------------------------------------
# —— 4.1 FigS2A AS_qNSC subtype Addmodulescore split by stage -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))
load("../result/astrocytes.combined.Celltype.RData")

astrocytes.combined$State<-astrocytes.combined$orig.ident

astrocytes.combined$State[which(astrocytes.combined$State %in% "4D")] <- "Neonatal"
astrocytes.combined$State[which(astrocytes.combined$State %in% c("31y","32y"))] <- "Middle-age"
astrocytes.combined$State[which(astrocytes.combined$State %in% c("50y","56y","60y","64f","64m","68y"))] <- "Aged"
astrocytes.combined$State[which(astrocytes.combined$State %in% "48y")] <- "Injury"

tmp_Plist2<-list()
for (tmp_age in c("Neonatal","Middle-age","Aged","Injury")) {
  tmp_Sobj<-subset(astrocytes.combined,subset=State%in%tmp_age)
  tmp_Idents<-"Celltype"
  Idents(tmp_Sobj)<-tmp_Idents
  DefaultAssay(tmp_Sobj) <- "RNA"
  
  all.markers <- read.delim(paste("data/Fig2EF.AS_qNSC.subtype.subtype.res0.1.findallmarker.xls",sep = ""), stringsAsFactors=FALSE)
  
  ortholog <- read.delim("data/Homo_MFA_MMU_Pig_Mouse.Gid_Gname.xls", header=T,stringsAsFactors=FALSE,na.strings = "") %>% 
    .[!is.na(.$Homo_Gname),]
  diff_gene<-all.markers 
  Marker<-read.delim("data/Fig2EF.Gene_set.sheet3.xls") %>%
    .[.$Celltype%in%c("RGL_progenitor","astrocytes"),] %>% 
    left_join(.,ortholog[,c("Homo_Gname","Mouse_Gname")],by=c("Marker"="Mouse_Gname")) %>% 
    {.[is.na(.$Homo_Gname),]$Homo_Gname<-.[is.na(.$Homo_Gname),]$Marker;.} %>% 
    unique() %>% subset(.,select=-Marker) %>% {colnames(.)<-c("Celltype","Marker");.} %>% 
    .[.$Marker%in%diff_gene$gene,] %>% unique() %>% .[.$Marker!="GFAP",] %>% 
    {.[.$Celltype%in%"RGL_progenitor",]$Celltype<-"RGL gene set";.} %>% 
    {.[.$Celltype%in%"astrocytes",]$Celltype<-"Astrocyte gene set";.}
  
  setdiff(Marker$Marker,rownames(tmp_Sobj))
  # "Glast"  "CD49f"  "GLT1"   "Cx43"   "SPOT14"
  
  lt<-split(Marker$Marker, Marker$Celltype)
  
  tmp_Plist<-list()
  for (i in 1:length(lt)) {
    # i<-1
    gene <- lt[i]
    tmp_Sobj2 <- AddModuleScore(
      object = tmp_Sobj,
      features = gene,
      assay = "RNA",
      name = names(lt[i]),
      ctrl = 100)
    Idents(tmp_Sobj2)<-tmp_Idents
    pf<-FeaturePlot(tmp_Sobj2, features = paste(gsub(" ",".",names(lt[i])),"1",sep = ""),repel = T,pt.size = 3,raster = T,
                    reduction = "umap",label =T,label.size = 4)+
      scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
      tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                       arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
      theme(panel.grid = element_blank(),plot.title = element_blank(),
            # legend.position = "none",
            axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
            axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold"))&
      NoLegend()
    
    
    pp<-VlnPlot(tmp_Sobj2, features = paste(gsub(" ",".",names(lt[i])),"1",sep = ""),split.plot = F,
                cols = paletteer_d("khroma::soil"),
                pt.size = 0, combine = T)&
      coord_flip()&
      theme(plot.title = element_blank(),axis.title = element_blank(),
            legend.position = "none")&
      geom_boxplot(alpha=.2,outlier.shape = NA)
    my_comparisons=list(c("qNSC2","AS"),c("AS","qNSC1"))
    pv<-pp$data %>% {colnames(.)<-c("Cell","ident");.} %>% 
      ggviolin(., x='ident', y='Cell', fill = 'ident')+
      geom_boxplot(alpha=0,outlier.shape = NA) +
      # stat_compare_means(comparisons = my_comparisons,method = "t.test") + 
      ggpubr::stat_compare_means(comparisons = my_comparisons,
                                 symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                  symbols = c("****", "***", "**", "*", "ns")),
                                 method = "wilcox.test",
                                 label = "p.signif",
                                 size = 5)+
      scale_fill_manual(values = c("#795548","#40407a","#DC0000B2"))+
      theme(axis.title = element_blank())+
      coord_flip()+NoLegend()
    
    p<-((pf|pv)+plot_layout(widths = c(1.5, 1))) %>% as.ggplot()
    p<-p+labs(title = paste(tmp_age,": ",names(lt[i]),sep = ""))+
      theme(panel.border = element_rect(fill=NA,color="black", size=0.75, linetype="solid"),
            plot.title = element_text(hjust=0.5,size = 16))
    tmp_Plist[[i]]<-p
  }
  p<-(tmp_Plist[[2]]+tmp_Plist[[1]])+plot_layout(ncol = 1)
  tmp_Plist2[[tmp_age]]<-as.ggplot(p)
}

ppp<-Reduce("+",tmp_Plist2)+plot_layout(ncol = 4)

ggsave(paste("figure/FigS2A.AS_qNSC_subtype.Addmodulescore.split_by_stage.pdf",sep = ""),ppp,width = 18,height = 16)

# —— 4.2 FigS2B AS_qNSC subtype cellpercent split by stage -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_cluster<-c("qNSC1","qNSC2","Astrocytes","M-AS")
tmp_final<-subset(tmp_Sobj,subset=Minor3%in%tmp_cluster) %>% 
  {.$Minor3<-factor(.$Minor3,levels = tmp_cluster);.}
DimPlot(tmp_final,group.by = "Minor3",label = T)


DefaultAssay(tmp_final) <- "RNA"
cluster_split<-split(tmp_final$Minor3,tmp_final$orig.ident)
cluster_split<-lapply(X = cluster_split, FUN = function(x) {
  as.data.frame(prop.table(table(x)))
})
df<-bind_rows(cluster_split, .id = "tmp_final$orig.ident")
colnames(df)<-c("Sample","Cluster","Percent")
df <- df[order(df$Cluster),]
age<-c("4D","31y","32y","48y","50y","56y","60y","64f","64m","68y")
df$Sample<-factor(df$Sample,levels = age)
subdf <- df %>% {.$Group<-as.character(.$Sample);.} %>% 
  {.$Group[which(.$Sample %in% "4D")] <- "Neonatal";.} %>%
  {.$Group[which(.$Sample %in% c("31y","32y"))] <- "Middle-age";.} %>%
  {.$Group[which(.$Sample %in% c("50y","56y","60y","64f","64m","68y"))] <- "Aged";.} %>%
  {.$Group[which(.$Sample %in% "48y")] <- "Injury";.} %>%
  .[.$Group%in%c("Neonatal","Middle-age","Aged","Injury"),] %>%
  .[.$Cluster%in%tmp_cluster,] %>% 
  tidyr::unite(., "Group2",Cluster, Group,remove=F) 

subdf<-subdf %>% group_by(Group2) %>% dplyr::summarise(mean = mean(Percent),sd=sd(Percent)) %>%
  left_join(subdf,.,by=c("Group2"="Group2")) %>%
  .[,c("Group","Cluster","mean","sd")] %>% unique() %>%
  {.$Group<-factor(.$Group,levels = c("Neonatal","Middle-age","Aged","Injury"));.} %>%
  {.$Cluster<-factor(.$Cluster,levels = tmp_cluster);.}

levels(subdf$Group) %<>% {.[.%in%"Neonatal"]<-"N";.}
levels(subdf$Group) %<>% {.[.%in%"Middle-age"]<-"M";.}
levels(subdf$Group) %<>% {.[.%in%"Aged"]<-"A";.}
levels(subdf$Group) %<>% {.[.%in%"Injury"]<-"I";.}

# subdf[is.na(subdf$sd),]$sd<-0

(p <- ggplot(data = subdf, mapping = aes(x = Cluster, y = mean,fill=Group)) + 
    geom_bar(stat= 'identity', position = 'dodge',color="black") +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    # scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")) +
    scale_fill_manual(values=c("#f5d662","#EB5230","#800026","#2A4E9D")) +
    ylab(label = "Cell Number percent") + xlab(label = "")+
    scale_x_discrete(expand = c(0, 0.6)) + 
    scale_y_continuous(expand = c(0, 0.002),limits=c(0, 0.6)) + 
    # scale_fill_npg()+
    #ylim(0,0.35) + 
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          axis.text.x = element_text(angle = 0,hjust = 0.5),
          legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.grid = unit(1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")) +
    facet_grid(~Group)+NoLegend()+RotatedAxis())
ggsave("figure/FigS2B.AS_qNSC_subtype.cellpercent.split_by_stage.pdf",p,width = 8,height = 4)

# —— 4.3 FigS2C qNSC2_OLG corr heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_sub<-subset(tmp_Sobj,subset=Minor2%in%c("qNSC1","qNSC2","Astrocytes","pNSC","aNSC","NB","OLG","OPC"))

cluster.averages <- AverageExpression(tmp_sub,assays = "integrated",group.by="Minor2")
ht_df<-cluster.averages[["integrated"]] %>% cor()
(p<-pheatmap::pheatmap(ht_df[c("Astrocytes","qNSC1","qNSC2","pNSC","aNSC","NB","OPC","OLG"),
                             c("Astrocytes","qNSC1","qNSC2","pNSC","aNSC","NB","OPC","OLG")],
                       show_rownames = T,show_colnames = T,cluster_rows = F,cluster_cols = F,use_raster=F,border_color="white",
                       display_numbers = TRUE,number_color = "black",fontsize_number = 15,fontsize=15,legend = F,
                       color = colorRampPalette(c("#72bcd5", "#F7F5F2", "#FD5D5D"))(30)))

ggsave("figure/FigS2C.qNSC2_OLG.cor_heatmap.pdf",p,width = 6,height = 6)
# —— 4.4 FigS2D pNSC vs aNSC DEGs heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

DefaultAssay(tmp_Sobj) <- "RNA"
Idents(tmp_Sobj)<-"Minor3"

tmp_Sub<-subset(tmp_Sobj,subset=Minor3%in%c("aNSC","pNSC"))
all.markers <- FindMarkers(tmp_Sobj, ident.1 = "aNSC", ident.2 = "pNSC", verbose = T) %>% 
  {.$gene<-rownames(.);.} %>% .[.$p_val<0.05,] %>% 
  {.$Diff<-"UP";.$Diff[.$avg_log2FC<0]<-"DOWN";.} %>% 
  .[grep("^MT-",.$gene,invert = T),]
tmp_G<-all.markers %>% .[order(.$Diff,abs(.$avg_log2FC),decreasing = T),] %>% group_by(Diff) %>% 
  do(head(.,20)) %>% .[order(.$Diff,abs(.$avg_log2FC),decreasing = T),]
(p<-DoHeatmap(tmp_Sub, features = unique(tmp_G$gene),raster = F)+
    scale_fill_gradientn(colors = colorRampPalette(c("#1a2a6c", "white", "#c21e20"))(20),na.value = "white"))
ggsave(paste("figure/FigS2D.pNSC_vs_aNSC.DEGs_heatmap.pdf",sep = ""),p,width = 10,height = 6)

# —— 4.5 FigS2E pNSC vs aNSC DEGs GOBP barplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")

rm(list=setdiff(ls(),"tmp_Sobj"))
tmp_df <- read.delim("data/FigS2E.aNSC_vs_pNSC.DEGs_GOBP.xls")

(p1<-tmp_df[tmp_df$Group%in%"UP",] %>% 
    .[rev(order(.$pvalue)),] %>% 
    {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=5,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_text(size=12,colour = "black"),
          axis.text = element_text(size=12,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=12),
          legend.title = element_text(size=12),
          strip.background = element_rect(color="#CF5051",fill = "#CF5051"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v"))
(p2<-tmp_df[tmp_df$Group%in%"DOWN",] %>% 
    .[rev(order(.$pvalue)),] %>% 
    {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=5,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_text(size=12,colour = "black"),
          axis.text = element_text(size=12,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=12),
          legend.title = element_text(size=12),
          strip.background = element_rect(color="#4C598C",fill = "#4C598C"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v"))

(pp<-p1/p2)

ggsave(paste("figure/FigS2E.pNSC_vs_aNSC.DEGs_GOBP.barplot.pdf",sep = ""),pp,width = 5,height = 6)

# — 5 Fig3: Discovery of novel markers distinguishing various types of NSCs and NBs in the human hippocampus -----------------------------------------------------------
# —— 5.1 Fig3AB find new marker by scHPF & Findallmarker -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

subseurat<-subset(tmp_Sobj,subset=orig.ident%in%"4D"&Minor2%in%c("qNSC1","qNSC2","pNSC","aNSC","NB"))
subseurat$seurat_cluster_cell_typer<-subseurat$Minor2
subseurat$seurat_cluster_cell_typer<-factor(subseurat$seurat_cluster_cell_typer,levels = c("qNSC1","qNSC2","pNSC","aNSC","NB"))


DefaultAssay(subseurat) <- "RNA"
Idents(subseurat)<-"seurat_cluster_cell_typer"

all.markers <- read.delim("data/Fig3A.findallmarkers.all.markers.v1.xls", stringsAsFactors=FALSE)

top20 <- all.markers[all.markers$cluster%in%c("Primed NSC","aNSC","Neuroblast"),] %>% .[.$p_val<0.05,] %>% 
  {.$cluster<-factor(.$cluster,levels = c("Primed NSC","aNSC","Neuroblast"));.} %>% 
  .[order(.$cluster),] %>% 
  group_by(cluster) %>% do(head(.,n=15))

tmp_gene<-c("LRRC3B","GNA14","RHOJ","SLC4A4","GLI3","CNTNAP4","ELMO1","IL1RAPL1","ST18",top20$gene) %>% 
  setdiff(.,c("HNRNPA2B1","AC124312.1","CEP170","CEP170","HSPE1","HNRNPH1","SRSF11","MYT1L"))
pp<-DotPlot(object = subseurat, features = tmp_gene, group.by  = "seurat_cluster_cell_typer",scale = T,scale.by = "size",
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()

result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
# result$mean_Expr[result$mean_Expr>1]=1

(p<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(2,10))+
    scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    theme_classic()+
    theme(axis.text=element_text(size=15, color="black"),legend.position = "top") + RotatedAxis())

ggsave(filename = "figure/Fig3A.find_new_marker.by_findallmarkers.bubble.pdf",p,width = 20,height = 4)

top=100
rank <- read.delim("data/Fig3A.sample4D_4types_K7.ranked_genes.txt", header=FALSE) %>% 
  {colnames(.)<-paste("V",seq(2,ncol(.)+1),sep = "");.} %>% .[1:top,] %>% rownames_to_column() %>% 
  gather(colname,gene,-rowname) %>% 
  {.$celltype<-.$colname;.} %>% 
  {.$celltype[.$celltype%in%c("V2")]<-"pNSC-2";.} %>% 
  {.$celltype[.$celltype%in%c("V3")]<-"AS-qNSC-2";.} %>%
  {.$celltype[.$celltype%in%c("V4")]<-"pNSC-1";.} %>%
  {.$celltype[.$celltype%in%c("V5")]<-"AS-qNSC-1";.} %>%
  {.$celltype[.$celltype%in%c("V6")]<-"AS-qNSC-3";.} %>%
  {.$celltype[.$celltype%in%c("V7")]<-"NB";.} %>%
  {.$celltype[.$celltype%in%c("V8")]<-"aNSC";.} %>% 
  {.$Type<-.$celltype;.} %>% 
  {.$Type[.$Type%in%c("AS-qNSC-1","AS-qNSC-2","AS-qNSC-3")]<-"AS-qNSC";.} %>%
  {.$Type[.$Type%in%c("pNSC-1","pNSC-2")]<-"pNSC";.} %>%
  {.$Type[.$Type%in%c("NB")]<-"NB";.} %>%
  {.$Type[.$Type%in%c("aNSC")]<-"aNSC";.} 


markers<-rank %>% group_by(celltype) %>% do(head(.,n=10)) %>% .[.$Type!="AS-qNSC",] %>% 
  {.$celltype<-factor(.$celltype,
                      levels = c("pNSC-1","pNSC-2","aNSC","NB"));.} %>% 
  .[order(.$celltype),]

tmp_gene<-c("LRRC3B","GNA14","RHOJ","SLC4A4","GLI3","CNTNAP4","ELMO1","IL1RAPL1","ST18",markers$gene) %>% 
  setdiff(.,c("VIT","TFAP2C","NPFFR","FMO1","CPXM1","GDF15","APLNR","IRX1","AGBL1","KCNJ16",
              "IL17B","AC092111.3","OR5AU1","CWH43","C9orf57","TMEM236","KYNU","LIF"))
pp<-DotPlot(object = subseurat,
            features = tmp_gene,
            group.by  = "seurat_cluster_cell_typer",scale = T,scale.by = "size",
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()

result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
# result$mean_Expr[result$mean_Expr>1]=1

(p<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(2,10))+
    scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    # scale_color_manual(values=color16)+
    theme_classic()+
    theme(axis.text=element_text(size=15, color="black"),legend.position = "top") + RotatedAxis())

ggsave(filename = "figure/Fig3A.find_new_marker.by_scHPF.bubble.pdf",p,width = 20,height = 4)



# —— 5.3 Fig3C New marker Featureplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

subseurat<-subset(tmp_Sobj,subset = orig.ident %in% c("4D") & Major %in% setdiff(levels(tmp_Sobj$Major),c("UN1","UN2"))) %>%
  {.$orig.ident<-factor(.$orig.ident,levels = c("4D"));.}
DefaultAssay(subseurat) <- "RNA"
Idents(subseurat)<-"Minor"

tmp_G <- read.delim("data/Fig3C.gene.for_feature.v6.list")

# devtools::install_github('junjunlab/scRNAtoolVis')
# 
# library(scRNAtoolVis)

library(tidydr)
tmp_Sobj<-subseurat
(p<-FeaturePlot(tmp_Sobj, features = tmp_G$Gene,pt.size = 1,raster = T,
               ncol = 4,
               reduction = "umap",label = F,label.size = 5)&
  scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")))&
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
  theme(panel.grid = element_blank(),legend.position = "none",
        axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5,size = 15,face = "bold")))

ggsave(filename = "figure/Fig3C.New_marker.Featureplot.pdf",p,width = 12,height = 8)

# —— 5.4 Fig3DF NB_GC markers heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

cell<-c("AS/qNSC","pNSC","aNSC","NB","GC","GABA-IN")
subseurat<-subset(tmp_Sobj,subset = Major %in% cell & orig.ident %in% c("4D")) %>%
  {.$orig.ident<-factor(.$orig.ident,levels = c("4D"));
  .$Major<-factor(.$Major,levels = cell);.}
DefaultAssay(subseurat) <- "RNA"
tmp_Idents<-"Major"
Idents(subseurat)<-tmp_Idents

tmp_G<-c("DLX1","GABRG3","CCK","SLC6A11","SLC6A1","GAD1","GAD2","CNR1","GRM1","RELN","VIP",
         "CALB2","SEMA6D","TENM1","DCX","PROX1","GRIP2","GRIK1","MACROD1","GRIK3","CALB1","IL16",
         "NRGN","FXYD7","NRN1","GNG3","TCEAL5","RPS13","TMSB10","NNAT","SLC17A7","MARCKSL1","THY1","NEUROD2","AC005258.1","AC004949.1","GADD45G","AL627171.2","H2AFZ","HSPA1B","TERF2IP","CALM3","DNAJA1","TTC9B")

cluster.averages <- AverageExpression(subseurat,assays = "RNA",group.by=tmp_Idents)
(p<-pheatmap::pheatmap(t(cluster.averages[["RNA"]][tmp_G,]),scale = "column",border_color = "white",gaps_col = c(11,22,22,22,22,22,22,22,22),
                       show_rownames = T,show_colnames = T,cluster_rows = F,cluster_cols = F,use_raster=F,
                       # annotation_col = ann_col[,"Type1",drop=F],annotation_colors = ann_colors,
                       display_numbers = F,number_color = "black",fontsize_number = 15,fontsize=15,
                       color = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(30)))
ggsave(filename = "figure/Fig3DF.NB_GC.markers_heatmap.pdf",p,width = 16,height = 3)

# —— 5.5 Fig3E Percent & Expression pointplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

DefaultAssay(tmp_Sobj) <- "RNA"


tmp_GF<-read.delim("../data/NB_marker.Findallmarkers.xls") %>% .[.$p_val_adj<0.01,]
tmp_GS<-read.delim("../data/NB_marker.scHPF.xls")

intersect(tmp_GF$gene,tmp_GS$gene)

pp<-DotPlot(object = tmp_Sobj, features = unique(tmp_GF$gene), group.by  = "Celltype",scale = T)+
  RotatedAxis()
result_GF<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>% 
  .[.$Cluster%in%"IN",]


pp<-DotPlot(object = tmp_Sobj, features = unique(tmp_GS$gene), group.by  = "Celltype",scale = T)+
  RotatedAxis()
result_GS<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>% 
  .[.$Cluster%in%"IN",]


data<-rbind(result_GS[,c("Percent Expressed","Average Expression")],
            result_GF[,c("Percent Expressed","Average Expression")]
) %>% 
  {.$Group<-"<50%";.} %>% 
  {.[.$`Percent Expressed`>50,]$Group<-">50%";.}
p<-ggplot(data, aes(x=`Percent Expressed`, y=`Average Expression`,color=`Percent Expressed`)) +
  geom_point()+
  scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")))+
  theme_bw()+
  theme(legend.position="none",
        panel.background = element_blank(),panel.grid = element_blank())
pp<-ggMarginal(p, type="histogram")

ggsave(filename = "figure/Fig3E.Percent_Expression.pointplot.pdf",pp,width = 6,height = 6)

# — 6 FigS3 Related to Fig3: Reported neuroblast genes were widely distributed in the adult human interneurons. -----------------------------------------------------------
# —— 6.1 FigS3A corr with ref heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

DefaultAssay(SeuratObj) <- "RNA"
table(SeuratObj$Major)


sample<-"Homo_MMU.filterObj"
load(file = paste("../Interation1234/figure/",sample,".seurat.combined.RData",sep = ""))
cell_type<-c("Granule cells","Interneuron","Neuroblast","aNSC","Primed NSC",
             "Astrocytes-qNSC","Mature astrocytes","OPC","Oligodendrocytes","Microglia",
             "Endothelial cells","Pyramidal neuron","CR","Pericytes","Unknown1","Unknown2")
seurat.combined$Celltype<-factor(seurat.combined$Celltype,levels = c(cell_type,sort(setdiff(unique(seurat.combined$Celltype),cell_type))))
levels(seurat.combined$Celltype) %<>% {.[.%in%"Mature astrocytes"]<-"M-AS";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Astrocytes-qNSC"]<-"AS/qNSC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Primed NSC"]<-"pNSC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Neuroblast"]<-"NB";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Granule cells"]<-"GC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Interneuron"]<-"GABA-IN";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Oligodendrocytes"]<-"OLG";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Microglia"]<-"MG";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Endothelial cells"]<-"EC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Pyramidal neuron"]<-"Pyr";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Pericytes"]<-"Per";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Unknown1"]<-"UN1";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Unknown2"]<-"UN2";.}
seurat.combined$Celltype<-as.character(seurat.combined$Celltype)
seurat.combined@meta.data %<>% {.[.$Species%in%"Rhesus"&.$Celltype%in%"NB",]$Celltype<-"ref_NB";.}
seurat.combined@meta.data %<>% {.[.$Species%in%"Rhesus"&.$Celltype%in%"GC",]$Celltype<-"ref_GC";.}
seurat.combined@meta.data %<>% {.[.$Species%in%"Rhesus"&.$Celltype%in%"OPC",]$Celltype<-"ref_OPC";.}
seurat.combined@meta.data %<>% {.[.$Species%in%"Rhesus"&.$Celltype%in%"CR",]$Celltype<-"ref_CR";.}
table(seurat.combined$Species,seurat.combined$Celltype)
cluster.averages <- AverageExpression(seurat.combined,assays = "integrated",group.by="Celltype")
own_cell<-c("AS/qNSC","pNSC","aNSC","NB","GC","GABA-IN")
ref_cell<-c("Astro","RGL","nIPC","ref_NB","ref_GC","InN")
ht_df<-cluster.averages[["integrated"]] %>% cor()
(p1<-pheatmap::pheatmap(ht_df[own_cell,ref_cell],
                        show_rownames = T,show_colnames = T,cluster_rows = F,cluster_cols = F,use_raster=F,border_color="white",
                        display_numbers = TRUE,number_color = "black",fontsize_number = 15,fontsize=15,legend = F,
                        color = colorRampPalette(c("#72bcd5", "#F7F5F2", "#FD5D5D"))(30)))


sample<-"Homo_MouseC.filterObj"
load(file = paste("../Interation1234/figure/",sample,".seurat.combined.RData",sep = ""))
cell_type<-c("Granule cells","Interneuron","Neuroblast","aNSC","Primed NSC",
             "Astrocytes-qNSC","Mature astrocytes","OPC","Oligodendrocytes","Microglia",
             "Endothelial cells","Pyramidal neuron","CR","Pericytes","Unknown1","Unknown2")
seurat.combined$Celltype<-factor(seurat.combined$Celltype,levels = c(cell_type,sort(setdiff(unique(seurat.combined$Celltype),cell_type))))
table(seurat.combined$Species,seurat.combined$Celltype)
levels(seurat.combined$Celltype) %<>% {.[.%in%"Mature astrocytes"]<-"M-AS";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Astrocytes-qNSC"]<-"AS/qNSC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Primed NSC"]<-"pNSC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Neuroblast"]<-"NB";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Granule cells"]<-"GC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Interneuron"]<-"GABA-IN";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Oligodendrocytes"]<-"OLG";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Microglia"]<-"MG";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Endothelial cells"]<-"EC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Pyramidal neuron"]<-"Pyr";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Pericytes"]<-"Per";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Unknown1"]<-"UN1";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Unknown2"]<-"UN2";.}
seurat.combined$Celltype<-as.character(seurat.combined$Celltype)
seurat.combined@meta.data %<>% {.[.$Species%in%"Mouse"&.$Celltype%in%"NB",]$Celltype<-"ref_NB";.}
seurat.combined@meta.data %<>% {.[.$Species%in%"Mouse"&.$Celltype%in%"OPC",]$Celltype<-"ref_OPC";.}
table(seurat.combined$Species,seurat.combined$Celltype)
cluster.averages <- AverageExpression(seurat.combined,assays = "integrated",group.by="Celltype")
own_cell<-c("AS/qNSC","pNSC","aNSC","NB","GC","GABA-IN")
ref_cell<-c("Astro-juv","Astro-adult","RGL","RGL_young","ref_NB","Immature-GC",
            "GC-juv","GC-adult","Immature-GABA","GABA")
ht_df<-cluster.averages[["integrated"]] %>% cor()
(p2<-pheatmap::pheatmap(ht_df[own_cell,ref_cell],
                        show_rownames = T,show_colnames = T,cluster_rows = F,cluster_cols = F,use_raster=F,border_color="white",
                        display_numbers = TRUE,number_color = "black",fontsize_number = 15,fontsize=15,legend = F,
                        color = colorRampPalette(c("#72bcd5", "#F7F5F2", "#FD5D5D"))(30)))

All<-"refHuman_Human"
load(paste("../Interation1234/figure/",All,".filterObj.seurat.combined.RData",sep = ""))
cell_type<-c("Granule cells","Interneuron","Neuroblast","aNSC","Primed NSC",
             "Astrocytes-qNSC","Mature astrocytes","OPC","Oligodendrocytes","Microglia",
             "Endothelial cells","Pyramidal neuron","CR","Pericytes","Unknown1","Unknown2")
seurat.combined$Celltype<-factor(seurat.combined$Celltype,levels = c(cell_type,sort(setdiff(unique(seurat.combined$Celltype),cell_type))))
table(seurat.combined$Species,seurat.combined$Celltype)
levels(seurat.combined$Celltype) %<>% {.[.%in%"Mature astrocytes"]<-"M-AS";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Astrocytes-qNSC"]<-"AS/qNSC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Primed NSC"]<-"pNSC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Neuroblast"]<-"NB";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Granule cells"]<-"GC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Interneuron"]<-"GABA-IN";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Oligodendrocytes"]<-"OLG";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Microglia"]<-"MG";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Endothelial cells"]<-"EC";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Pyramidal neuron"]<-"Pyr";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Pericytes"]<-"Per";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Unknown1"]<-"UN1";.}
levels(seurat.combined$Celltype) %<>% {.[.%in%"Unknown2"]<-"UN2";.}
seurat.combined$Celltype<-as.character(seurat.combined$Celltype)
table(seurat.combined$Species,seurat.combined$Celltype)
cluster.averages <- AverageExpression(seurat.combined,assays = "integrated",group.by="Celltype")
own_cell<-c("AS/qNSC","GC","GABA-IN")
ref_cell<-c("Astro","DG GC","InN")
ht_df<-cluster.averages[["integrated"]] %>% cor()
(p3<-pheatmap::pheatmap(ht_df[own_cell,ref_cell],
                        show_rownames = T,show_colnames = T,cluster_rows = F,cluster_cols = F,use_raster=F,border_color="white",
                        display_numbers = TRUE,number_color = "black",fontsize_number = 15,fontsize=15,legend = F,
                        color = colorRampPalette(c("#72bcd5", "#F7F5F2", "#FD5D5D"))(30)))

p<-as.ggplot(p1)|as.ggplot(p2)|as.ggplot(p3)
ggsave("figure/FigS3A.corr_with_ref.heatmap.pdf",p,width = 18,height = 4)

# —— 6.2 FigS3BC NB GABA-IN specific bubble split by age -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

DefaultAssay(tmp_Sobj) <- "RNA"
table(tmp_Sobj$Celltype)
levels(tmp_Sobj$Major) %<>% {.[. %in% c("IN")] <- "GABA-IN";.}

celltype<-c("NB","GC","GABA-IN")
age<-c("31y","32y","50y","56y","60y","64f","64m","68y")

gene<-c("NRGN","FXYD7","NRN1","GNG3","TCEAL5","RPS13","TMSB10","NNAT","SLC17A7","MARCKSL1","THY1")

for (i in 1:length(age)) {
  # i<-1
  sub<-subset(tmp_Sobj,subset = Major %in% celltype & orig.ident %in% age[i]) %>% 
    {.$Major<-factor(.$Major,levels=celltype);.} %>% 
    {.$orig.ident<-factor(.$orig.ident,levels=age[i]);.}
  DefaultAssay(sub) <- "RNA"
  
  if (i==1 | i==5) {
    pp<-DotPlot(object = sub, features = unique(gene), group.by  = "Major",scale = T,scale.by = "size",col.min = -1.5,col.max = 2.5,
                cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
    result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
    (p<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,9),limits = c(0,100))+
        scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
        labs(title = age[i])+
        theme_classic()+
        theme(axis.text=element_text(size=13, color="black"),legend.position = "top",
              axis.title.x = element_blank(),axis.title.y = element_blank()) + RotatedAxis())
    assign(paste("p",i,sep = ""),p)
  } else {
    pp<-DotPlot(object = sub, features = unique(gene), group.by  = "Major",scale = T,scale.by = "size",col.min = -1.5,col.max = 2.5,
                cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
    result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
    (p<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,9),limits = c(0,100))+
        scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
        # scale_color_manual(values=color16)+
        labs(title = age[i])+
        theme_classic()+
        theme(axis.text=element_text(size=13, color="black"),legend.position = "top",
              axis.title.x = element_blank(),axis.title.y = element_blank(),
              axis.line.y = element_blank(),axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) + RotatedAxis())
    assign(paste("p",i,sep = ""),p)
  }
}

p<-(p1+p2+p3+p4+p5+p6+p7+p8)+plot_layout(ncol = 4,guides = "collect")&theme(legend.position = "top")
ggsave("figure/FigS3BC.NB_GABA-IN_specific.bubble.split_by_age.part1.pdf",p,width = 20,height = 6,dpi = 300)


gene<-c("NEUROD2","AC005258.1","AC004949.1","GADD45G","AL627171.2","H2AFZ","HSPA1B","TERF2IP","CALM3","DNAJA1","TTC9B")
for (i in 1:length(age)) {
  # i<-1
  sub<-subset(tmp_Sobj,subset = Major %in% celltype & orig.ident %in% age[i]) %>% 
    {.$Major<-factor(.$Major,levels=celltype);.} %>% 
    {.$orig.ident<-factor(.$orig.ident,levels=age[i]);.}
  DefaultAssay(sub) <- "RNA"
  
  if (i==1 | i==5) {
    pp<-DotPlot(object = sub, features = unique(gene), group.by  = "Major",scale = T,scale.by = "size",col.min = -1.5,col.max = 2.5,
                cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
    result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
    (p<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,9),limits = c(0,100))+
        scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
        labs(title = age[i])+
        theme_classic()+
        theme(axis.text=element_text(size=13, color="black"),legend.position = "top",
              axis.title.x = element_blank(),axis.title.y = element_blank()) + RotatedAxis())
    assign(paste("p",i,sep = ""),p)
  } else {
    pp<-DotPlot(object = sub, features = unique(gene), group.by  = "Major",scale = T,scale.by = "size",col.min = -1.5,col.max = 2.5,
                cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
    result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
    (p<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
        geom_point()+
        scale_size_continuous(range=c(1,9),limits = c(0,100))+
        scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
        # scale_color_manual(values=color16)+
        labs(title = age[i])+
        theme_classic()+
        theme(axis.text=element_text(size=13, color="black"),legend.position = "top",
              axis.title.x = element_blank(),axis.title.y = element_blank(),
              axis.line.y = element_blank(),axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) + RotatedAxis())
    assign(paste("p",i,sep = ""),p)
  }
}

p<-(p1+p2+p3+p4+p5+p6+p7+p8)+plot_layout(ncol = 4,guides = "collect")&theme(legend.position = "top")
ggsave("figure/FigS3BC.NB_GABA-IN_specific.bubble.split_by_age.part2.pdf",p,width = 20,height = 6,dpi = 300)

# — 7 FigS4 Related to Fig3: Neuroblast marker DCX were expressed in interneuron in macaque hippocampus of 3 months. -----------------------------------------------------------
# — 8 Fig4: The transcriptional dynamics predicated by RNA velocity and pseudotime reconstruction revealed developmental potentials of NSC in the neonatal human hippocampus -----------------------------------------------------------
# —— 8.1 Fig4A neonatal NSCs RNA velocity -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

library(Seurat)
library(devtools)
# install_github("velocyto-team/velocyto.R")
# Shell : git clone https://github.com/velocyto-team/velocyto.R
# devtools::install_local("/data1/zhur/soft/velocyto.R/",dep=T,upgrade_dependencies=F)
library(velocyto.R)
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_local("/home/devdata/zr/proj/scRNA/ltq/seruat_6.16/sample8_v4/velocyto/seurat-disk")
library(SeuratDisk)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratObject)
library(SeuratWrappers)
library(dplyr)
library(ggplot2)
library(loomR)
library(hdf5r)
# unloadNamespace("SeuratDisk")
# remotes::install_version("Seurat", version = "2.3.4")
# remotes::install_github("mojaveazure/seurat-disk")
# devtools::install_github(repo = "hhoeflin/hdf5r")
# devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")

DefaultAssay(tmp_Sobj) <- "RNA"

cell<-c("qNSC1","Astrocytes","qNSC2","pNSC","aNSC","NB","GC1","GC2")
color7<-c("#F9D56E","#1E212D","#FFBCBC","#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
seurat_object<-subset(tmp_Sobj, subset = orig.ident %in% "4D"& Minor2 %in% cell)

seurat_object$Minor<-seurat_object$Minor2
seurat_object$Minor<-factor(seurat_object$Minor,levels = cell)

table(seurat_object$orig.ident,seurat_object$Minor)

DimPlot(seurat_object, reduction = "umap",label.size = 6,cols = color7,
        label = T, pt.size = 0.4,group.by = "Minor",repel = T)


ldat <- ReadVelocity(file = "../4D_48y_velocyto/4DHLH_20200313NA/velocyto/4DHLH_20200313NA.loom")

emat <- ldat$spliced
nmat <- ldat$unspliced


emb <- seurat_object@reductions$umap@cell.embeddings
id<-rownames(emb)

s4<-seurat_object
# Estimate the cell-cell distances
cell.dist <- as.dist(1-armaCor(t(s4@reductions$umap@cell.embeddings)))
cell.dist <- s4@reductions$umap@cell.embeddings %>% t()
cell.dist <- 1-armaCor(cell.dist) %>% as.dist()

colnames(emat) <- paste("4D_",substring(colnames(emat),18,33),"-1",sep="")
colnames(nmat) <- paste("4D_",substring(colnames(nmat),18,33),"-1",sep="")
emat<-emat[,id]
nmat<-nmat[,id]
ncol(nmat)
nrow(emb)
nrow(nmat)
# nrow(cell.dist)

# I'm not sure what this parameter does to be honest. 0.02 default
# perform gamma fit on a top/bottom quantiles of expression magnitudes
fit.quantile <- 0.02
# Main velocity estimation
rvel.cd <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=20,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile,
                                            n.cores=20,
                                            verbose=T)

# This section gets the colors out of the seurat tSNE object so that my seurat and velocyto plots use the same color scheme.
Idents(seurat_object)<-'Minor'

gg <- UMAPPlot(seurat_object,cols=color7)
ggplot_build(gg)$data
colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(emb)

Vresult <- show.velocity.on.embedding.cor(emb,rvel.cd,n=100,scale='sqrt',
                                          cell.colors=ac(colors,alpha=0.8),
                                          cex=0.8,arrow.scale=2,show.grid.flow=T,
                                          min.grid.cell.mass=1.0,grid.n=50,arrow.lwd=1,
                                          do.par=F,cell.border.alpha = 0.1,
                                          # cc=p1$cc,
                                          return.details=T,
                                          n.cores=30,main="4D Cell Velocity")

umap_data <- emb %>%
  as.data.frame() %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  mutate(X0 = Vresult$arrows[, "x0"],
         X1 = Vresult$arrows[, "x1"],
         Y0 = Vresult$arrows[, "y0"],
         Y1 = Vresult$arrows[, "y1"]) %>%
  mutate(X2 = X0 + (X1 - X0) * 1,
         Y2 = Y0 + (Y1 - Y0) * 1) %>%
  mutate(Cluster = seurat_object@meta.data[rownames(emb),c("orig.ident","Minor")]$Minor)

ggplot(umap_data) +
  geom_point(aes(x = UMAP1, y = UMAP2, colour = Cluster)) +
  geom_segment(aes(x = X0, xend = X2, y = Y0, yend = Y2),
               arrow = arrow(length = unit(3, "points"), type = "closed"),
               colour = "grey20", alpha = 0.8) +
  theme_minimal()


umap_arrows <- Vresult$garrows %>%
  as.data.frame() %>%
  mutate(x2 = x0 + (x1 - x0) * 5,
         y2 = y0 + (y1 - y0) * 5)

(p<-ggplot(umap_data) +
    geom_point(aes(x = UMAP1, y = UMAP2, colour = Cluster),alpha=1,size=0.5) +
    geom_curve(data = umap_arrows,curvature = 0.1,
               aes(x = x0, xend = x2, y = y0, yend = y2),
               size = 0.5,
               arrow = arrow(length = unit(5, "points"), type = "closed"),
               colour = "grey30", alpha = 0.8) +
    scale_color_manual(values = color7)+
    guides(color = guide_legend(override.aes = list(size = 6)))+
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.title = element_blank(),
          axis.text = element_blank()))


ggsave("figure/Fig4A.neonatal_NSCs.RNA_velocity.pdf",p,width = 10,height = 9)
# —— 8.2 Fig4B qNSC1_qNSC2 vs pNSC DEGs enrichment -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

sub<-tmp_Sobj %>% subset(.,orig.ident %in% "4D")
DefaultAssay(sub) <- "RNA"
Idents(sub)<-"Subtype"

DimPlot(sub,group.by = "Subtype")

library( "clusterProfiler")
load(file = '/data1/zhur/proj/scRNA/enrichment_annotation/human_enrichment.RData')
gtf <- read.delim("data/cellranger.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
id2name<-gtf %>% .[.$V3=="gene",] %>% {.$V9<-gsub(".*gene_id ","",.$V9);.$V9<-gsub(";.*gene_name ",",",.$V9);.$V9<-gsub(";.*","",.$V9);.} %>%
  separate(data = ., col = V9, into = c("gene_id","gene_name"), sep = ",") %>% .[,c("gene_id","gene_name")] %>% unique()
minGSSize=3
maxGSSize=500

for (i in c("qNSC1","qNSC2")) {
  # i<-"qNSC1"
  diff<-FindMarkers(sub, ident.1 = i, ident.2 = "pNSC", verbose = T) %>% 
    {.$gene<-rownames(.);.} %>% .[.$p_val<0.05,] %>% 
    {.$Diff<-"UP";.$Diff[.$avg_log2FC<0]<-"DOWN";.}
  
  ## BP
  term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
  term2gene<-term2gene[!is.na(term2gene$gene_name),]
  term2name <- gobp[,c('GO','Name')]
  ego_BP <- enricher(gene = diff[diff$Diff%in%"UP",]$gene,pvalueCutoff = 0.05,pAdjustMethod = "fdr",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
  query_BP_result<-as.data.frame(ego_BP@result)
  
  ## BP
  term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
  term2gene<-term2gene[!is.na(term2gene$gene_name),]
  term2name <- gobp[,c('GO','Name')]
  ego_BP <- enricher(gene = diff[diff$Diff%in%"DOWN",]$gene,pvalueCutoff = 0.05,pAdjustMethod = "fdr",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
  query_BP_result<-as.data.frame(ego_BP@result)

}

tmp_df <- read.delim("data/Fig4B.qNSC1_qNSC2_vs_pNSC.GOBP.xls")

(p1<-tmp_df[tmp_df$Group%in%"qNSC1",] %>% 
    .[rev(order(.$pvalue)),] %>% 
    {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_text(size=13,colour = "black"),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#4C598C",fill = "#4C598C"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v"))
(p2<-tmp_df[tmp_df$Group%in%"qNSC2",] %>% 
    .[rev(order(.$pvalue)),] %>% 
    {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_text(size=13,colour = "black"),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#4C598C",fill = "#4C598C"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v"))

p<-p1|p2
ggsave("figure/Fig4B.qNSC1_qNSC2_vs_pNSC.DEGs_enrichment.bar_plot.pdf",p,width = 13.5,height = 4.5,dpi = 500)

# —— 8.3 Fig4CD neonatal pseudotime trajectory & gene expression -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

cell<-c("pNSC","aNSC","NB","GC1","GC2")
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
scRNAsub<-subset(tmp_Sobj, subset = Minor %in% cell & orig.ident %in% "4D")
scRNAsub$Minor <- factor(scRNAsub$Minor,levels = cell)
scRNAsub$orig.ident <- factor(scRNAsub$orig.ident,levels = "4D")
scRNAsub$Minor %>% table()
# pNSC   aNSC     NB    GC1    GC2 
# 1470    944   1045   3592   1852

color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(pp<-DimPlot(scRNAsub,group.by = "Minor",cols = color7))


data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(as.matrix(data),
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=50, relative_expr = TRUE)

disp_table <- dispersionTable(mycds)

disp.genes <- rownames(scRNAsub@assays$integrated)
mycds <- setOrderingFilter(mycds, disp.genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
# Error in graph.adjacency.dense(adjmatrix, mode = mode, weighted = weighted,  : 
#                                  long vectors not supported yet: ../../src/include/Rinlinedfuns.h:522
# ncells > 50k
# https://github.com/cole-trapnell-lab/monocle-release/issues/138
mycds <- orderCells(mycds)
mycds$type<-factor(mycds$Minor,levels = cell)
disp.genes<-as.character(disp.genes)
filename<-paste("result/4D.noroot.expr.Doublet_remove.v6.RData",sep = "")
save(mycds,disp.genes,file = filename)

GM_state <- function(mycds){
  if (length(unique(pData(mycds)$State)) > 1){
    T0_counts <- table(pData(mycds)$State, pData(mycds)$type)[,"pNSC"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
mycds <- orderCells(mycds, root_state = GM_state(mycds))
plot_cell_trajectory(mycds, color_by = "Pseudotime")

filename<-paste("result/4D.root-pNSC.expr.Doublet_remove.v6.RData",sep = "")
save(disp.genes,mycds,file = filename)
list<-load(filename)
mycds@phenoData@data %<>% {.[.$State%in%"1"&.$Celltype%in%"GC",]$Minor<-"GC2";.} %<>% 
  {.[.$State%in%"2"&.$Celltype%in%"GC",]$Minor<-"GC1";.}

plot_cell_trajectory(mycds, color_by = "Minor",cell_size = 1)

(p1<-plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 1)+
    scale_color_gradientn(colors=c("#916BBF","#54BAB9","#2EB086","#A7D129"),name = "Pseudotime"))
(p2<-plot_cell_trajectory(mycds, color_by = "Minor",cell_size = 1)+
    scale_color_manual(values = color7))
(p3<-pp$data %>% cbind(.,mycds$Pseudotime) %>% {colnames(.)<-c("UMAP_1","UMAP_2","Celltype","Pseudotime");.} %>% 
    ggplot(.,aes(x = UMAP_1,y=UMAP_2,colour=Pseudotime,ylab=''))+
    geom_point(size=0.5)+
    scale_size_continuous(range=c(0,6))+
    # labs(title = i)+
    scale_color_gradientn(colors=c("#916BBF","#54BAB9","#2EB086","#A7D129"),name = "Pseudotime")+
    # scale_color_manual(values=color16)+
    # scale_color_viridis_c()+
    theme_classic()+
    theme(axis.text=element_text(size=10, color="black"),legend.position = "right",
          axis.title=element_text(size=10),title = element_text(size=12),
          plot.title = element_text(hjust = 0.5),
          legend.text=element_text(size=10)))

gene<-c("ALDH1L1","HOPX","VIM","SOX2","PROX1","DCX","SYT1","SV2B")
for (n in 1:length(gene)) {
  # n<-1
  (pd<-plot_genes_in_pseudotime(mycds[gene[n],], color_by = "type",ncol = 3,
                                vertical_jitter=0.14,relative_expr = T)+
     theme(text = element_text(size=15, color="black"),
           axis.text=element_text(size=13, color="black"))+
     NoLegend())
  mycds$Expression<-pd$data$expectation
  (p<-plot_cell_trajectory(mycds, color_by = "Expression",cell_size = 0.5)+
      scale_color_gradientn(colors=c("black","#7986C7","#FFEA85","#F73F52"))+
      labs(title = gene[n])+
      theme(plot.title = element_text(hjust = 0.5),
            axis.text=element_text(size=10, color="black"),legend.position = "right",
            axis.title=element_text(size=10),title = element_text(size=12),
            legend.text=element_text(size=10),legend.title = element_text(size=12)))
  assign(paste("p", n, sep=""),p)
}

tmp_GC1<-mycds@phenoData@data %>% .[.$Minor%in%"GC1",] %>% rownames()
tmp_GC2<-mycds@phenoData@data %>% .[.$Minor%in%"GC2",] %>% rownames()

SeuratObj@meta.data %<>% {.[tmp_GC1,]$Minor<-"GC1";.[tmp_GC2,]$Minor<-"GC2";.}

mycds@phenoData@data %<>% {.[.$State%in%"1"&.$Celltype%in%"GC",]$Minor<-"GC2";.} %<>% 
  {.[.$State%in%"2"&.$Celltype%in%"GC",]$Minor<-"GC1";.} %>% 
  {.$type<-.$Minor;.}
save(mycds,file = "result/4D.root-pNSC.expr.Doublet_remove.final.RData")

color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")

(p2<-plot_cell_trajectory(mycds, color_by = "type",cell_size = 1)+
    scale_color_manual(values=color7)+
    theme(axis.text = element_text(size = 15,colour = "black"),axis.title = element_text(size = 15))+
    guides(color=guide_legend(override.aes = list(size=5))))

ggsave("figure/Fig4C.Neonatal_pseudotime_trajectory.png",p2,width = 8,height = 6,dpi = 300)
ggsave("figure/Fig4C.Neonatal_pseudotime_trajectory.pdf",p2,width = 8,height = 6)

gene<-c("HOPX","VIM","SOX2","DCX","STMN2","SYT1")
for (i in 1:length(gene)) {
  # i<-1
  (pd<-plot_genes_in_pseudotime(mycds[gene[i],], color_by = "type",ncol = 3,
                                vertical_jitter=0.14,relative_expr = T)+
     theme(text = element_text(size=15, color="black"),
           axis.text=element_text(size=13, color="black"))+
     NoLegend())
  mycds$Expression<-pd$data$expectation
  (p<-plot_cell_trajectory(mycds, color_by = "Expression",cell_size = 0.5)+
      scale_color_gradientn(colors = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e")))+
      labs(title = gene[i])+
      theme(plot.title = element_text(hjust = 0.5),
            axis.line.x = element_blank(),axis.line.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text=element_blank(),legend.position = "none",
            axis.title=element_blank(),title = element_blank(),
            legend.text=element_blank(),legend.title = element_blank()))
  assign(paste("p", i, sep=""),p)
}
p<-(p1+p2+p3+p4+p5+p6)+plot_layout(ncol = 3)
ggsave("figure/Fig4C.Neonatal_8feature.expr_trajectory.png",p,width = 18,height = 8,dpi = 300)
ggsave("figure/Fig4C.Neonatal_8feature.expr_trajectory.pdf",p,width = 18,height = 8,dpi = 300)



# —— 8.4 Fig4F neonatal pseudotime branch heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
filename<-paste("result/4D.root-pNSC.expr.Doublet_remove.v6.RData",sep = "")
load(file = filename)
load("result/SeuratObj.Celltype_Minor3.Stage.v11.RData")
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")

(p1<-plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 1)+
    scale_color_gradientn(colors=c("#916BBF","#54BAB9","#2EB086","#A7D129"),name = "Pseudotime"))

beam_res <- BEAM(mycds, branch_point = 1,progenitor_method = "duplicate",cores = 10)
# add progenitor_method = "duplicate"("sequential_split")
# Warning messages:
#   1: In if (progenitor_method == "duplicate") { :
#       the condition has length > 1 and only the first element will be used
#     2: In if (progenitor_method == "sequential_split") { :
#         the condition has length > 1 and only the first element will be used
beam_res <- beam_res[order(beam_res$qval),] %>% .[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
# save(beam_res,mycds_sub,mycds_sub_beam,file = "result/4D.AS_RGL.no_OPC.root-qNSC.expr0.1.Branch.RData")
save(beam_res,mycds_sub,mycds_sub_beam,file = "result/4D.root-pNSC.expr.Doublet_remove.Branch.RData")
# load("result/4D.root-pNSC.expr.Doublet_remove.Branch.RData")

branch_point=1
num_clusters=4
hmcols = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(2000)
ph<-plot_genes_branched_heatmap(mycds_sub_beam,
                                branch_point = branch_point,
                                num_clusters = num_clusters,
                                cores = 30,hmcols = hmcols,
                                use_gene_short_name = T,
                                return_heatmap = T,
                                show_rownames = F)
ggsave("figure/Fig4F.neonatal.Branch.heatmap.png",ph$ph_res,width = 10,height = 8)


top_ann = ph$annotation_col %>% as.data.frame() %>% {colnames(.)<-"State";.} %>%
  HeatmapAnnotation(df = .,show_legend = T,show_annotation_name = F,annotation_label=NULL,
                    col = list(`State` = c("Pre-branch"="#f55c47","Cell fate 1"="#564a4a","Cell fate 2"="#4aa96c")))
p1<-ph$heatmap_matrix %>% .[rownames(mycds_sub_beam[ph$ph$tree_row[["order"]],]),] %>% 
  Heatmap(.,cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = F,
          top_annotation = top_ann,row_km = num_clusters,
          column_split = c(rep(c("A"),100),rep(c("B"),100)),column_title = NULL,
          heatmap_legend_param =list(title="Average Expressed",at=seq(-3,3)),
          col = hmcols)

png("figure/Fig4F.neonatal.Branch.heatmap.png",width = 1000,height = 900)
p1
dev.off()

pdf("figure/Fig4F.neonatal.Branch.heatmap.pdf",width = 10,height = 8)
p1
dev.off()


save(ph,file = "result/Fig4F.neonatal.Branch.heatmap.RData")


allgene<-rownames(mycds_sub_beam[ph$ph_res$tree_row[["order"]],])

tmp_genedf<-cutree(ph$ph_res$tree_row,k=num_clusters) %>% as.data.frame() %>% .[allgene,,drop=F] %>% 
  {colnames(.)<-"Cluster";.}

library( "clusterProfiler")
load(file = '/data1/zhur/proj/scRNA/enrichment_annotation/human_enrichment.RData')
gtf <- read.delim("data/cellranger.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
id2name<-gtf %>% .[.$V3=="gene",] %>% {.$V9<-gsub(".*gene_id ","",.$V9);.$V9<-gsub(";.*gene_name ",",",.$V9);.$V9<-gsub(";.*","",.$V9);.} %>%
  separate(data = ., col = V9, into = c("gene_id","gene_name"), sep = ",") %>% .[,c("gene_id","gene_name")] %>% unique()
minGSSize=3
maxGSSize=500

for (i in unique(tmp_genedf$Cluster)) {
  # i=1
  gene <- tmp_genedf[tmp_genedf$Cluster%in%i,,drop=F] %>% rownames()
  df<-data.frame(gene=gene,cluster=paste("cluster",i,sep = ""))
  filename<-paste("result/Fig4F.4D.branch_heatmap.cluster",i,".gene.xls",sep = "")
  write.table(df,file = filename,quote = F,sep = "\t",row.names = F,col.names = T)
  
  mark_gene<-gene
  ## BP
  term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
  term2gene<-term2gene[!is.na(term2gene$gene_name),]
  term2name <- gobp[,c('GO','Name')]
  
  ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
  query_BP_result<-as.data.frame(ego_BP@result)
  filename<-paste("figure/Fig4F.4D.branch_heatmap.cluster",i,".GOBP.xls",sep = "")
  write.table(query_BP_result,file = filename,quote = F,sep = "\t",row.names = F,col.names = T)
  
  (p2<-plot_genes_branched_pseudotime(mycds[head(intersect(rownames(beam_res),gene),n=15),],color_by = "type",ncol = 5,
                                      branch_point = 1,relative_expr = T,label_by_short_name = F)+
      scale_color_manual(values=color7)+
      labs(title = paste("cluster",i,sep = ""))+
      theme(text = element_text(size=15, color="black"),
            axis.text=element_text(size=13, color="black")))
  filename<-paste("figure/Fig4F.4D.branch_heatmap.cluster",i,".Gene.png",sep = "")
  ggsave(filename = filename,p2,width = 18,height = 9)
}


# —— 8.5 Fig4G neonatal pseudotime branch gene expression -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
load("result/4D.root-pNSC.expr.Doublet_remove.final.RData")
cluster1<-c("HOPX","VIM","CHI3L1","TNC","NRGN","STMN2","SV2B")

gene<-c(cluster1)
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(p<-plot_genes_branched_pseudotime(mycds[gene,],color_by = "type",ncol = 1,
                                   branch_point = 1,relative_expr = T,label_by_short_name = F)+
    scale_color_manual(values=color7)+
    theme(text = element_text(size=15, color="black"),
          axis.text=element_text(size=13, color="black")))

ggsave("figure/Fig4G.neonatal.pseudotime_branch.gene_expression.png",p,width = 5,height = 9,dpi = 300)
ggsave("figure/Fig4G.neonatal.pseudotime_branch.gene_expression.pdf",p,width = 5,height = 9)

# — 9 FigS5 Related to Fig4: Pseudotime reconstruction of the neurogenic lineage development in the neonatal Day 4 human hippocampus. -----------------------------------------------------------
# —— 9.1 FigS5AB Neurogenic_lineage pseudotime trajectory -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

load("result/4D.root-pNSC.expr.Doublet_remove.final.RData")
table(mycds$type)

color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(p1<-plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 1)+
    scale_color_gradient2(low="#1C1C1C",mid="#007965",high ="#ffcc29",midpoint = mean(mycds$Pseudotime)))
(p2<-plot_cell_trajectory(mycds, color_by = "type",cell_size = 1)+
    scale_color_manual(values=color7)+
    NoLegend()+
    # scale_color_gradient2(low="#1C1C1C",mid="#007965",high ="#ffcc29",midpoint = mean(mycds$Pseudotime))+
    facet_wrap(~type,ncol = 5))
(p<-(p1+p2)+plot_layout(widths = c(1,5)))

ggsave("figure/FigS5AB.Neurogenic_lineage.pseudotime_trajectory.png",p,width = 25,height = 5,dpi = 300)
ggsave("figure/FigS5AB.Neurogenic_lineage.pseudotime_trajectory.pdf",p,width = 25,height = 5)

# —— 9.2 FigS5C N1 N2 specific Marker heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
SeuratObj<-tmp_Sobj
DefaultAssay(SeuratObj) <- "RNA"

load("result/4D.root-pNSC.expr.Doublet_remove.v6.RData")
color7<-c("#40407a","#00A087B2","#DC0000B2","#9196f7","#4DBBD5B2")
(p1<-plot_cell_trajectory(mycds, color_by = "State",cell_size = 1))

cell1<-pData(mycds) %>% .[.$type%in%"GC"&.$State%in%"1",] %>% {.$GCtype<-"GC1";.}
cell2<-pData(mycds) %>% .[.$type%in%"GC"&.$State%in%"3",] %>% {.$GCtype<-"GC2";.}
table(mycds$type)
Sub<-SeuratObj %>% subset(.,subset=orig.ident%in%"4D")
meta1<-rbind(cell1[,"GCtype",drop=F],cell2[,"GCtype",drop=F]) %>% {colnames(.)<-"Minor";.}
meta2<-Sub %>% .@meta.data %>% .[.$Celltype%in%setdiff(levels(SeuratObj$Celltype),"GC"),] %>% 
  .[,"Celltype",drop=F] %>% {colnames(.)<-"Minor";.}
Sub<-AddMetaData(Sub,metadata = rbind(meta1,meta2))

color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2","#40407a","#FF82AB",
           "#FFA500","#FF7256","#4682B4","#8B658B","#FFFF00","#00F5FF","#000000",
           "#7E6148B2","#7FFF00","#7FFF00")
(pa<-DimPlot(Sub, reduction = "umap",cols = color16,label.size = 6,
             label = F, pt.size = 0.1,group.by = "Minor",repel = T))

GCobj<-subset(Sub,subset=Celltype%in%"GC")
(pa<-DimPlot(GCobj, reduction = "umap",cols = c("red","black"),label.size = 6,
             label = F, pt.size = 0.1,group.by = "Minor",repel = T))

DefaultAssay(GCobj) <- "RNA"
Idents(GCobj)<-"Minor"

diff1<-FindMarkers(GCobj, ident.1 = "GC1", ident.2 = "GC2", verbose = T) %>% 
  {.$gene<-rownames(.);.} %>% .[.$p_val<0.05,] %>% 
  {.$Diff<-"UP";.$Diff[.$avg_log2FC<0]<-"DOWN";.} %>% 
  .[rev(order(.$Diff,abs(.$avg_log2FC))),]
diff1$pct<-abs(diff1$pct.1-diff1$pct.2)
diff1$sort<-abs(diff1$pct*diff1$avg_log2FC)

G_df<-data.frame()
for (i in c(1,2,3,4)) {
  filename<-paste("result/Fig4F.4D.branch_heatmap.cluster",i,".gene.xls",sep = "")
  G_df<- read.delim(file = filename) %>% rbind(.,G_df)
}
row.names(mycds_sub_beam)[unlist(row_order(p1))] %>% setdiff(G_df$gene)

load("result/Fig4F.neonatal.Branch.heatmap.RData")
pp<-p1
mtx<-pp@ht_list[[" "]]@matrix

row<-rownames(mycds_sub_beam)[unlist(row_order(pp))]
col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,colnames(GCobj))
ht_mtx<-mtx[row,col]
Heatmap(ht_mtx,
        use_raster=T,
        show_column_names = F,show_row_names = F,
        cluster_columns = F,cluster_rows = T,row_dend_reorder = T,
        column_title = "no reordering")

# GC1 UP
Change<-"UP"
tmp_G<-diff1[diff1$Diff%in%Change,]
row<-rownames(mycds_sub_beam)[unlist(row_order(pp))] %>% intersect(tmp_G$gene,.)
col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,colnames(GCobj))
ht_mtx<-mtx[row,col]
ph1<-Heatmap(ht_mtx,
             row_km = 4,
             use_raster=T,
             show_column_names = F,show_row_names = F,
             cluster_columns = F,cluster_rows = T,row_dend_reorder = T,
             column_title = "no reordering") %>% draw()

rownames(ht_mtx)[unlist(row_order(ph1)[2:4])]

tmp_G<-rownames(ht_mtx)[unlist(row_order(ph1)[c("2","3","4")])]
row<-rownames(mycds_sub_beam)[unlist(row_order(pp))] %>% intersect(tmp_G,.)
col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,colnames(GCobj))
ht_mtx<-mtx[row,col]

top_ann = HeatmapAnnotation(df = GCobj@meta.data[col,"Minor",drop=F],
                            show_annotation_name = T,show_legend = T,
                            col = list(Minor = c(`GC1`="#CE1212",`GC2`="#1B1717")))
ph2<-Heatmap(ht_mtx,
             row_km = 3,
             # col = colorRampPalette(c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"))(50),
             # col = c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"),
             # col = c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff",
             #         "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b"),
             # col = rev(c("#c969a1", "#ce4441", "#ee8577", "#eb7926", "#ffbb44", "#859b6c", "#62929a", "#004f63", "#122451")),
             # col = rev(c("#591c19", "#9b332b", "#b64f32", "#d39a2d", "#f7c267", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")),
             # col = colorRampPalette(rev(c("#591c19", "#9b332b", "#b64f32", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")))(8),
             col = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e")),
             top_annotation = top_ann,
             column_split = GCobj@meta.data[col,"Minor",drop=F],
             use_raster=T,
             show_column_names = F,show_row_names = T,
             row_names_gp = gpar(fontsize = 5),
             cluster_columns = F,cluster_rows = T,row_dend_reorder = T,
             column_title = Change) %>% draw()

tmp_up<-data.frame()
for (i in names(row_order(ph2))) {
  # i<-1
  tmp_up<-diff1[rownames(ht_mtx)[unlist(row_order(ph2)[i])],] %>% 
    {.$Cluster<-paste("Cluster",i,sep = "");.} %>% 
    rbind(tmp_up,.)
}

library( "clusterProfiler")
load(file = '/data1/zhur/proj/scRNA/enrichment_annotation/human_enrichment.RData')
gtf <- read.delim("data/cellranger.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
id2name<-gtf %>% .[.$V3=="gene",] %>% {.$V9<-gsub(".*gene_id ","",.$V9);.$V9<-gsub(";.*gene_name ",",",.$V9);.$V9<-gsub(";.*","",.$V9);.} %>%
  separate(data = ., col = V9, into = c("gene_id","gene_name"), sep = ",") %>% .[,c("gene_id","gene_name")] %>% unique()
minGSSize=3
maxGSSize=500
for (i in unique(tmp_up$Cluster)) {
  mark_gene<-tmp_up[tmp_up$Cluster%in%i,]$gene
  ## BP
  term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
  term2gene<-term2gene[!is.na(term2gene$gene_name),]
  term2name <- gobp[,c('GO','Name')]
  
  ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
  query_BP_result<-as.data.frame(ego_BP@result)

}

mark_gene<-tmp_up$gene
## BP
term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
term2gene<-term2gene[!is.na(term2gene$gene_name),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
query_BP_result<-as.data.frame(ego_BP@result)

# GC1 DOWN
Change<-"DOWN"
tmp_G<-diff1[diff1$Diff%in%Change,]
row<-rownames(mycds_sub_beam)[unlist(row_order(pp))] %>% intersect(tmp_G$gene,.)
col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,colnames(GCobj))
ht_mtx<-mtx[row,col]
ph1<-Heatmap(ht_mtx,
             row_km = 4,
             use_raster=T,
             show_column_names = F,show_row_names = F,
             cluster_columns = F,cluster_rows = T,row_dend_reorder = T,
             column_title = "no reordering") %>% draw()

rownames(ht_mtx)[unlist(row_order(ph1)[2:4])]

tmp_G<-rownames(ht_mtx)[unlist(row_order(ph1)[c("2","3","4")])]
row<-rownames(mycds_sub_beam)[unlist(row_order(pp))] %>% intersect(tmp_G,.)
col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,colnames(GCobj))
ht_mtx<-mtx[row,col]

top_ann = HeatmapAnnotation(df = GCobj@meta.data[col,"Minor",drop=F],
                            show_annotation_name = T,show_legend = T,
                            col = list(Minor = c(`GC1`="#CE1212",`GC2`="#1B1717")))
ph2<-Heatmap(ht_mtx,
             row_km = 3,
             # col = colorRampPalette(c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"))(50),
             # col = c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"),
             # col = c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff",
             #         "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b"),
             # col = rev(c("#c969a1", "#ce4441", "#ee8577", "#eb7926", "#ffbb44", "#859b6c", "#62929a", "#004f63", "#122451")),
             # col = rev(c("#591c19", "#9b332b", "#b64f32", "#d39a2d", "#f7c267", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")),
             # col = colorRampPalette(rev(c("#591c19", "#9b332b", "#b64f32", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")))(8),
             col = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e")),
             top_annotation = top_ann,
             column_split = GCobj@meta.data[col,"Minor",drop=F],
             use_raster=T,
             show_column_names = F,show_row_names = T,
             row_names_gp = gpar(fontsize = 5),
             cluster_columns = F,cluster_rows = T,row_dend_reorder = T,
             column_title = Change) %>% draw()


tmp_down<-data.frame()
for (i in names(row_order(ph2))) {
  # i<-1
  tmp_down<-diff1[rownames(ht_mtx)[unlist(row_order(ph2)[i])],] %>% 
    {.$Cluster<-paste("Cluster",i,sep = "");.} %>% 
    rbind(tmp_down,.)
}

library( "clusterProfiler")
load(file = '/data1/zhur/proj/scRNA/enrichment_annotation/human_enrichment.RData')
gtf <- read.delim("data/cellranger.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
id2name<-gtf %>% .[.$V3=="gene",] %>% {.$V9<-gsub(".*gene_id ","",.$V9);.$V9<-gsub(";.*gene_name ",",",.$V9);.$V9<-gsub(";.*","",.$V9);.} %>%
  separate(data = ., col = V9, into = c("gene_id","gene_name"), sep = ",") %>% .[,c("gene_id","gene_name")] %>% unique()
minGSSize=3
maxGSSize=500
for (i in unique(tmp_down$Cluster)) {
  mark_gene<-tmp_down[tmp_down$Cluster%in%i,]$gene
  ## BP
  term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
  term2gene<-term2gene[!is.na(term2gene$gene_name),]
  term2name <- gobp[,c('GO','Name')]
  
  ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
  query_BP_result<-as.data.frame(ego_BP@result)
}

mark_gene<-tmp_down$gene
## BP
term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
term2gene<-term2gene[!is.na(term2gene$gene_name),]
term2name <- gobp[,c('GO','Name')]

ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
query_BP_result<-as.data.frame(ego_BP@result)


tmp_G<-rbind(tmp_up,tmp_down)
row<-rownames(mycds_sub_beam)[unlist(row_order(pp))] %>% intersect(tmp_G$gene,.)
col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,colnames(GCobj))
ht_mtx<-mtx[row,col]

top_ann = HeatmapAnnotation(df = GCobj@meta.data[col,"Minor",drop=F],
                            show_annotation_name = T,show_legend = T,
                            col = list(Minor = c(`GC1`="#CE1212",`GC2`="#1B1717")))
ph2<-Heatmap(ht_mtx,
             row_km = 2,
             # col = colorRampPalette(c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"))(50),
             # col = c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"),
             # col = c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff",
             #         "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b"),
             # col = rev(c("#c969a1", "#ce4441", "#ee8577", "#eb7926", "#ffbb44", "#859b6c", "#62929a", "#004f63", "#122451")),
             # col = rev(c("#591c19", "#9b332b", "#b64f32", "#d39a2d", "#f7c267", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")),
             # col = colorRampPalette(rev(c("#591c19", "#9b332b", "#b64f32", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")))(8),
             col = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e")),
             top_annotation = top_ann,
             column_split = GCobj@meta.data[col,"Minor",drop=F],
             use_raster=T,
             show_column_names = F,show_row_names = T,
             row_names_gp = gpar(fontsize = 5),
             cluster_columns = F,cluster_rows = T,row_dend_reorder = T,
             column_title = Change) %>% draw()

load("result/4D.root-pNSC.expr.Doublet_remove.v6.RData")
load("result/Fig4F.neonatal.Branch.heatmap.RData")

pp<-p1
mtx<-pp@ht_list[[" "]]@matrix

col<-colnames(mtx)[unlist(column_order(pp))] %>% intersect(.,rownames(pData(mycds)[pData(mycds)$Celltype%in%"GC",]))
ht_mtx<-mtx[row,col]

top_ann = HeatmapAnnotation(df = pData(mycds)[col,"Minor",drop=F],
                            show_annotation_name = T,show_legend = T,
                            col = list(Minor = c(`GC1`="#CE1212",`GC2`="#1B1717")))
ph2<-Heatmap(ht_mtx,
             row_km = 2,
             # col = colorRampPalette(c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"))(50),
             # col = c("#b1615c","#d88782","#e3aba7","#edd7d9","#c9c9dd","#9d9dc7","#8282aa","#5a5a83"),
             # col = c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff",
             #         "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b"),
             # col = rev(c("#c969a1", "#ce4441", "#ee8577", "#eb7926", "#ffbb44", "#859b6c", "#62929a", "#004f63", "#122451")),
             # col = rev(c("#591c19", "#9b332b", "#b64f32", "#d39a2d", "#f7c267", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")),
             # col = colorRampPalette(rev(c("#591c19", "#9b332b", "#b64f32", "#b9b9b8", "#8b8b99", "#5d6174", "#41485f", "#262d42")))(8),
             col = rev(c("#810000", "#CE1212", "#F05454", "#ffd06f", "#ffe6b7", "#aadce0", "#72bcd5", "#528fad", "#376795", "#1e466e")),
             top_annotation = top_ann,
             column_split = pData(mycds)[col,"Minor",drop=F],
             use_raster=T,
             show_column_names = F,show_row_names = T,
             row_names_gp = gpar(fontsize = 5),
             cluster_columns = F,cluster_rows = T,row_dend_reorder = T) %>% draw()

pdf(file = paste("figure/FigS5C.N1_N2.specific_Marker.heatmap.pdf",sep = ""),width = 10,height = 9)
ph2
dev.off()
# —— 9.3 FigS5D N1 N2 specific Marker Featureplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample")
rm(list=setdiff(ls(),"tmp_Sobj"))
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
load("result/4D.root-pNSC.expr.Doublet_remove.v6.RData")

DefaultAssay(tmp_Sobj)<-"RNA"
GCobj<-tmp_Sobj[,rownames(pData(mycds)[pData(mycds)$type%in%c("GC1","GC2"),])]
Idents(GCobj)<-"Minor"
gene<-c("NCKAP5","SGCZ","DCC","FAM19A2","FLRT2","RIMS2","NKAIN2","XKR4")

p1<-FeaturePlot(GCobj, features = gene,pt.size = 0.01,
                ncol = 2,reduction = "umap",label = F,label.size = 5)&
  scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))&
  theme_dr(xlength = 0.2,
           ylength = 0.2,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
  theme(panel.grid = element_blank(),legend.position = "none",
        axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
        plot.title = element_text(hjust = 0.5,size = 15,face = "bold"))

ggsave("figure/FigS5D.N1_N2.specific_Marker_Featureplot.png",p1,width = 9,height = 8,dpi = 300)
ggsave("figure/FigS5D.N1_N2.specific_Marker_Featureplot.pdf",p1,width = 9,height = 8,dpi = 300)


# —— 9.4 FigS5E N1 N2 pseudotime branch heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample")
rm(list=ls())

load("result/4D.root-pNSC.expr.Doublet_remove.Branch.RData")
load("result/Fig4F.neonatal.Branch.heatmap.RData")
TFlist <- read.csv("data/motifAnnotations.tf.xls", sep="") %>% {colnames(.)<-"Gene";.} %>% 
  {.$Type<-"TF";.}
allgene<-rownames(mycds_sub_beam[ph$ph_res$tree_row[["order"]],]) %>% 
  intersect(.,TFlist$Gene)

branch_point=1
num_clusters=4
hmcols = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(70)
ph<-plot_genes_branched_heatmap(mycds_sub_beam[allgene,],
                                branch_point = branch_point,
                                num_clusters = num_clusters,
                                hmcols = hmcols,scale_max = 3,scale_min = -3,
                                use_gene_short_name = T,
                                return_heatmap = T,
                                show_rownames = T) 
ggsave("figure/FigS5E.N1_N2.pseudotime_branch.heatmap.png",ph$ph_res,width = 10,height = 8)

# — 10 FigS6 Related to Fig4: Differentially expressed genes along the pseudotime of neurogenic lineage in the neonatal human hippocampus. -----------------------------------------------------------
# —— 10.1 FigS6ABCD Neurogenic_lineage pseudotime branch gene expression -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample")
rm(list=ls())
load("result/4D.root-pNSC.expr.Doublet_remove.final.RData")
cluster1<-c("VIM","GFAP","SOX6","GPC6","CD44","TNC","CHI3L1","EGFR","HOPX")

gene<-c(cluster1)
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(p1<-plot_genes_branched_pseudotime(mycds[gene,],color_by = "type",ncol = 5,cell_size = 0.3,
                                    branch_point = 1,relative_expr = T,label_by_short_name = F)+
    scale_color_manual(values=color7)+
    theme(text = element_text(size=15, color="black"),
          axis.text=element_text(size=13, color="black")))

cluster1<-c("BMP6","CCK","CALB2","NEUROD6","NRGN","STMN2","CALB1")

gene<-c(cluster1)
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(p2<-plot_genes_branched_pseudotime(mycds[gene,],color_by = "type",ncol = 5,cell_size = 0.3,
                                    branch_point = 1,relative_expr = T,label_by_short_name = F)+
    scale_color_manual(values=color7)+
    theme(text = element_text(size=15, color="black"),
          axis.text=element_text(size=13, color="black")))


cluster1<-c("SATB2","SEMA3A","PROM1","SV2B","POU6F2","FOXP2","TNR","SYT6","SYN2")
gene<-c(cluster1)
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(p3<-plot_genes_branched_pseudotime(mycds[gene,],color_by = "type",ncol = 5,cell_size = 0.3,
                                    branch_point = 1,relative_expr = T,label_by_short_name = F)+
    scale_color_manual(values=color7)+
    theme(text = element_text(size=15, color="black"),
          axis.text=element_text(size=13, color="black")))


cluster1<-c("MYT1","DCX","GRIK2","SEMA6D","NRG1","CASC15")
gene<-c(cluster1)
color7<-c("#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
(p4<-plot_genes_branched_pseudotime(mycds[gene,],color_by = "type",ncol = 5,cell_size = 0.3,
                                    branch_point = 1,relative_expr = T,label_by_short_name = F)+
    scale_color_manual(values=color7)+
    theme(text = element_text(size=15, color="black"),
          axis.text=element_text(size=13, color="black")))
p<-p1/p2/p3/p4
ggsave("figure/FigS6ABCD.Neurogenic_lineage.pseudotime_branch.gene_expression.png",p,width = 10,height = 8,dpi = 300)
ggsave("figure/FigS6ABCD.Neurogenic_lineage.pseudotime_branch.gene_expression.pdf",p,width = 10,height = 8)

# — 11 Fig5: Age-dependent molecular alterations of the hippocampal NSCs and NBs -----------------------------------------------------------
# —— 11.1 Fig5A NSCs_NB Age-dependent umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

sub<-subset(tmp_Sobj,Stage %in% c("Neonatal","Middle-age","Aged") & Celltype%in%setdiff(unique(tmp_Sobj$Celltype),c("UN1","UN2"))) %>% 
  {.$Stage<-factor(.$Stage,levels = c("Neonatal","Middle-age","Aged"));.}


cell<-c("GC","IN","NB","aNSC","RGL","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")
color16<-c("#bdc3c7","#bdc3c7","#9196f7","#DC0000B2","#00A087B2",
           "#FFD460","black","#40407a","#bdc3c7","#bdc3c7",
           "#bdc3c7","#bdc3c7","#bdc3c7","#bdc3c7","#bdc3c7","#bdc3c7")

sub$Subtype<-sub$Minor2

sub$Subtype<-factor(sub$Subtype,levels = c("GC","GC1","GC2","IN","NB","aNSC","pNSC","qNSC1","qNSC2","Astrocytes","M-AS",
                                           "OPC","OLG","MG","EC","Pyr","CR","Per"))
levels(sub$Subtype) %<>% {.[.%in%c("GC1","GC2")]<-"GC";.}  
(p1<-DimPlot(sub, reduction = "umap",cols = color16,label.size = 6,raster = T,
             label = F, pt.size = 1,group.by = "Subtype",repel = T) + 
    # theme_bw()+
    theme(text = element_text(size = 15),legend.position = "none",
          plot.title = element_blank(),title = element_blank(),
          axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_blank(),panel.grid = element_blank(),
          strip.switch.pad.grid = element_blank())+
    facet_wrap(~sub$Stage,nrow = 1)+NoLegend())
ggsave("figure/Fig5A.NSCs_NB.Age-dependent.umap.pdf",p1,width = 14,height = 6,dpi = 500)

# —— 11.2 Fig5B NSCs_NB Age-dependent cellpercent -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

SubSobj <- subset(tmp_Sobj,subset=Minor2%in%c("qNSC1","qNSC2","pNSC","aNSC","NB"))


DefaultAssay(SubSobj) <- "RNA"
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")
age<-c("4D","31y","32y","48y","50y","56y","60y","64f","64m","68y")
cluster_split<-split(SubSobj$Minor2,SubSobj$orig.ident)
cluster_split<-lapply(X = cluster_split, FUN = function(x) {
  as.data.frame(prop.table(table(x)))
})
df<-bind_rows(cluster_split, .id = "SubSobj$orig.ident")
colnames(df)<-c("Sample","Cluster","Percent")
df <- df[order(df$Cluster),]
df$Sample<-factor(df$Sample,levels = age)
subdf <- df %>% {.$Group<-as.character(.$Sample);.} %>% 
  {.$Group[which(.$Sample %in% "4D")] <- "Neonatal";.} %>%
  {.$Group[which(.$Sample %in% c("31y","32y"))] <- "Middle-age";.} %>%
  {.$Group[which(.$Sample %in% c("50y","56y","60y","64f","64m","68y"))] <- "Aged";.} %>%
  {.$Group[which(.$Sample %in% "48y")] <- "Injury";.} %>%
  .[.$Group%in%c("Neonatal","Middle-age","Aged"),] %>%
  .[.$Cluster%in%c("qNSC1","qNSC2","pNSC","aNSC","NB"),] %>% 
  tidyr::unite(., "Group2",Cluster, Group,remove=F) 

subdf<-subdf %>% group_by(Group2) %>% dplyr::summarise(mean = mean(Percent),sd=sd(Percent)) %>%
  left_join(subdf,.,by=c("Group2"="Group2")) %>%
  .[,c("Group","Cluster","mean","sd")] %>% unique() %>%
  {.$Group<-factor(.$Group,levels = c("Neonatal","Middle-age","Aged"));.} %>%
  {.$Cluster<-factor(.$Cluster,levels = c("qNSC1","qNSC2","pNSC","aNSC","NB"));.}

(p <- ggplot(data = subdf, mapping = aes(x = Group, y = mean,fill=Group)) + 
    geom_bar(stat= 'identity', position = 'dodge') +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    # scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")) +
    scale_fill_manual(values=c("#f5d662","#EB5230","#800026","#2A4E9D")) +
    ylab(label = "Cell Number percent") + xlab(label = "")+
    scale_x_discrete(expand = c(0, 0.6)) + 
    scale_y_continuous(expand = c(0, 0.002),limits=c(0, 0.7)) + 
    # scale_fill_npg()+
    #ylim(0,0.35) + 
    theme_bw()+
    theme(axis.title.x = element_text(size = 15,colour = "black"),
          axis.text=element_text(size=15,colour = "black"),
          axis.text.x = element_text(angle = 45,hjust = 1),
          axis.title.y = element_text(size=15,colour = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          strip.background = element_blank(), strip.placement = "outside",
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title = element_text(size=18,hjust = 0.5),
          strip.text =  element_text(size = 10)) +
    facet_wrap(~Cluster,nrow = 1)+NoLegend())
ggsave("figure/Fig5B.NSCs_NB.Age-dependent.cellpercent.pdf",p,width = 8,height = 5)

# —— 11.3 Fig5C NSCs_NB Age-dependent marker bubble -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

sub<-subset(tmp_Sobj,subset = Minor2%in%c("qNSC1","qNSC2","pNSC","aNSC","NB")&orig.ident%in%c(setdiff(unique(tmp_Sobj$orig.ident),"48y")))
DefaultAssay(sub)<-"RNA"

gene<-c("LRRC3B","RHOJ","SLC4A4","HOPX","SOX2","VIM","NES","CHI3L1","ASCL1","EOMES","MKI67","STMN2","DCX")
p<-DotPlot(object = sub, features = gene, group.by  = "Stage")
(p3<-p$data %>% 
    {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>%
    # {.$Gene[.$Gene%in%"ENSMMUG00000032158"]<-"NANOG";.} %>% 
    # {.$Gene<-factor(.$Gene,levels = gene2);.} %>% 
    # {.[.$Cluster%in%c("Neonatal","Middle-age","Aged"),];} %>%
    {.$Cluster<-factor(.$Cluster,levels = c("Aged","Middle-age","Neonatal"));.} %>%
    ggplot(.,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(1,10))+
    scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    # scale_color_manual(values=color16)+
    theme_classic()+
    theme(axis.text=element_text(size=12, color="black"),axis.title = element_blank(),
          legend.position = "top") + RotatedAxis())

plotfile = paste('figure/Fig5C.NSCs_NB.Age-dependent_marker.bubble.pdf',sep = "")
ggsave(filename = plotfile,p3,width = 4.8,height = 2.5)
# —— 11.5 Fig5E qNSC1_qNSC2 Age-dependent DEGs violinplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

levels(tmp_Sobj$Stage) %<>% {.[.%in%"Neonatal"]<-"N";.}
levels(tmp_Sobj$Stage) %<>% {.[.%in%"Middle-age"]<-"Ad";.}
levels(tmp_Sobj$Stage) %<>% {.[.%in%"Aged"]<-"Ag";.}

tmp_obj<-subset(tmp_Sobj,subset=Stage%in%c("N","Ad","Ag")&Minor2%in%"qNSC1") %>% 
  {.$Stage<-factor(.$Stage,levels = rev(c("N","Ad","Ag")));.}

DefaultAssay(tmp_obj) <- "RNA"
Idents(tmp_obj)<-"Stage"
gene<-c("LPAR1","TNC","SQSTM1","GPC6","CD44","CHI3L1","SPP1","CASC15","CABLES1","TENM2","CNTN1","RANBP3L","PTGDS","MAP3K5","NDRG2","BRINP1")
(p1<-VlnPlot(tmp_obj, features = gene, group.by = "Stage", split.plot = F,cols = rev(c("#f5d662","#EB5230","#800026")),ncol = 4,
             pt.size = 0, combine = T)&
    theme(axis.title = element_blank(),
          # axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(), 
          axis.text.x = element_text(size = 15,angle = 0,hjust = 0.5))&
    coord_flip()&
    geom_boxplot(alpha=.2,outlier.shape = NA))

tmp_obj<-subset(tmp_Sobj,subset=Stage%in%c("N","Ad","Ag")&Minor2%in%"qNSC2") %>% 
  {.$Stage<-factor(.$Stage,levels = rev(c("N","Ad","Ag")));.}

DefaultAssay(tmp_obj) <- "RNA"
Idents(tmp_obj)<-"Stage"
gene<-c("CASC15","TNC","GAP43","CCND2","SOX2","SOX4","HOPX","VIM","MBP","PLP1","MOBP","SPOCK1","ENPP2","CNTNAP4","CABLES1","CNTN1")
(p2<-VlnPlot(tmp_obj, features = gene, group.by = "Stage", split.plot = F,cols = rev(c("#f5d662","#EB5230","#800026")),ncol = 4,
             pt.size = 0, combine = T)&
    theme(axis.title = element_blank(),
          # axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(), 
          axis.text.x = element_text(size = 15,angle = 0,hjust = 0.5))&
    coord_flip()&
    geom_boxplot(alpha=.2,outlier.shape = NA))


(p<-(p1/p2))
ggsave(filename = "figure/Fig5E.qNSC1_qNSC2.Age-dependent_DEGs.violinplot.pdf",p,width = 8,height = 12)

# —— 11.6 Fig5F qNSC1_qNSC2 Age-dependent GOBP barplot  -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/")
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_df <- read.delim("data/Fig5F.qNSC1_2.GOBP.xls")

tmp_Plist<-list()
for (i in unique(tmp_df$Celltype)) {
  p1<-tmp_df[tmp_df$Celltype%in%i&tmp_df$Group%in%"UP",] %>% 
    {.$Description<-factor(.$Description,levels = rev(unique(.$Description)));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#CF5051",fill = "#CF5051"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =16,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v")
  p2<-tmp_df[tmp_df$Celltype%in%i&tmp_df$Group%in%"DOWN",] %>% 
    {.$Description<-factor(.$Description,levels = rev(unique(.$Description)));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#4C598C",fill = "#4C598C"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =16,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v")
  p<-p1/p2
  tmp_Plist[[i]]<-p
}

p<-Reduce("|",tmp_Plist)+plot_layout(ncol = 2)
ggsave("figure/Fig5F.qNSC1_qNSC2.Age-dependent_GOBP.barplot.pdf",p,width = 12.5,height = 4.5,dpi = 500)
# — 12 FigS7 Related to Fig5: Alterations of the neurogenic lineage related genes in human hippocampus during aging. -----------------------------------------------------------
# —— 12.1 FigS7A pNSC_aNSC_NB marker bubbleplot split by age -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

cell<-c("Astrocytes-qNSC","Primed NSC","aNSC")
color7<-c("#40407a","#00A087B2","#DC0000B2")
scRNAsub<-subset(tmp_Sobj, subset = seurat_cluster_cell_typer  %in% cell)
DefaultAssay(scRNAsub) <- "RNA"

aNSC<-c("HOPX","VIM","SOX2","CCND2","SOX4","ATAD2") %>% rev()
NB<-c("PROX1","GAD1","STMN2","SEMA3C","SSTR2","TACC2","INPP5F","BASP1","NRGN","TERF2IP") %>% rev()
pNSC<-c("GFAP","VIM","PAX6","SOX2","CCND2","SYNE2","CHI3L1") %>% rev()
lt<-list(aNSC,NB,pNSC) %>% {names(.)<-c("aNSC","NB","pNSC");.}

tmp_Plist<-list()
for (i in c(3,1,2)) {
  # i<-1
  (pp<-DotPlot(object = scRNAsub, features = lt[[i]], group.by  = "orig.ident",scale = T,scale.by = "size",
               cols = c("#21aba5", "#e84a5f")) + RotatedAxis())
  result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.}
  (p<-ggplot(result,aes(x =Cluster,y=Gene,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
      geom_point()+
      scale_size_continuous(range=c(3,12))+
      scale_color_gradient2(low="grey90",mid="#FFB830",high ="#FF2442",midpoint = 0.5)+
      # scale_color_manual(values=color16)+
      theme_classic()+
      theme(axis.title.x = element_blank(),
            axis.text=element_text(size=15,colour = "black"),
            axis.title.y = element_text(size = 15,colour = "black"),
            legend.text=element_text(size=15),
            legend.title = element_text(size=15)) + RotatedAxis())
  tmp_Plist[[as.character(i)]]<-p
}

p<-Reduce("+",tmp_Plist)+plot_layout(nrow = 1,guides = "collect")
ggsave("figure/FigS7A.pNSC_aNSC_NB.marker.bubbleplot.split_by_age.pdf",p,width = 16,height = 5)

# — 13 FigS8 Related to Fig5: Differentially expressed genes and enrichment GO terms in pNSC, aNSC, and NB during aging, respectively. -----------------------------------------------------------
# —— 13.1 FigS8ADG pNSCs aNSCs NB Age-dependent DEGs -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

cell<-c("Primed NSC","aNSC","Neuroblast")
color4<-c("#40407a","#00A087B2","#DC0000B2","#9196f7")
samples<-c("Neonatal","Middle-age","Aged")

for (tmp_cell in cell) {
  sub<-tmp_Sobj %>% subset(.,seurat_cluster_cell_typer %in% tmp_cell & Stage %in% samples) %>%
    {.$seurat_cluster_cell_typer<-factor(.$seurat_cluster_cell_typer,levels = tmp_cell);.}
  DefaultAssay(sub) <- "RNA"
  Idents(sub)<-"Stage"
  
  diff1<-FindMarkers(sub, ident.1 = "Aged", ident.2 = "Neonatal", verbose = T) %>% 
    {.$gene<-rownames(.);.} %>% .[.$p_val<0.05,] %>% 
    {.$Diff<-"UP";.$Diff[.$avg_log2FC<0]<-"DOWN";.}
  
  TFlist <- read.csv("data/motifAnnotations.tf.xls", sep="") %>% {colnames(.)<-"Gene";.} %>% 
    {.$Type<-"TF";.}
  UPgene<-diff1 %>% .[.$Diff%in%"UP",] %>% left_join(.,TFlist,by=c("gene"="Gene"))
  DOWNgene<-diff1 %>% .[.$Diff%in%"DOWN",] %>% left_join(.,TFlist,by=c("gene"="Gene"))
  # write.table(UPgene,file = paste("figure/20221027.Review/",tmp_cell,".Neurogenic_niche.no_AS.A_vs_N.Gene_UP.xls",sep = ""),quote = FALSE,sep = "\t",row.names=TRUE, col.names=NA)
  # write.table(DOWNgene,file = paste("figure/20221027.Review/",tmp_cell,".Neurogenic_niche.no_AS.A_vs_N.Gene_DOWN.xls",sep = ""),quote = FALSE,sep = "\t",row.names=TRUE, col.names=NA)
  
  library( "clusterProfiler")
  load(file = '/data1/zhur/proj/scRNA/enrichment_annotation/human_enrichment.RData')
  gtf <- read.delim("data/cellranger.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
  id2name<-gtf %>% .[.$V3=="gene",] %>% {.$V9<-gsub(".*gene_id ","",.$V9);.$V9<-gsub(";.*gene_name ",",",.$V9);.$V9<-gsub(";.*","",.$V9);.} %>%
    separate(data = ., col = V9, into = c("gene_id","gene_name"), sep = ",") %>% .[,c("gene_id","gene_name")] %>% unique()
  minGSSize=3
  maxGSSize=500
  
  for (i in unique(diff1$Diff)) {
    mark_gene<-diff1[diff1$Diff%in%i,]$gene
    ## BP
    term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
    term2gene<-term2gene[!is.na(term2gene$gene_name),]
    term2name <- gobp[,c('GO','Name')]
    
    ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
    query_BP_result<-as.data.frame(ego_BP@result)
    
  }
  up_gene<-diff1[diff1$Diff%in%"UP",] %>% .[grep("^RP[SL]",.$gene,invert = T),] %>% .[grep("^MT-",.$gene,invert = T),] %>% 
    .$gene %>% head(.,n=10)
  down_gene<-diff1[diff1$Diff%in%"DOWN",] %>% .[grep("^RP[SL]",.$gene,invert = T),] %>% .[grep("^MT-",.$gene,invert = T),] %>% 
    .$gene %>% head(.,n=10)
  (pp<-DoHeatmap(sub, features = c(up_gene,down_gene),angle = 0,hjust = 0.5,
                 group.colors =c("#4682B4","#40407a","#DC0000B2"),raster = F)+
      scale_fill_gradientn(colors = colorRampPalette(c("#1a2a6c", "white", "#c21e20"))(10),na.value = "white"))
  ggsave(paste("figure/FigS8.",tmp_cell,".Age-dependent_DEGs.heatmap.pdf",sep = ""),pp,width = 9,height = 6)
  # ggsave(paste("figure/FigS8.",tmp_cell,".Age-dependent_DEGs.heatmap.png",sep = ""),pp,width = 9,height = 6)
}
# —— 13.2 FigS8BCEFHI pNSCs aNSCs NB Age-dependent GOBP -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample")
rm(list=ls())

tmp_df <- read.delim("data/FigS8.3Celltype.GOBP.xls")

for (tmp_cell in unique(tmp_df$Celltype)) {
  p1<-tmp_df[tmp_df$Celltype%in%tmp_cell&tmp_df$Group%in%"UP",] %>% 
    {.$Description<-factor(.$Description,levels = rev(unique(.$Description)));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#CF5051",fill = "#CF5051"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =16,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v")
  
  p2<-tmp_df[tmp_df$Celltype%in%tmp_cell&tmp_df$Group%in%"DOWN",] %>% 
    {.$Description<-factor(.$Description,levels = rev(unique(.$Description)));.} %>% 
    ggplot(., aes(x=Description,y=-log10(pvalue))) + 
    geom_bar(stat = "identity",fill="grey80") +
    geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0))+
    # labs(title = "UP")+
    theme_bw()+
    theme(axis.title.x = element_text(size = 15),
          axis.text = element_text(size=13,colour = "black"),
          axis.text.y = element_blank(),axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#4C598C",fill = "#4C598C"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.1, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =16,face = "bold",colour = "white"))+
    facet_wrap(Group~.,strip.position = "left",dir="v")
  p<-p1/p2
  ggsave(paste("figure/FigS8.",tmp_cell,".Age-dependent_GOBP.bar.pdf",sep = ""),p,width = 7,height = 6)
  
}
# — 14 Fig6: The transcriptomic signatures of the activated neurogenic lineage in the adult human injured hippocampus induced by stroke -----------------------------------------------------------
# —— 14.1 Fig6A NSCs_NB Injury-related umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

sub<-subset(tmp_Sobj,Celltype%in%setdiff(unique(tmp_Sobj$Celltype),c("UN1","UN2")))


cell<-c("GC","IN","NB","aNSC","RGL","AS/qNSC","M-AS",
        "OPC","OLG","MG","EC","Pyr","CR","Per",
        "UN1","UN2")
color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#FFA500","#FF7256","#4682B4",
           "#8B658B","#FFFF00","#00F5FF","#000000","#7E6148B2","#7FFF00")
color16<-c("#bdc3c7","#bdc3c7","#9196f7","#DC0000B2","#00A087B2",
           "#FFD460","black","#40407a","#bdc3c7","#bdc3c7",
           "#bdc3c7","#bdc3c7","#bdc3c7","#bdc3c7","#bdc3c7","#bdc3c7")

sub$Subtype<-sub$Minor2

sub$Subtype<-factor(sub$Subtype,levels = c("GC","GC1","GC2","IN","NB","aNSC","pNSC","qNSC1","qNSC2","Astrocytes","M-AS",
                                           "OPC","OLG","MG","EC","Pyr","CR","Per"))
levels(sub$Subtype) %<>% {.[.%in%c("GC1","GC2")]<-"GC";.}  
(p1<-DimPlot(sub, reduction = "umap",cols = color16,label.size = 6,raster = T,
             label = F, pt.size = 1,group.by = "Subtype",repel = T) + 
    # theme_bw()+
    theme(text = element_text(size = 15),legend.position = "none",
          plot.title = element_blank(),title = element_blank(),
          axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),
          strip.text = element_text(size = 15),strip.background = element_blank(),panel.grid = element_blank(),
          strip.switch.pad.grid = element_blank())+
    facet_wrap(~sub$Stage,nrow = 1)+NoLegend())
ggsave("figure/Fig6A.NSCs_NB.Injury-related.umap.pdf",p1,width = 18,height = 6,dpi = 500)
# —— 14.2 Fig6C NSCs_NB Injury-related Addmodulescore -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))


subseurat<- subset(tmp_Sobj,subset = seurat_cluster_cell_typer %in% c("aNSC","Primed NSC") & orig.ident %in% c("4D","48y")) %>%
  {.$seurat_cluster_cell_typer<-factor(.$seurat_cluster_cell_typer,levels = c("aNSC","Primed NSC"));.} %>% 
  {.$orig.ident<-factor(.$orig.ident,levels = c("4D","48y"));.}
table(subseurat$seurat_cluster_cell_typer,subseurat$orig.ident)
# Astrocytes-qNSC 
# 11071 

load("result/seurat_filtered.v1.RData")
data.filt<-subset(data.filt,cells=Cells(subseurat))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.filt.list <- SplitObject(data.filt, split.by = "orig.ident")
data.filt.list <- lapply(X = data.filt.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst",
                            nfeatures = 2000)
})

NSC.anchors <- FindIntegrationAnchors(object.list = data.filt.list , dims = 1:20)
NSC.combined <- IntegrateData(anchorset = NSC.anchors, dims = 1:20)
DefaultAssay(NSC.combined) <- "RNA"
NSC.combined <- ScaleData(NSC.combined,verbose = T)
NSC.combined <- CellCycleScoring(NSC.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(NSC.combined) <- "integrated"
NSC.combined <- ScaleData(NSC.combined,verbose = T)
# NSC.combined <- ScaleData(NSC.combined, vars.to.regress = "CC.Difference")
NSC.combined <- RunPCA(NSC.combined, npcs = 30, verbose = T)
# t-SNE and Clustering
NSC.combined <- RunUMAP(NSC.combined, seed.use = 42, reduction = "pca", dims = 1:20)
NSC.combined <- RunTSNE(NSC.combined, seed.use = 2, reduction = "pca", dims = 1:20)
NSC.combined <- FindNeighbors(NSC.combined, reduction = "pca", dims = 1:20)
# save(NSC.combined, file="result/NSC.combined.before_cluster.v4.RData")
NSC.combined <- FindClusters(NSC.combined, resolution = seq(0.1,1,0.1))

save(NSC.combined, file="result/NSC.combined.after_cluster.res0.1-1.v4.RData")
load("result/NSC.combined.after_cluster.res0.1-1.v4.RData")

tmp_Idents<-"integrated_snn_res.0.4"

color11<-c("#4682B4","#40407a","#DC0000B2","#FFA500","#9196f7","#6ab04c",
           "#4DBBD5B2","#00A087B2","#7FFF00","#8B658B","#4682B4",
           "#FF82AB","#00F5FF","#000000","#7E6148B2")
color3<-c("#4682B4","#ffd66b","#ec4646")

Idents(NSC.combined)<-tmp_Idents

ortholog <- read.delim("data/Homo_MFA_MMU_Pig_Mouse.Gid_Gname.xls", header=T,stringsAsFactors=FALSE,na.strings = "") %>% 
  .[!is.na(.$Homo_Gname),]
diff_gene<-all.markers 
Marker<-read.delim("data/Fig2EF.Gene_set.sheet3.xls") %>%
  .[.$Celltype%in%c("Reactive_astrocytes_all"),] %>% 
  left_join(.,ortholog[,c("Homo_Gname","Mouse_Gname")],by=c("Marker"="Mouse_Gname")) %>% 
  {.[is.na(.$Homo_Gname),]$Homo_Gname<-.[is.na(.$Homo_Gname),]$Marker;.} %>% 
  unique() %>% subset(.,select=-Marker) %>% {colnames(.)<-c("Celltype","Marker");.} %>% 
  .[.$Marker%in%diff_gene$gene,] %>% unique()

setdiff(Marker$Marker,rownames(tmp_Sobj))
# "Glast"  "CD49f"  "GLT1"   "Cx43"   "SPOT14"

lt<-split(Marker$Marker, Marker$Celltype)

tmp_Plist<-list()

for (i in 1:length(lt)) {
  # i<-1
  gene <- lt[i]
  tmp_Sobj2 <- AddModuleScore(
    object = NSC.combined,
    features = gene,
    assay = "RNA",
    name = names(lt[i]),
    ctrl = 100)
  Idents(tmp_Sobj2)<-tmp_Idents
  pf<-FeaturePlot(tmp_Sobj2, features = paste(names(lt[i]),"1",sep = ""),repel = T,pt.size = 0.001,
                  reduction = "umap",label =T,label.size = 6)+
    scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    theme(plot.title = element_blank())
  pv<-VlnPlot(tmp_Sobj2, features = paste(names(lt[i]),"1",sep = ""),split.plot = F,
              cols = paletteer_d("khroma::soil"),
              pt.size = 0, combine = T)&
    coord_flip()&
    theme(plot.title = element_blank(),axis.title = element_blank(),
          legend.position = "none")&
    geom_boxplot(alpha=.2,outlier.shape = NA)
  p<-((pf|pv)+plot_layout(widths = c(2, 1))) %>% as.ggplot()
  p<-p+labs(title = names(lt[i]))+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.75, linetype="solid"),
          plot.title = element_text(hjust=0.5,size = 16))
  tmp_Plist[[i]]<-p
}
p2<-Reduce("+",tmp_Plist)

ortholog <- read.delim("data/Homo_MFA_MMU_Pig_Mouse.Gid_Gname.xls", header=T,stringsAsFactors=FALSE,na.strings = "") %>% 
  .[!is.na(.$Homo_Gname),]
diff_gene<-all.markers 
Marker<-read.delim("data/Fig2EF.Gene_set.sheet3.xls") %>%
  .[.$Celltype%in%c("RGL_progenitor"),] %>% 
  left_join(.,ortholog[,c("Homo_Gname","Mouse_Gname")],by=c("Marker"="Mouse_Gname")) %>% 
  {.[is.na(.$Homo_Gname),]$Homo_Gname<-.[is.na(.$Homo_Gname),]$Marker;.} %>% 
  unique() %>% subset(.,select=-Marker) %>% {colnames(.)<-c("Celltype","Marker");.} %>% 
  .[.$Marker%in%diff_gene$gene,] %>%
  unique() %>% .[.$Marker!="GFAP",]

setdiff(Marker$Marker,rownames(tmp_Sobj))
# "Glast"  "CD49f"  "GLT1"   "Cx43"   "SPOT14"

lt<-split(Marker$Marker, Marker$Celltype)

tmp_Plist<-list()

for (i in 1:length(lt)) {
  # i<-1
  gene <- lt[i]
  tmp_Sobj2 <- AddModuleScore(
    object = NSC.combined,
    features = gene,
    assay = "RNA",
    name = names(lt[i]),
    ctrl = 100)
  Idents(tmp_Sobj2)<-tmp_Idents
  pf<-FeaturePlot(tmp_Sobj2, features = paste(names(lt[i]),"1",sep = ""),repel = T,pt.size = 0.001,
                  reduction = "umap",label =T,label.size = 6)+
    scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    theme(plot.title = element_blank())
  pv<-VlnPlot(tmp_Sobj2, features = paste(names(lt[i]),"1",sep = ""),split.plot = F,
              cols = paletteer_d("khroma::soil"),
              pt.size = 0, combine = T)&
    coord_flip()&
    theme(plot.title = element_blank(),axis.title = element_blank(),
          legend.position = "none")&
    geom_boxplot(alpha=.2,outlier.shape = NA)
  p<-((pf|pv)+plot_layout(widths = c(2, 1))) %>% as.ggplot()
  p<-p+labs(title = names(lt[i]))+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.75, linetype="solid"),
          plot.title = element_text(hjust=0.5,size = 16))
  tmp_Plist[[i]]<-p
}
p1<-Reduce("+",tmp_Plist)+plot_layout(ncol = 1)


(p<-(p1|p2))
ggsave("figure/Fig6C.NSCs_NB.Injury-related.Addmodulescore.pdf",p,width = 12,height = 5)

# —— 14.3 Fig6D NSCs_NB Injury-related cellpercent -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_final<-subset(tmp_Sobj,subset=Minor3%in%c("qNSC1","qNSC2","pNSC","aNSC","reactive astrocytes","NB")) %>% 
  {.$Minor3<-factor(.$Minor3,levels = c("qNSC1","qNSC2","pNSC","aNSC","reactive astrocytes","NB"));.}
DimPlot(tmp_final,group.by = "Minor3",label = T)


DefaultAssay(tmp_final) <- "RNA"
cluster_split<-split(tmp_final$Minor3,tmp_final$orig.ident)
cluster_split<-lapply(X = cluster_split, FUN = function(x) {
  as.data.frame(prop.table(table(x)))
})
df<-bind_rows(cluster_split, .id = "tmp_final$orig.ident")
colnames(df)<-c("Sample","Cluster","Percent")
df <- df[order(df$Cluster),]
age<-c("4D","31y","32y","48y","50y","56y","60y","64f","64m","68y")
df$Sample<-factor(df$Sample,levels = age)
subdf <- df %>% {.$Group<-as.character(.$Sample);.} %>% 
  {.$Group[which(.$Sample %in% "4D")] <- "Neonatal";.} %>%
  {.$Group[which(.$Sample %in% c("31y","32y"))] <- "Middle-age";.} %>%
  {.$Group[which(.$Sample %in% c("50y","56y","60y","64f","64m","68y"))] <- "Aged";.} %>%
  {.$Group[which(.$Sample %in% "48y")] <- "Injury";.} %>%
  .[.$Group%in%c("Neonatal","Middle-age","Aged","Injury"),] %>%
  .[.$Cluster%in%c("qNSC1","qNSC2","pNSC","aNSC","reactive astrocytes","NB"),] %>% 
  tidyr::unite(., "Group2",Cluster, Group,remove=F) 

subdf<-subdf %>% group_by(Group2) %>% dplyr::summarise(mean = mean(Percent),sd=sd(Percent)) %>%
  left_join(subdf,.,by=c("Group2"="Group2")) %>%
  .[,c("Group","Cluster","mean","sd")] %>% unique() %>%
  {.$Group<-factor(.$Group,levels = c("Neonatal","Middle-age","Aged","Injury"));.} %>%
  {.$Cluster<-factor(.$Cluster,levels = c("qNSC1","qNSC2","pNSC","aNSC","reactive astrocytes","NB"));.}

levels(subdf$Group) %<>% {.[.%in%"Neonatal"]<-"N";.}
levels(subdf$Group) %<>% {.[.%in%"Middle-age"]<-"M";.}
levels(subdf$Group) %<>% {.[.%in%"Aged"]<-"A";.}
levels(subdf$Group) %<>% {.[.%in%"Injury"]<-"I";.}

(p <- ggplot(data = subdf, mapping = aes(x = Group, y = mean,fill=Group)) + 
    geom_bar(stat= 'identity', position = 'dodge') +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    # scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")) +
    scale_fill_manual(values=c("#f5d662","#EB5230","#800026","#2A4E9D")) +
    ylab(label = "Cell Number percent") + xlab(label = "")+
    scale_x_discrete(expand = c(0, 0.6)) + 
    scale_y_continuous(expand = c(0, 0.002),limits=c(0, 0.7)) + 
    # scale_fill_npg()+
    #ylim(0,0.35) + 
    theme_bw()+
    theme(axis.title.x = element_text(size = 15,colour = "black"),
          axis.text=element_text(size=15,colour = "black"),
          axis.text.x = element_text(hjust = 1),
          axis.title.y = element_text(size=15,colour = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          strip.background = element_blank(), strip.placement = "outside",
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.title = element_text(size=18,hjust = 0.5),
          strip.text =  element_text(size = 10)) +
    facet_wrap(~Cluster,nrow = 1)+NoLegend())
ggsave("figure/Fig6D.NSCs_NB.Injury-related.cellpercent.pdf", plot = p,width = 7.8,height = 4)

# —— 14.4 Fig6E NSCs_NB Injury-related pseudotime trajectory -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

load("result/NSC.combined.after_cluster.res0.1-1.v4.RData")
tmp_meta3<-NSC.combined@meta.data %>% .[.$integrated_snn_res.0.4%in%c("0","1","3"),]

SubSobj <- subset(tmp_Sobj,subset=Minor2%in%c("qNSC1","qNSC2","pNSC","NB")) %>% 
  {.@meta.data$Subtype<-as.character(.@meta.data$Minor2);.} %>% 
  {.@meta.data[rownames(tmp_meta3[tmp_meta3$integrated_snn_res.0.4%in%"0",]),]$Subtype<-"reactive NSC1";.} %>%
  {.@meta.data[rownames(tmp_meta3[tmp_meta3$integrated_snn_res.0.4%in%"1",]),]$Subtype<-"reactive NSC2";.} %>%
  {.@meta.data[rownames(tmp_meta3[tmp_meta3$integrated_snn_res.0.4%in%"3",]),]$Subtype<-"reactive NSC3";.}

color7<-c("#40407a","#00A087B2","#DC0000B2","#9196f7","#E7D4B5","#603601")
cell<-c("qNSC1","qNSC2","reactive NSC1","reactive NSC2","reactive NSC3","NB")
scRNAsub<-subset(SubSobj, subset = orig.ident %in% "48y")
scRNAsub$Minor <- factor(scRNAsub$Subtype,levels = cell) 
scRNAsub$orig.ident <- factor(scRNAsub$orig.ident,levels = "48y")
scRNAsub$Minor %>% table()
# pNSC   aNSC     NB    GC1    GC2 
# 1470    944   1045   3592   1852

color7<-c("#E7D4B5","#603601","#7f8fa6","#273c75","#0097e6","#DC0000B2")
(pp<-DimPlot(scRNAsub,group.by = "Minor",cols = color7))

data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(as.matrix(data),
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=50, relative_expr = TRUE)

disp_table <- dispersionTable(mycds)

for (i in c(0.1)) {
  # i<-0.1
  disp.genes <- subset(disp_table, mean_expression >= i & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  mycds <- setOrderingFilter(mycds, disp.genes)
  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
  # Error in graph.adjacency.dense(adjmatrix, mode = mode, weighted = weighted,  : 
  #                                  long vectors not supported yet: ../../src/include/Rinlinedfuns.h:522
  # ncells > 50k
  # https://github.com/cole-trapnell-lab/monocle-release/issues/138
  mycds <- orderCells(mycds)
  mycds$type<-factor(mycds$Minor,levels = cell)
  disp.genes<-as.character(disp.genes)
  filename<-paste("result/48y.noroot.expr",i,".v4.RData",sep = "")
  save(mycds,disp.genes,file = filename)
  
  GM_state <- function(mycds){
    if (length(unique(pData(mycds)$State)) > 1){
      T0_counts <- table(pData(mycds)$State, pData(mycds)$type)[,"qNSC1"]
      return(as.numeric(names(T0_counts)[which
                                         (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }
  mycds <- orderCells(mycds, root_state = GM_state(mycds))
  plot_cell_trajectory(mycds, color_by = "Pseudotime")
  
  filename<-paste("result/48y.qNSC_reactive_NSC.root-pNSC.expr",i,".v4.RData",sep = "")
  save(disp.genes,mycds,file = filename)

}


color7<-c("#E7D4B5","#603601","#FFC312","#0097e6","#DC0000B2")
levels(mycds$type) %<>% {.[.%in%"reactive NSC1"]<-"aNSC";.}
levels(mycds$type) %<>% {.[.%in%c("reactive NSC2","reactive NSC3")]<-"pNSC";.}

(p2<-plot_cell_trajectory(mycds, color_by = "type",cell_size = 1)+
    scale_color_manual(values=color7)+
    theme(axis.text = element_text(size = 15,colour = "black"),axis.title = element_text(size = 15))+
    guides(color=guide_legend(override.aes = list(size=5))))
ggsave("figure/Fig6E.NSCs_NB.Injury-related.pseudotime.trajectory.pdf",p2,width = 6,height = 4)


# —— 14.5 Fig6F NSCs_NB Injury-related pseudotime branch heatmap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

filename<-paste("result/48y.qNSC_reactive_NSC.root-pNSC.expr0.1.v4.RData",sep = "")
load(file = filename)

color7<-c("#E7D4B5","#603601","#7f8fa6","#273c75","#0097e6","#DC0000B2")

(p1<-plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size = 1)+
    scale_color_gradientn(colors=c("#916BBF","#54BAB9","#2EB086","#A7D129"),name = "Pseudotime"))

beam_res <- BEAM(mycds, branch_point = 1,progenitor_method = "duplicate",cores = 10)
# add progenitor_method = "duplicate"("sequential_split")
# Warning messages:
#   1: In if (progenitor_method == "duplicate") { :
#       the condition has length > 1 and only the first element will be used
#     2: In if (progenitor_method == "sequential_split") { :
#         the condition has length > 1 and only the first element will be used
beam_res <- beam_res[order(beam_res$qval),] %>% .[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds[row.names(subset(beam_res, qval < 1e-4)),]
# save(beam_res,mycds_sub,mycds_sub_beam,file = "result/4D.AS_RGL.no_OPC.root-qNSC.expr0.1.Branch.RData")
save(beam_res,mycds_sub_beam,file = "result/48y.qNSC_reactive_NSC.root-pNSC.expr.Branch.RData")
# load("result/4D.root-pNSC.expr.Doublet_remove.Branch.RData")

branch_point=1
num_clusters=4
hmcols = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(2000)

ph<-plot_genes_branched_heatmap(mycds_sub_beam,
                                branch_point = branch_point,
                                num_clusters = num_clusters,
                                cores = 30,
                                use_gene_short_name = T,
                                return_heatmap = T,
                                show_rownames = F)

top_ann = ph$annotation_col %>% as.data.frame() %>% {colnames(.)<-"State";.} %>%
  HeatmapAnnotation(df = .,show_legend = T,show_annotation_name = F,annotation_label=NULL,
                    col = list(`State` = c("Pre-branch"="#f55c47","Cell fate 1"="#564a4a","Cell fate 2"="#4aa96c")))
p1<-ph$heatmap_matrix %>% .[rownames(mycds_sub_beam[ph$ph$tree_row[["order"]],]),] %>% 
  Heatmap(.,cluster_columns = F,cluster_rows = F,show_column_names = F,show_row_names = F,
          top_annotation = top_ann,row_km = num_clusters,
          column_split = c(rep(c("A"),100),rep(c("B"),100)),column_title = NULL,
          heatmap_legend_param =list(title="Average Expressed",at=seq(-3,3)),
          col = hmcols)

png("figure/Fig6F.NSCs_NB.Injury-related.pseudotime_branch.heatmap.png",width = 1000,height = 900)
p1
dev.off()

pdf("figure/Fig6F.NSCs_NB.Injury-related.pseudotime_branch.heatmap.pdf",width = 10,height = 8)
p1
dev.off()


save(ph,file = "result/Fig6F.48y.qNSC_reactive_NSC.branch_heatmap.RData")

allgene<-rownames(mycds_sub_beam[ph$ph_res$tree_row[["order"]],])

tmp_genedf<-cutree(ph$ph_res$tree_row,k=num_clusters) %>% as.data.frame() %>% .[allgene,,drop=F] %>% 
  {colnames(.)<-"Cluster";.}

library( "clusterProfiler")
load(file = '/data1/zhur/proj/scRNA/enrichment_annotation/human_enrichment.RData')
gtf <- read.delim("data/cellranger.gtf", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
id2name<-gtf %>% .[.$V3=="gene",] %>% {.$V9<-gsub(".*gene_id ","",.$V9);.$V9<-gsub(";.*gene_name ",",",.$V9);.$V9<-gsub(";.*","",.$V9);.} %>%
  separate(data = ., col = V9, into = c("gene_id","gene_name"), sep = ",") %>% .[,c("gene_id","gene_name")] %>% unique()
minGSSize=3
maxGSSize=500

for (i in unique(tmp_genedf$Cluster)) {
  # i=1
  gene <- tmp_genedf[tmp_genedf$Cluster%in%i,,drop=F] %>% rownames()
  df<-data.frame(gene=gene,cluster=paste("cluster",i,sep = ""))

  mark_gene<-gene
  ## BP
  term2gene <- gobp[!is.na(gobp$ENSEMBL),] %>% left_join(.,id2name,by=c("ENSEMBL"="gene_id")) %>% .[,c('GO','gene_name')]
  term2gene<-term2gene[!is.na(term2gene$gene_name),]
  term2name <- gobp[,c('GO','Name')]
  
  ego_BP <- enricher(gene = mark_gene,pvalueCutoff = 0.05,pAdjustMethod = "BH",minGSSize = minGSSize,maxGSSize = maxGSSize,qvalueCutoff = 0.2,TERM2GENE = term2gene,TERM2NAME = term2name)
  query_BP_result<-as.data.frame(ego_BP@result)

}
# —— 14.6 Fig6G NSCs_NB Injury-related DEGs Bubble -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"SeuratObj"))
load("result/astrocytes.combined.after_cluster.res0.1-1.v4.RData")
tmp_meta1<-astrocytes.combined@meta.data %>% .[.$integrated_snn_res.0.1%in%c("0","2"),]
tmp_meta2<-SeuratObj@meta.data %>% .[.$Celltype%in%c("NB","aNSC","pNSC"),]

load("result/NSC.combined.after_cluster.res0.1-1.v4.RData")
tmp_meta3<-NSC.combined@meta.data %>% .[.$integrated_snn_res.0.4%in%c("0","1","3"),]

SubSobj <- subset(SeuratObj,cells=c(rownames(tmp_meta1),rownames(tmp_meta2),rownames(tmp_meta3))) %>% 
  {.@meta.data$Subtype<-as.character(.@meta.data$Celltype);.} %>% 
  {.@meta.data[rownames(tmp_meta1[tmp_meta1$integrated_snn_res.0.1%in%"0",]),]$Subtype<-"qNSC1";.} %>% 
  {.@meta.data[rownames(tmp_meta1[tmp_meta1$integrated_snn_res.0.1%in%"2",]),]$Subtype<-"qNSC2";.} %>% 
  {.@meta.data[rownames(tmp_meta3[tmp_meta3$integrated_snn_res.0.4%in%"0",]),]$Subtype<-"reactive NSC1";.} %>%
  {.@meta.data[rownames(tmp_meta3[tmp_meta3$integrated_snn_res.0.4%in%"1",]),]$Subtype<-"reactive NSC2";.} %>%
  {.@meta.data[rownames(tmp_meta3[tmp_meta3$integrated_snn_res.0.4%in%"3",]),]$Subtype<-"reactive NSC3";.}

tmp_cell1<-SubSobj@meta.data %>% 
  .[.$Subtype%in%c("qNSC1","qNSC2","reactive NSC1","reactive NSC2","reactive NSC3","NB")&.$Stage%in%"Injury",]
tmp_cell2<-SubSobj@meta.data %>% 
  .[.$Subtype%in%c("qNSC1","qNSC2","pNSC","aNSC","NB")&.$Stage%in%"Aged",]

tmpSobj <- subset(SeuratObj,cells=c(rownames(tmp_cell1),rownames(tmp_cell2))) %>% 
  {.@meta.data$Group<-"Group";.} %>% 
  {.@meta.data[rownames(tmp_cell1),]$Group<-"Injury";.} %>% 
  {.@meta.data[rownames(tmp_cell2),]$Group<-"Aged";.} %>% 
  {.$Group<-factor(.$Group,levels = c("Aged","Injury"));.}

tmp_dir<-"figure/20221027.Review/Fig6I.Neurogenic_lineage.I_vs_A.Gene.Bubble"
tmp_Idents<-"Stage"

DefaultAssay(tmpSobj) <- "RNA"
Idents(tmpSobj)<-tmp_Idents

gene<-c("TNC","GPC6","CRYAB","GBP2","VMP1","CHI3L1","SPP1","VIM","ITGB1","CYR61","CD63","ACTN4","ELL2","SPOCD1","IFI44L","RASGEF1B","STAT1","ACTN2","GBE1","NAMPT")
pp<-DotPlot(object = tmpSobj, features = gene, group.by  = tmp_Idents,scale = T,scale.by = "size",
            cols = c("#21aba5", "#e84a5f")) + RotatedAxis()
samples<-c("Neonatal","Middle-age","Aged","Injury") %>% rev()
result<-pp$data %>% {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>% 
  {.$Cluster<-factor(.$Cluster,levels = samples);.}
(p1<-ggplot(result,aes(x = Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(2,10))+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    theme_classic()+
    theme(axis.text=element_text(size=13, color="black"),legend.position = "top",axis.title = element_blank()) + RotatedAxis())
ggsave("figure/Fig6G.NSCs_NB.Injury-related.DEGs_Bubble.pdf",p1,width = 7.5,height = 3.5)

# —— 14.7 Fig6H NSCs_NB Injury-related DEGs enrichment -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_df <- read.delim("data/Fig6H.GOBP.xls")

p1<-tmp_df %>% 
  .[rev(order(.$pvalue)),] %>% 
  {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
  ggplot(., aes(x=Description,y=-log10(pvalue))) + 
  geom_bar(stat = "identity",fill="grey80") +
  geom_text(aes(y=0,label=Description),size=7,hjust = 0)+
  coord_flip()+
  scale_y_continuous(expand = c(0, 0))+
  # labs(title = "UP")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=13,colour = "black"),
        axis.text.y = element_blank(),axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),axis.line.y = element_blank(),
        legend.text= element_text(size=13),
        legend.title = element_text(size=13),
        strip.background = element_rect(color="#CF5051",fill = "#CF5051"),
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.wrap = unit(0.1, "mm"),
        plot.title = element_blank(),
        strip.text = element_text(size =16,face = "bold",colour = "white"))+
  facet_wrap(Group~.,strip.position = "left",dir="v")
ggsave("figure/Fig6H.NSCs_NB.Injury-related.DEGs_enrichment.pdf",p1,width = 9.2,height = 4)
# — 15 FigS9 Related to Fig6: Stroke injury induced hippocampal cell apoptosis, astrocyte reactivation and neuronal damages. -----------------------------------------------------------
# —— 15.1 FigS9A GC_IN Injury-related GOBP barplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

df<-read.delim(file = "data/FigS9.GC.select.GO.xls", stringsAsFactors=FALSE) %>% 
  na.omit() %>% {.$number = factor(rev(1:nrow(.)));.}

bubble<-df %>% .[order(.$Count),] %>% 
  {.$Description<-factor(.$Description,levels = unique(.$Description));.}

(p1<-ggplot(data = bubble,aes(x=Description,y=Count,fill=-log10(pvalue))) + 
    geom_bar(aes(x=Description,y=Count,fill=-log10(pvalue)),stat = "identity")+
    geom_text(aes(x=Description,y=0,label=Description),size=7,hjust = 0)+
    # geom_text(data = data[data$Group%in%"down",],aes(x=Description,y=0,label=Description),size=7,hjust = 1)+
    scale_y_continuous(expand = c(0, 0),limit=c(0,max(bubble$Count)))+
    # scale_fill_continuous(values = c("red","blue"))+
    scale_fill_gradient2(low="grey90",mid="#FFC069",high ="#FF2442",
                         limit=c(min(-log10(bubble$pvalue)),max(-log10(bubble$pvalue))),
                         midpoint = median(-log10(bubble$pvalue)))+
    # scale_color_gradientn(colors=c("#7986C7","#FFEA85","#F73F52"))+
    labs(title = "GOBP")+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),axis.text.y = element_blank(),
          axis.line.y = element_blank(),axis.text.x = element_text(size=15,color="black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          strip.background = element_rect(color="#FFD369",fill = "#FFD369"),
          panel.grid = element_blank(),
          # strip.placement = "outside",
          plot.title = element_text(size=18,hjust = 0.5),
          strip.text =  element_text(size =15,face = "bold",colour = "white"))+
    # facet_wrap(Change~.,strip.position = "left",dir="v")+
    coord_flip())


df<-read.delim(file = "data/FigS9.IN.select.GO.xls", stringsAsFactors=FALSE) %>% 
  na.omit() %>% {.$number = factor(rev(1:nrow(.)));.}

bubble<-df %>% .[order(.$Count),] %>% 
  {.$Description<-factor(.$Description,levels = unique(.$Description));.}

(p2<-ggplot(data = bubble,aes(x=Description,y=Count,fill=-log10(pvalue))) + 
    geom_bar(aes(x=Description,y=Count,fill=-log10(pvalue)),stat = "identity")+
    geom_text(aes(x=Description,y=0,label=Description),size=7,hjust = 0)+
    # geom_text(data = data[data$Group%in%"down",],aes(x=Description,y=0,label=Description),size=7,hjust = 1)+
    scale_y_continuous(expand = c(0, 0),limit=c(0,max(bubble$Count)))+
    # scale_fill_continuous(values = c("red","blue"))+
    scale_fill_gradient2(low="grey90",mid="#FFC069",high ="#FF2442",
                         limit=c(min(-log10(bubble$pvalue)),max(-log10(bubble$pvalue))),
                         midpoint = median(-log10(bubble$pvalue)))+
    # scale_color_gradientn(colors=c("#7986C7","#FFEA85","#F73F52"))+
    labs(title = "GOBP")+
    theme_bw()+
    theme(axis.title.x = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),axis.text.y = element_blank(),
          axis.line.y = element_blank(),axis.text.x = element_text(size=15,color="black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15),
          strip.background = element_rect(color="#FFD369",fill = "#FFD369"),
          panel.grid = element_blank(),
          # strip.placement = "outside",
          plot.title = element_text(size=18,hjust = 0.5),
          strip.text =  element_text(size =15,face = "bold",colour = "white"))+
    # facet_wrap(Change~.,strip.position = "left",dir="v")+
    coord_flip())

(p<-p1|p2)
ggsave("figure/FigS9A.GC_IN.Injury-related.GOBP.barplot.pdf",p,width = 18,height = 3)

# —— 15.2 FigS9B GC_IN Injury-related DEGs violinplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
(list<-load("data/SeuratObj.Celltype_Minor3.Stage.v11.RData"))
rm(list=setdiff(ls(),"tmp_Sobj"))

SubObj<-subset(tmp_Sobj,subset=Celltype%in%c("IN")&Stage%in%c("Injury","Aged"))
DimPlot(SubObj,group.by = "Celltype",split.by = "Stage")
(p3<-VlnPlot(SubObj, features = c("CASP3","PARP1","TP53BP1","OGG1","GABARAPL1","BECN1","PRKN","LAMP2"),
             group.by = "Stage", split.plot = F,ncol = 8,
             cols = c("#800026","#2A4E9D"),
             pt.size = 0, combine = T)&coord_flip()&theme(axis.title = element_blank())&geom_boxplot(alpha=.2,outlier.shape = NA))

ggsave("figure/FigS9B.GC_IN.Injury-related.DEGs.violinplot.pdf",p3,width = 18,height = 4)
# — 16 FigS10 Related to Fig6: Initially defined pNSCs and aNSCs from stroke-injured hippocampus contained reactive astrocytes and reactivated NSCs. -----------------------------------------------------------
# —— 16.1 FigS10ABC Injury-related reactive subtype umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

load("result/NSC.combined.after_cluster.res0.1-1.v4.RData")

tmp_Idents<-"integrated_snn_res.0.4"

color11<-c("#4682B4","#40407a","#DC0000B2","#FFA500","#9196f7","#6ab04c",
           "#4DBBD5B2","#00A087B2","#7FFF00","#8B658B","#4682B4",
           "#FF82AB","#00F5FF","#000000","#7E6148B2")
color3<-c("#4682B4","#ffd66b","#ec4646")

Idents(NSC.combined)<-tmp_Idents
(p1 <- DimPlot(NSC.combined, reduction = "umap",
               group.by = tmp_Idents,cols=color11,raster = T,
               pt.size = 1,label = T,label.size = 6))

DefaultAssay(NSC.combined) <- "RNA"
Idents(NSC.combined)<-tmp_Idents
all.markers <- FindAllMarkers(NSC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- all.markers %>% group_by(cluster) %>% do(head(.,n=10))

DefaultAssay(NSC.combined) <- "RNA"
Idents(NSC.combined)<-tmp_Idents
(p2<-DoHeatmap(NSC.combined, features = unique(top10$gene),angle = 0,hjust = 0.5,label = F,group.bar = F,
               group.colors =c("#4682B4","#40407a","#DC0000B2"),raster = T)+
    scale_fill_gradientn(colors = colorRampPalette(c("#1a2a6c", "white", "#c21e20"))(10),na.value = "white"))


cluster_split<-split(NSC.combined@meta.data[[tmp_Idents]],NSC.combined$orig.ident)
cluster_split<-lapply(X = cluster_split, FUN = function(x) {
  as.data.frame(prop.table(table(x)))
})
df<-bind_rows(cluster_split, .id = "NSC.combined$orig.ident")
colnames(df)<-c("Sample","Cluster","Percent")
df <- df[order(df$Cluster),]

(p3 <- ggplot(data = df, mapping = aes(x = Cluster, y = Percent,fill = Sample)) + 
    geom_bar(stat= 'identity', position = 'dodge') +
    scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#b15928")) +
    ylab(label = "Cell Number percent") +
    scale_x_discrete(expand = c(0, 0.001)) + 
    scale_y_continuous(expand = c(0, 0.002)) + 
    #ylim(0,0.35) + 
    theme_classic() +
    theme(axis.title.x = element_text(size = 15,colour = "black"),
          axis.text=element_text(size=15,colour = "black"),
          axis.title.y = element_text(size = 15,colour = "black"),
          legend.text=element_text(size=15),
          legend.title = element_text(size=15)))

(p<-(p1|p2|p3)+plot_layout(widths = c(1,2,1.5)))

ggsave("figure/FigS10ABC.Injury-related.reactive_subtype.umap_heatmap_cellpercent.pdf",p,width = 16,height = 4.5)
# —— 16.2 FigS10D Injury-related reactive subtype Featureplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())
load("result/NSC.combined.after_cluster.res0.1-1.v4.RData")
DefaultAssay(NSC.combined) <- "RNA"
Idents(NSC.combined)<-"integrated_snn_res.0.4"
Gene<-c("VIM","HOPX","LPAR1","SOX2","STMN1","OSMR","TIMP1","LGALS3")
(p<-FeaturePlot(NSC.combined, features = Gene,repel = T,pt.size = 1,ncol = 4,raster = T,
                reduction = "umap",label = F)&
    scale_color_gradientn(colors=rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")))&
    scale_y_continuous(breaks=NULL)&
    scale_x_continuous(breaks=NULL)&
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
          axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank()))

ggsave("figure/FigS10D.Injury-related.reactive_subtype.Featureplot.pdf",p,width = 12,height = 6)
# —— 16.5 FigS10E NSCs_NB Injury-related pseudotime gene expression -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/")
rm(list=setdiff(ls(),"SeuratObj"))

load("result/48y.qNSC_reactive_NSC.root-pNSC.expr0.1.v4.RData")

tmp_title<-"FigS9D.Injury.plot_genes_in_pseudotime"
tmp_dir<-"figure/20221027.Review/FigS9D.Injury.plot_genes_in_pseudotime"

gene<-c("HOPX", "PAX6",  "VIM",  "CD44", "TNC",
        "CHI3L1", "SOX2", "CKAP5", "RANGAP1", "STMN2")
for (i in 1:length(gene)) {
  # i<-1
  (pd<-plot_genes_in_pseudotime(mycds[gene[i],], color_by = "type",ncol = 3,
                                vertical_jitter=0.14,relative_expr = T)+
     theme(text = element_text(size=15, color="black"),
           axis.text=element_text(size=13, color="black"))+
     NoLegend())
  mycds$Expression<-pd$data %>% {rownames(.)<-.$Cell;.} %>% {.[colnames(mycds),]$expectation}
  (p<-plot_cell_trajectory(mycds, color_by = "Expression",cell_size = 0.5)+
      scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
      labs(title = gene[i])+
      theme(title = element_text(size = 25),legend.position = 'none',plot.title = element_text(hjust = 0.5),
            legend.text = element_blank(),axis.line.x = element_blank(),axis.line.y = element_blank(),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
            axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank()))
  assign(paste("p", i, sep=""),p)
}
p<-(p1+p2+p3+p4+p5+p6+p7+p8+p9+p10)+plot_layout(ncol = 5)

ggsave("figure/FigS10E.NSCs_NB.Injury-related.pseudotime.gene_expression.pdf",p,width = 18,height = 5)
# — 17 FigS11 Related to Fig6: Integration of our snRNA-seq dataset with other published data. -----------------------------------------------------------
# —— 17.1 FigS11A Check NSC_NB new marker by Wang's data -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
load("result/SeuratObj.Celltype_Minor3.Stage.v11.RData")

rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_Sub<-subset(tmp_Sobj,subset=Minor3%in%c("M-AS","Astrocytes","reactive astrocytes","qNSC1","qNSC2","pNSC","aNSC","NB"))
DefaultAssay(tmp_Sub) <- "RNA"

gene<-c("ETNPPL","STMN1","STMN2") %>% rev()
cell<-c("M-AS","Astrocytes","reactive astrocytes","qNSC1","qNSC2","pNSC","aNSC","NB") %>% rev()

p<-DotPlot(object = tmp_Sub, features = gene, group.by  = "Minor3")
(p3<-p$data %>% 
    {colnames(.)<-c("avg.exp","Percent Expressed","Gene","Cluster","Average Expression");.} %>%
    {.$Cluster<-as.character(.$Cluster);.} %>%
    {.$Cluster<-factor(.$Cluster,levels = cell);.} %>%
    ggplot(.,aes(x =Gene,y=Cluster,size=`Percent Expressed`,colour=`Average Expression`,ylab=''))+
    geom_point()+
    scale_size_continuous(range=c(1,12))+
    scale_color_gradientn(colours = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    # scale_color_manual(values=color16)+
    theme_classic()+
    theme(axis.text=element_text(size=11, color="black")) + RotatedAxis())
ggsave(paste("figure/FigS11A.Check_NSC_NB.new_marker.by_Wang_data.pdf",sep = ""),p3,width = 4.5,height = 4)
# —— 17.1 FigS11B Integration Zhou's refdata umap -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

load(file = paste("result/Song_hongjun_et_al.v2.data.filt.RData",sep = ""))
refmeta<-data.filt@meta.data %>% .[,c("orig.ident"),drop=F] %>% {colnames(.)<-c("Samples");.} %>%
  {.$Group<-"Song_hongjun_et_al";.$Celltype<-"Song_hongjun_et_al";.}
refCounts<-data.filt@assays[["RNA"]]@counts

load("result/SeuratObj.Celltype.Stage.v11.RData")
meta<-tmp_Sobj@meta.data %>% .[,c("orig.ident","Celltype")] %>% {colnames(.)<-c("Samples","Celltype");.} %>%
  {.$Group<-"Own";.}
Counts<-tmp_Sobj@assays[["RNA"]]@counts

tmp_gene<-intersect(rownames(Counts),rownames(refCounts))

color3<-c("#1a2a6c","#83ccd2")
gene_list <- list(`Song_hongjun_et_al` = rownames(refCounts),
                  `This study` = rownames(Counts))

(p1<-plot(euler(gene_list, shape = "ellipse"), quantities = TRUE,
          fills = list(fill = color3, alpha = 0.6)) %>% as.ggplot())


final_count<-cbind(Counts[tmp_gene,],refCounts[tmp_gene,])
final_meta<-rbind(meta,refmeta)

final_obj<-CreateSeuratObject(final_count,meta.data = final_meta)

data.list <- SplitObject(final_obj, split.by = "Samples")
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = data.list)
seurat.anchors <- FindIntegrationAnchors(object.list = data.list,reference = 1)
# save(seurat.anchors,file = paste("figure/",tmp_out,".Integrate_own.seurat.anchors.RData",sep = ""))
seurat.combined <- IntegrateData(anchorset = seurat.anchors)
DefaultAssay(seurat.combined) <- "integrated"
seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined, npcs = 30, verbose = FALSE)
seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:30)
seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:30)
seurat.combined <- FindClusters(seurat.combined, resolution = 0.1)
save(seurat.combined,file = paste("figure/Song_hongjun_et_al.v2.Integrate_own.seurat.combined.RData",sep = ""))

tmp_celltype<-c("GC","IN","NB","aNSC","pNSC","AS/qNSC","M-AS",
                "OPC","OLG","MG","EC","Pyr","CR","Per",
                "UN1","UN2","Song_hongjun_et_al")
seurat.combined$Celltype<-factor(seurat.combined$Celltype,levels = tmp_celltype)

color<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2","#40407a","#FF82AB",
         "#FFA500","#FF7256","#4682B4","#8B658B","#FFFF00","#00F5FF","#000000",
         "#7E6148B2","#7FFF00","grey40")
p1<-DimPlot(seurat.combined, reduction = "umap", group.by = "Celltype",
            label = F,repel = T,cols = color,label.size = 6,raster=FALSE)
# color<-c("grey40","#DC0000B2")
p2<-DimPlot(seurat.combined, reduction = "umap", group.by = "Group",
            label = F,repel = T,label.size = 6,raster=FALSE)

(p<-p1|p2)
ggsave(paste("figure/Song_hongjun_et_al.v2.display.cca.umap.png",sep = ""),p,width = 18,height = 9)
ggsave(paste("figure/FigS11B.Integration.Zhou_refdata.umap.png",sep = ""),p,width = 18,height = 9)

# —— 17.1 FigS11CD Integration Zhou's refdata featureplot -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample")
rm(list=setdiff(ls(),"SeuratObj"))

load(file = paste("figure/Song_hongjun_et_al.v2.Integrate_own.seurat.combined.RData",sep = ""))

tmp_Sobj<-seurat.combined
tmp_Idents<-"Celltype"
DefaultAssay(tmp_Sobj) <- "RNA"

Idents(tmp_Sobj)<-tmp_Idents
tmp_gene<-c("VIM","TNC","CHI3L1","SQSTM1")
p<-FeaturePlot(tmp_Sobj, features = tmp_gene,repel = T,pt.size = 0.01,ncol = 2,raster = F,split.by = "Group",
               cols = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")),reduction = "umap",label = F)&
  # scale_y_continuous(breaks=NULL)&
  # scale_x_continuous(breaks=NULL)&
  theme(legend.key.size = unit(0.8,'cm'),legend.text = element_blank(),
        # title = element_text(size = 0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())&
  NoLegend()

ggsave(filename = "figure/FigS11CD.Integration.Zhou_refdata.featureplot.part1.png",p,width = 12,height = 22)


Idents(tmp_Sobj)<-tmp_Idents
tmp_gene<-c("SOX4","SOX11","STMN1","NRGN")
p<-FeaturePlot(tmp_Sobj, features = tmp_gene,repel = T,pt.size = 0.01,ncol = 2,raster = F,split.by = "Group",
               cols = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","grey90")),reduction = "umap",label = F)&
  # scale_y_continuous(breaks=NULL)&
  # scale_x_continuous(breaks=NULL)&
  theme(legend.key.size = unit(0.8,'cm'),legend.text = element_blank(),
        # title = element_text(size = 0),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())&
  NoLegend()

ggsave(filename = "figure/FigS11CD.Integration.Zhou_refdata.featureplot.part2.png",p,width = 12,height = 22)

# —— 17.1 FigS11E Multimodal reference mapping wang_et_al -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=setdiff(ls(),"tmp_Sobj"))

tmp_title<-"Wang_et_al"

load("result/SeuratObj.Celltype_Minor3.Stage.v11.RData")
ref_obj<-subset(tmp_Sobj,subset=Minor3%in%c("M-AS","Astrocytes","qNSC1","qNSC2","pNSC","aNSC","NB"))
# ref_obj<-tmp_Sobj

color16<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2","#40407a","#FF82AB",
           "#FFA500","#FF7256","#4682B4","#8B658B","#FFFF00","#00F5FF","#000000",
           "#7E6148B2","#7FFF00")
color_list<-data.frame(celltype=levels(ref_obj$Major),color=color16,row.names = levels(ref_obj$Major))
DimPlot(ref_obj,group.by = "Major",cols = color_list[levels(ref_obj$Major),]$color)
DimPlot(ref_obj,group.by = "Minor3",cols = color_list$color)

ref_obj[["umap.new"]] <- CreateDimReducObject(embeddings = ref_obj[["umap"]]@cell.embeddings, key = "UMAPnew_")

# set UMAP models
umap.new.model <- list()
umap.new.model$n_epochs <- 500
umap.new.model$alpha <-1
umap.new.model$method <- "umap"
umap.new.model$negative_sample_rate <- 5
umap.new.model$gamma <- 1
umap.new.model$approx_pow <- 0
umap.new.model$metric$cosine <- list()
umap.new.model$embedding <- ref_obj[["umap.new"]]@cell.embeddings
ab_param <- uwot:::find_ab_params(spread = 1, min_dist = 0.3)
umap.new.model$a <- ab_param["a"]
umap.new.model$b <- ab_param["b"]
ref_obj[["umap.new"]]@misc$model <- umap.new.model

DimPlot(ref_obj,group.by = "Minor3",reduction = "umap.new")

load(file = paste("figure/Wang_xiaoqun_et_al.Integrate_own.seurat.combined.RData",sep = ""))

tmp_celltype<-c(grep("^ref_",unique(seurat.combined$Celltype),value = T,invert = T),
                grep("^ref_",unique(seurat.combined$Celltype),value = T))
seurat.combined$Celltype<-factor(seurat.combined$Celltype,levels = tmp_celltype)
query_obj<-subset(seurat.combined,subset=Group%in%"Wang_xiaoqun_et_al"&Celltype%in%c("ref_NSC/As","ref_ImmN"))

color<-c("#4DBBD5B2","#6ab04c","#9196f7","#DC0000B2","#00A087B2","#40407a","#FF82AB",
         "#FFA500","#FF7256","#4682B4","#8B658B","#FFFF00","#00F5FF","#000000",
         "#7E6148B2","#7FFF00",
         rep("grey80",length(grep("^ref_",unique(query_obj$Celltype),value = T))))
DimPlot(query_obj, reduction = "umap", group.by = "Celltype",
        label = T,repel = T,cols = color,label.size = 6,raster=FALSE)


anchors <- FindTransferAnchors(
  reference = ref_obj,
  query = query_obj,
  k.filter = NA ,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:20)

# save(anchors,query_obj,ref_obj,file = "result/tmp.RData")
# load("result/tmp.RData")

test_obj <- MapQuery(
  anchorset = anchors,
  query = query_obj,
  reference = ref_obj,
  refdata = list(
    Major = "Major",
    Minor = "Minor3"),
  reference.reduction = "pca",
  reduction.model = "umap.new"
)

save(test_obj,ref_obj,file = paste("result/",tmp_title,"_mapping.RData",sep = ""))

ref_obj$id <- 'reference'
test_obj$id <- 'query'
# refquery <- merge(subset(ref_obj,subset=Minor3%in%c("pNSC","qNSC1","qNSC2","Astrocytes","M-AS")),
#                   subset(test_obj,subset=Type%in%c("Astro")&predicted.Minor%in%c("pNSC","qNSC1","Astrocytes","M-AS")))
# refquery[["umap"]] <- merge(subset(ref_obj,subset=Minor3%in%c("pNSC","qNSC1","qNSC2","Astrocytes","M-AS"))[["umap"]],
#                             subset(test_obj,subset=Type%in%c("Astro")&predicted.Minor%in%c("pNSC","qNSC1","Astrocytes","M-AS"))[["ref.umap"]])

refquery <- merge(ref_obj,test_obj)
refquery[["umap"]] <- merge(ref_obj[["umap"]],
                            test_obj[["ref.umap"]])


Idents(refquery)<-"predicted.Minor"
(p1<-DimPlot(refquery, group.by = 'id', shuffle = TRUE,label=T,cols = c("#f39c12","grey60")))
(p2<-FeaturePlot(refquery, features = "predicted.Minor.score",label = T)+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e"))))

p1|p2

Idents(refquery)<-"RNA"
tmp_cell<-p2$data
# tmp_cell<-p2$data %>% .[.$UMAP_1>2&.$UMAP_2<5&.$UMAP_2>-5,]
(p4<-subset(refquery,cells=rownames(tmp_cell)) %>% 
    FeaturePlot(., features = "predicted.Minor.score",label = T,raster = T)+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

(p5<-subset(refquery,cells=rownames(tmp_cell)) %>% 
    DimPlot(., group.by = 'id', shuffle = TRUE,label=T,cols = c("#f39c12","grey50"),raster = T)+
    # scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

(pp1<-(p5/p4))
ggsave(filename = paste("figure/FigS11E.Multimodal_reference_mapping.wang_et_al.pdf",sep = ""),
       pp1,width = 5,height = 8)
# —— 17.1 FigS11F Multimodal reference mapping Franjic_et_al -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

tmp_title<-"Franjic_Astro"
load(paste("result/",tmp_title,".Integrat.seurat.combined.RData",sep = ""))

seurat.combined@meta.data %<>% 
  mutate(Group=case_when(
    Batch %in% c("31y","32y","48y","4D","50y","56y","60y","64f","64m","68y") ~ "This Study",
    Batch %in% c("HSB179","HSB181","HSB231","HSB237","HSB282","HSB628") ~ "Franjic_et_al.Neuron.2022"
  ))
seurat.combined$Group<-factor(seurat.combined$Group,levels = c("This Study","Franjic_et_al.Neuron.2022"))
DimPlot(subset(seurat.combined,subset=Raw_Celltype%in%c("pNSC","qNSC1","qNSC2","Astrocytes","M-AS","Astro")),
        group.by = "Raw_Celltype",label = T,split.by = "Group")


load("result/Franjic_mapping.RData")

ref_obj$id <- 'reference'
test_obj$id <- 'query'
refquery <- merge(subset(ref_obj,subset=Minor3%in%c("pNSC","qNSC1","qNSC2","Astrocytes","M-AS","aNSC","NB")),
                  subset(test_obj,subset=Type%in%c("Astro")&predicted.Minor%in%c("pNSC","qNSC1","Astrocytes","M-AS","aNSC","NB")))
refquery[["umap"]] <- merge(subset(ref_obj,subset=Minor3%in%c("pNSC","qNSC1","qNSC2","Astrocytes","M-AS","aNSC","NB"))[["umap"]],
                            subset(test_obj,subset=Type%in%c("Astro")&predicted.Minor%in%c("pNSC","qNSC1","Astrocytes","M-AS","aNSC","NB"))[["ref.umap"]])

Idents(refquery)<-"predicted.Minor"
(p1<-DimPlot(refquery, group.by = 'id', shuffle = TRUE,label=T,cols = c("#f39c12","grey60")))
(p2<-FeaturePlot(refquery, features = "predicted.Minor.score",label = T)+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e"))))

p1|p2

Idents(refquery)<-"RNA"
tmp_cell<-p2$data
# tmp_cell<-p2$data %>% .[.$UMAP_1>2&.$UMAP_2<5&.$UMAP_2>-5,]
(p4<-subset(refquery,cells=rownames(tmp_cell)) %>% 
    FeaturePlot(., features = "predicted.Minor.score",label = T,raster = T)+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

(p5<-subset(refquery,cells=rownames(tmp_cell)) %>% 
    DimPlot(., group.by = 'id', shuffle = TRUE,label=T,cols = c("#f39c12","grey50"),raster = T)+
    # scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

(pp1<-(p5/p4))
ggsave(filename = paste("figure/FigS11F.Multimodal_reference_mapping_Franjic_et_al.pdf",sep = ""),
       pp1,width = 5,height = 8)
# —— 17.1 FigS11G Multimodal reference mapping Ayhan_et_al -----------------------------------------------------------
setwd("/data1/zhur/proj/scRNA/ltq/sample10_v9_aggr_change_sample/Code_upload")
rm(list=ls())

tmp_title<-"Fatma_Ayhan"
load(paste("result/",tmp_title,".Integrat.seurat.combined.RData",sep = ""))

seurat.combined@meta.data %<>% 
  mutate(Group=case_when(
    Batch %in% c("31y","32y","48y","4D","50y","56y","60y","64f","64m","68y") ~ "This Study",
    Batch %in% c("Batch2","Batch3","Batch5","Batch6") ~ "Ayhan_et_al.Neuron.2021"
  ))
seurat.combined$Group<-factor(seurat.combined$Group,levels = c("This Study","Ayhan_et_al.Neuron.2021"))
DimPlot(seurat.combined,group.by = "Raw_Celltype",label = T,split.by = "Group")


load("result/Fatma_Ayhan_mapping.RData")

ref_obj$id <- 'reference'
test_obj$id <- 'query'
refquery <- merge(subset(ref_obj,subset=Minor3%in%c("qNSC1","qNSC2","pNSC","Astrocytes","M-AS","aNSC","NB")),
                  subset(test_obj,subset=predicted.Minor%in%c("Astrocytes","qNSC1","aNSC")&Cluster%in%c("Astro1","Astro2")))
refquery[["umap"]] <- merge(subset(ref_obj,subset=Minor3%in%c("qNSC1","qNSC2","pNSC","Astrocytes","M-AS","aNSC","NB"))[["umap"]],
                            subset(test_obj,subset=predicted.Minor%in%c("Astrocytes","qNSC1","aNSC")&Cluster%in%c("Astro1","Astro2"))[["ref.umap"]])

Idents(refquery)<-"predicted.Minor"
(p1<-DimPlot(refquery, group.by = 'id', shuffle = TRUE,label=T,cols = c("#f39c12","grey60")))
(p2<-FeaturePlot(refquery, features = "predicted.Minor.score",label = T)+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e"))))

p1|p2

Idents(refquery)<-"RNA"
tmp_cell<-p2$data
# tmp_cell<-p2$data %>% .[.$UMAP_1>2&.$UMAP_2<5&.$UMAP_2>-5,]
(p4<-subset(refquery,cells=rownames(tmp_cell)) %>% 
    FeaturePlot(., features = "predicted.Minor.score",label = T,raster = T)+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

(p5<-subset(refquery,cells=rownames(tmp_cell)) %>% 
    DimPlot(., group.by = 'id', shuffle = TRUE,label=T,cols = c("#f39c12","grey50"),raster = T)+
    # scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    tidydr::theme_dr(xlength = 0.2,ylength = 0.2,
                     arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed"))&
    theme(panel.grid = element_blank(),plot.title = element_blank(),
          # legend.position = "none",
          axis.title.x = element_text(hjust = 0),axis.title.y = element_text(hjust = 0),
          axis.title.y.right = element_text(size = 15,hjust = 0.5,face = "bold")))

(pp1<-(p5/p4))
ggsave(filename = paste("figure/FigS11G.Multimodal_reference_mapping.Ayhan_et_al.pdf",sep = ""),
       pp1,width = 5,height = 8)
