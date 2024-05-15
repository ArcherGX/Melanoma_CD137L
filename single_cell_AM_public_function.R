# ========================================================================================================
# Fig. 1 G and H
# ========================================================================================================
# 确保目录存在
checkDir(paste0(res_home, "Figure/phase5.single_spatial/"))
checkDir(paste0(res_home, "Table/phase5.single_spatial/"))
checkDir(paste0(res_home, "RData/phase5.single_spatial/"))

seurat_obj_AM <-  readRDS("~/Data/Melanoma_SC_public/Acral_Melanoma/AcralMelanoma_WholeTissueHarmony_GSE189889.rds.gz")
# ================================================================================================================
# 原文方法描述为
# For those, we merged these clusters from different curated cell types into a single
# "main category". Each cell group (or main category) is then assigned a
# broad cell category based on the combination of SingleR predictions,
# mostly using the BlueprintEncode panel and followed by curation.
#
# Two small clusters of cells (of a total of 220 cells) with potential doublet signatures
# identified were further excluded from the analyses, resulting in a final 35,988 cells included in the analyses.
# ================================================================================================================
#
library(Seurat)
sampletype <- seurat_obj_AM@meta.data[, "sample"]
sampletype <- gsub("Acral-01", "metastatic_node", sampletype)
sampletype <- gsub("MOF-AM[2-4]", "primary", sampletype)
sampletype <- gsub("MOF-AM5", "metastatic", sampletype)
sampletype <- gsub("MOF-AM6-RNA", "primary", sampletype)
sampletype <- gsub("MOF-AM7-RNA", "metastatic_node", sampletype)
sampletype <- gsub("MOF-AM8-Node-RNA", "metastatic_node", sampletype)
sampletype <- gsub("MOF-AM8-Toe-RNA", "primary", sampletype)
seurat_obj_AM@meta.data[, "SampleType"] <- sampletype

# ================================================================================================================
# 根据结果进行细胞cluster删除，以及重新注释
# ================================================================================================================
# 可视化blueprint.main注释
outfile_tmp <- paste0(res_home, "8.single_spatial/Figure/seurat_obj_AM.RNA_snn_res.0.5.pdf")
ggobj <- DimPlot(seurat_obj_AM, label = TRUE,  group.by ="RNA_snn_res.0.5", raster=TRUE)
ggsave(filename=outfile_tmp, plot=ggobj, width = 9, height =7)
# 可视化blueprint.main注释
outfile_tmp <- paste0(res_home, "8.single_spatial/Figure/seurat_obj_AM.blueprint.main.pdf")
ggobj <- DimPlot(seurat_obj_AM, label = TRUE, cols=color_space_main, group.by ="blueprint.main", raster=TRUE)
ggsave(filename=outfile_tmp, plot=ggobj, width = 9, height =7)
# 删除细胞数<100的簇
anno_final <- as.character(seurat_obj_AM@meta.data[, "RNA_snn_res.0.5"])
anno_final[which(anno_final=="0")] <- "Melanocytes"
anno_final[which(anno_final=="1")] <- "Melanocytes"
anno_final[which(anno_final=="2")] <- "CD4+ T cells"
anno_final[which(anno_final=="3")] <- "Melanocytes"
anno_final[which(anno_final=="4")] <- "CD8+ T cells"
anno_final[which(anno_final=="5")] <- "Melanocytes"
anno_final[which(anno_final=="6")] <- "Fibroblasts"
anno_final[which(anno_final=="7")] <- "Endothelial cells"
anno_final[which(anno_final=="8")] <- "Monocytes"
anno_final[which(anno_final=="9")] <- "CD4+ T cells"
anno_final[which(anno_final=="10")] <- "Melanocytes"
anno_final[which(anno_final=="11")] <- "B cells"
anno_final[which(anno_final=="12")] <- "Melanocytes"
anno_final[which(anno_final=="13")] <- "NK"
anno_final[which(anno_final=="14")] <- "Fibroblasts"
anno_final[which(anno_final=="15")] <- "Melanocytes"
anno_final[which(anno_final=="16")] <- "B cells"
anno_final[which(anno_final=="17")] <- "Melanocytes"
anno_final[which(anno_final=="18")] <- "Endothelial cells"
anno_final[which(anno_final=="19")] <- "Melanocytes"
anno_final[which(anno_final=="20")] <- "Melanocytes"
anno_final[which(anno_final=="21")] <- "RM"
anno_final[which(anno_final=="22")] <- "RM"
anno_final[which(anno_final=="23")] <- "Melanocytes"
seurat_obj_AM@meta.data[, "anno_final"] <- anno_final
# ex_cells <- rownames(seurat_obj_AM@meta.data)[!(seurat_obj_AM@meta.data[, "blueprint.main"] %in% names(which(anno_main<100)))]
seurat_obj_AM_sub <- subset(seurat_obj_AM, anno_final!="RM")
save(seurat_obj_AM_sub, file=paste0(res_home, "RData/phase5.single_spatial/seurat_obj_AM_sub.RData"))
# 重注释
# seurat_obj_AM_sub@meta.data[, "blueprint.main"] <- gsub("Adipocytes", "Fibroblasts", seurat_obj_AM_sub@meta.data[, "blueprint.main"])
# anno_final <- seurat_obj_AM_sub@meta.data[, "blueprint.main"]


anno_main <- table(seurat_obj_AM@meta.data[, "blueprint.main"])
color_space_main <- getColor_palette(pal_custom="jimmy", col_num=length(anno_main), pal_name=NULL)
color_space_main <- color_space_main[1:length(anno_main)]
names(color_space_main) <- names(anno_main)
# 可视化blueprint.main注释
outfile_tmp <- paste0(res_home, "Figure/8.single_spatial/GSE189889.blueprint.main.pdf")
ggobj <- DimPlot(seurat_obj_AM, label = TRUE, cols=color_space_main, group.by ="blueprint.main", split.by="SampleType")
ggsave(filename=outfile_tmp, plot=ggobj, width = 10, height =7)
# 可视化blueprint.main注释下，基因表达分布 -- VlnPlot
outfile_tmp <- paste0(res_home, "Figure/8.single_spatial/seurat_obj_AM.interest_genes.blueprint.main.VlnPlot.pdf")
ggobj <- VlnPlot(object=seurat_obj_AM, features = interest_genes, pt.size = 1, cols = color_space_main, group.by="blueprint.main")
ggsave(filename=outfile_tmp, plot=ggobj, width = 15, height =7)
# 可视化blueprint.main注释下，基因表达分布 -- VlnPlot
outfile_tmp <- paste0(res_home, "Figure/8.single_spatial/GSE189889.interest_genes.AllTypes.VlnPlot.pdf")
ggobj <- VlnPlot(object=seurat_obj_AM, features = interest_genes, pt.size = 1, 
  cols = color_space_main, group.by="blueprint.main", split.by="SampleType", raster=TRUE)
ggsave(filename=outfile_tmp, plot=ggobj, width = 15, height =7)
# 可视化blueprint.main注释下，基因表达分布 -- DotPlot
outfile_tmp <- paste0(res_home, "Figure/8.single_spatial/GSE189889.interest_genes.blueprint.main.DotPlot.pdf")
ggobj <- DotPlot(object=seurat_obj_AM, features = interest_genes, cols = c("lightgrey", "blue"), group.by="blueprint.main") +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5))
ggsave(filename=outfile_tmp, plot=ggobj, width = 7, height =7)
# # 可视化blueprint.main注释下，基因表达分布 -- FeaturePlot
# outfile_tmp <- paste0(res_home, "Figure/8.single_spatial/GSE189889.interest_genes.dimplot.pdf")
# ggobj <- FeaturePlot(seurat_obj_AM, features = interest_genes, cols=c("lightgray","orange", "red"), split.by="SampleType")
# ggsave(filename=outfile_tmp, plot=ggobj, width = 14, height =14)

tnsf9_grp <- rep("neg", dim(seurat_obj_AM_sub@meta.data)[1])
tnsf9_grp[which(seurat_obj_AM_sub@assays$RNA@data["TNFSF9", ]>0)] <- "pos"

seurat_obj_AM_sub@meta.data[, "tnsf9_grp"] <- tnsf9_grp

table(seurat_obj_AM_sub@meta.data[, c("tnsf9_grp", "anno_final")])

seurat_obj_AM_sub@meta.data[, "com2grp"] <- paste(seurat_obj_AM_sub@meta.data[, "anno_final"],
    seurat_obj_AM_sub@meta.data[, "tnsf9_grp"],sep=".")

Idents(seurat_obj_AM_sub) <- "com2grp"
diff_genes <- FindMarkers(object=seurat_obj_AM_sub, ident.1 = "Melanocytes.pos", ident.2 = "Melanocytes.neg")

diff_genes_ENTREZ <- id_convert(input=rownames(diff_genes), from="SYMBOL", to="ENTREZID")

# 功能富集
fun_res_GO <- clusterProfiler_enrich(inputgenes=diff_genes_ENTREZ$ENTREZID, min_size=5, cutoff=0.05, adjustMethod="BH", enrich_type="GO", simplify=FALSE)
fun_res_KEGG <- clusterProfiler_enrich(inputgenes=diff_genes_ENTREZ$ENTREZID, min_size=5, cutoff=0.05, adjustMethod="BH", enrich_type="KEGG", simplify=FALSE)
save(fun_res_GO, file=paste0(res_home, "RData/phase5.single_spatial/fun_res_GO.RData"))
save(fun_res_KEGG, file=paste0(res_home, "RData/phase5.single_spatial/fun_res_KEGG.RData"))




diff_gene_total <- FindMarkers(object=seurat_obj_AM_sub, ident.1 = "Melanocytes.pos", ident.2 = "Melanocytes.neg",
    logfc.threshold = 0)
save(diff_gene_total, file=paste0(res_home, "RData/phase5.single_spatial/diff_gene_total.RData"))


rank_gene_tmp <- diff_gene_total[, "avg_log2FC"]
names(rank_gene_tmp) <- rownames(diff_gene_total)
rank_gene_tmp <- sort(rank_gene_tmp, decreasing=TRUE)
gsea_res_immunesigdb <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="~/Data/Gene_Set/MsigDB/version7.5/c7.immunesigdb.v7.5.1.symbols.gmt")
gsea_res_biocarta <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="~/Data/Gene_Set/MsigDB/version7.5/c2.cp.biocarta.v7.5.1.symbols.gmt")
gsea_res_GO <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="~/Data/Gene_Set/MsigDB/version7.5/c5.go.bp.v7.5.1.symbols.gmt")
gsea_res_KEGG <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="~/Data/Gene_Set/MsigDB/version7.5/c2.cp.kegg.v7.5.1.symbols.gmt")
gsea_res_Reactome <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="~/Data/Gene_Set/MsigDB/version7.5/c2.cp.reactome.v7.5.1.symbols.gmt")


save(gsea_res_immunesigdb, file=paste0(res_home, "RData/phase5.single_spatial/gsea_res_immunesigdb.RData"))
save(gsea_res_biocarta, file=paste0(res_home, "RData/phase5.single_spatial/gsea_res_biocarta.RData"))
save(gsea_res_GO, file=paste0(res_home, "RData/phase5.single_spatial/gsea_res_GO.RData"))
save(gsea_res_KEGG, file=paste0(res_home, "RData/phase5.single_spatial/gsea_res_KEGG.RData"))
save(gsea_res_Reactome, file=paste0(res_home, "RData/phase5.single_spatial/gsea_res_Reactome.RData"))


# ex_idx <- c(2:3, 5:6, 9, 12:13, 15, 24, 35, 37:40, 45, 54:55, 84, 90:91, 102, 115, 159)
ex_idx <- c(83, 96, 103, 131)
#  [83] "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION"
#  [96] "REACTOME_TCR_SIGNALING"
# [103] "REACTOME_ADAPTIVE_IMMUNE_SYSTEM"
# [131] "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION"

gsea_res_sub <- gsea_res_Reactome[ex_idx, , asis=TRUE]
outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_reactome.tnfs9high.ridgeplot.pdf")
ggobj <- ridgeplot(gsea_res_sub, label_format=100, orderBy = "NES", fill = "NES", decreasing = TRUE)
ggobj <- ggobj + scale_fill_gradient(low = "#fca33c", high ="#d90017", space = "Lab", na.value = "grey50")
ggsave(plot=ggobj, filename=outfile_tmp, width=14)



gsea_res_sub <- gsea_res_GO[ex_idx, , asis=TRUE]


# gsea_res_Reactome2 <- gsea_res_Reactome[ex_idx, , asis=TRUE]
# gsea_res_Reactome2 <- gsea_res_Reactome2[grep("CELL_CYCLE", gsea_res_Reactome2[,"ID"], invert=TRUE), , asis=TRUE]

# gsea_res_Reactome2 <- gsea_res_Reactome2[grep("MAPK", gsea_res_Reactome2[,"ID"], invert=TRUE), , asis=TRUE]
# gsea_res_Reactome2 <- gsea_res_Reactome2[grep("REACTOME_SIGNALING_BY_NOTCH4", gsea_res_Reactome2[,"ID"], invert=TRUE), , asis=TRUE]
# # subset(gsea_res_Reactome, 1:dim(gsea_res_Reactome)[1] %in% ex_idx)
outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.pdf")
res <- clusterProfiler_plot(enrich_obj=gsea_res_sub, outfile=outfile_tmp, 
  showCategory=50, gseaplot=TRUE, geneSetID=1, child_width=7)






ex_idx <- c(grep("LYMPHOCYTE", gsea_res_GO[,1]),
  grep("IMMUNE", gsea_res_GO[,1]),
  grep("_T_CELL", gsea_res_GO[,1]))

ex_idx <- setdiff(ex_idx, grep("T_CELL_SELECTION", gsea_res_GO[,1]))






gsea_res_sub <- gsea_res_GO[ex_idx, , asis=TRUE]


# gsea_res_Reactome2 <- gsea_res_Reactome[ex_idx, , asis=TRUE]
# gsea_res_Reactome2 <- gsea_res_Reactome2[grep("CELL_CYCLE", gsea_res_Reactome2[,"ID"], invert=TRUE), , asis=TRUE]

# gsea_res_Reactome2 <- gsea_res_Reactome2[grep("MAPK", gsea_res_Reactome2[,"ID"], invert=TRUE), , asis=TRUE]
# gsea_res_Reactome2 <- gsea_res_Reactome2[grep("REACTOME_SIGNALING_BY_NOTCH4", gsea_res_Reactome2[,"ID"], invert=TRUE), , asis=TRUE]
# # subset(gsea_res_Reactome, 1:dim(gsea_res_Reactome)[1] %in% ex_idx)
outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.pdf")
res <- clusterProfiler_plot(enrich_obj=gsea_res_sub, outfile=outfile_tmp, 
  showCategory=50, gseaplot=TRUE, geneSetID=1, child_width=7)


data_barplot <- gsea_res_sub[, c("ID", "NES", "pvalue", "p.adjust")]

data_barplot[, "padj_grp"] <- "Sig."
ns_idx <- which(data_barplot[, "p.adjust"]>0.05)
if(length(ns_idx)>0)
data_barplot[ns_idx, "padj_grp"] <- "NS."

data_barplot[, "NES_grp"] <- "Enriched"
depleted_idx <- which(data_barplot[, "NES_grp"]<0)
if(length(depleted_idx)>0)
data_barplot[depleted_idx, "NES_grp"] <- "Depleted"


outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.barplot.pdf")
fill_space <- c("red", "blue")
names(fill_space) <- c("Enriched", "Depleted")
color_space <- c("black", "white")
names(color_space) <- c("Sig.", "NS.")
#' @author ZGX
ggobj <- barplot_ggplot(inputData_ggplot=data_barplot, outfile=outfile_tmp, xcol="ID", ycol="NES",
    fill_col='NES_grp', color_col="padj_grp", scale_fill_manual=fill_space, scale_fill_gradient = NULL, 
    mid_value=NULL, scale_colour_manual = color_space, logyaxis=NULL, 
    title=NULL, stack_barplot=FALSE, order=TRUE, add_text=FALSE, horizontal=TRUE, width = 12, height = 15)


  ggobj_list <- NULL
  # 计算emapplot
  edo <- pairwise_termsim(gsea_res_sub)
  ggobj_list[["star"]] <- emapplot(edo, layout="star", cex_line=0.6, cex_label_category=0.6)
  ggobj_list[["nicely"]] <- emapplot(edo, layout="nicely", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["circle"]] <- emapplot(edo, layout="circle", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["gem"]] <- emapplot(edo, layout="gem", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["dh"]] <- emapplot(edo, layout="dh", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["graphopt"]] <- emapplot(edo, layout="graphopt", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["grid"]] <- emapplot(edo, layout="grid", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["mds"]] <- emapplot(edo, layout="mds", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["randomly"]] <- emapplot(edo, layout="randomly", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["fr"]] <- emapplot(edo, layout="fr", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["kk"]] <- emapplot(edo, layout="kk", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["drl"]] <- emapplot(edo, layout="drl", cex_line=0.6, cex_label_category=0.6)
    ggobj_list[["lgl"]] <- emapplot(edo, layout="lgl", cex_line=0.6, cex_label_category=0.6)




outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.emapplot.pdf")
  multiplot_ggarrange(ggobj_list=ggobj_list, outfile=outfile_tmp, labels=names(ggobj_list), ncol=3, nrow=5, legend="bottom",
      common.legend=FALSE, width=15, height=25)


outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.ridgeplot.pdf")
ggobj <- ridgeplot(gsea_res_sub, label_format=100, orderBy = "NES", fill = "NES", decreasing = TRUE)
ggobj <- ggobj + scale_fill_gradient(low = "#fca33c", high ="#d90017", space = "Lab", na.value = "grey50")
ggsave(plot=ggobj, filename=outfile_tmp, width=14)


 


outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.treeplot.pdf")

ggobj_list <- NULL

ggobj_list[["ward.D"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="ward.D")
ggobj_list[["ward.D2"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="ward.D2")
ggobj_list[["average"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="average")
ggobj_list[["mcquitty"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="mcquitty")
ggobj_list[["median"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="median")
ggobj_list[["centroid"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="centroid")
ggobj_list[["complete"]] <- treeplot(edo, label_format=1000, nClusters=7, offset=20, offset_tiplab=1, clust_method="complete")

outfile_tmp <- paste0(res_home, "Figure/phase5.single_spatial/gsea_res_go.tnfs9high.treeplot.pdf")
  multiplot_ggarrange(ggobj_list=ggobj_list, outfile=outfile_tmp, labels=names(ggobj_list), ncol=2, nrow=4, legend="bottom",
      common.legend=FALSE, width=14, height=28)
ggsave(plot=ggobj_list[["ward.D2"]], filename=outfile_tmp, width=12)
