# ========================================================================================================
# Fig. 4 L
# ========================================================================================================
# 设置目录
checkDir(paste0(res_home, "RData/fig6a.diff_gene_mouse_AMPK/"))
checkDir(paste0(res_home, "Figure/fig6a.diff_gene_mouse_AMPK/"))
checkDir(paste0(res_home, "Table/fig6a.diff_gene_mouse_AMPK/"))

data_spec_home <- "/home/pub/Projects/Melanoma_AMPK/Data/"

# 导入数据 
data_diff <- read.csv(paste0(data_spec_home, "AICAR_ShExp_DEseq2_Result.csv"))


data_fc <- data_diff[which(data_diff[,1] %in% ex_genes), ]
data_fc[, "grp"] <- rep("UP", dim(data_fc)[1])
data_fc[which(data_fc[, "log2FoldChange"]<0), "grp"] <- "DN"
data_fc[which(data_fc[, "padj"]>0.05), "grp"] <- "NS"
color_space_tmp <- c("blue", "lightgrey", "red")
names(color_space_tmp) <- c("DN", "NS", "UP")


# data_diff2 <- data_diff
colnames(data_fc)[c(1,9,13)] <- c("SYMBOL", "logFC", "adj.P.Val")
outfile_tmp <- paste0(res_home, "fig6a.diff_gene_mouse_AMPK/Figure/diff_volcano_better.pdf")
res <- diff_volcano_plot(input=data_fc, output=outfile_tmp, log2fc_cutoff=1 , padj_cutoff=0.05, 
  genes_highlight=ex_genes, width=7, height=7)
  

# data_diff2 <- data_diff
# colnames(data_diff2)[c(9,13)] <- c("logFC", "adj.P.Val")
# res <- volcano_plot(input=data_diff2, output=outfile_tmp, log2fc_cutoff=1 , padj_cutoff=0.05, 
#   genes_highlight=NULL, width=7, height=7)
  
data_exprs <- read.csv(paste0(data_spec_home, "AMPK and Adora1 manipulated samples ResultsGeneExpression.csv"))

# C10orf54==VSIR 
ex_genes <- "TNFSF9,TNFRSF12A,C10ORF54,CEACAM1,PVR,TNFRSF9,CD274,PDCD1,CD86,TNFRSF4,HLA-DRB1,IL2RB,ADORA2A,PDCD1LG2,TNFRSF18,CD276,TNFSF4,CD200,LAG3,CD27,LGALS3,CTLA4,CD40,TNFRSF14,CD28"
# ex_genes <- "TNFSF9,TNFRSF12A,VSIR,CEACAM1,PVR,TNFRSF9,CD274,PDCD1,CD86,TNFRSF4,HLA-DRB1,IL2RB,ADORA2A,PDCD1LG2,TNFRSF18,CD276,TNFSF4,CD200,LAG3,CD27,LGALS3,CTLA4,CD40,TNFRSF14,CD28"
ex_genes <- unique(unlist(strsplit(ex_genes, split=",")))
# setdiff(ex_genes, data_exprs[,1])
# data_exprs[grep("^C10", data_exprs[,1]), 1]


data_fc <- data_exprs[which(data_exprs[,1] %in% ex_genes), c("symbol", "AMPK_FC", "AMPK_pval")]
data_fc[, "grp"] <- rep("UP", dim(data_fc)[1])
data_fc[which(data_fc[, "AMPK_FC"]<0), "grp"] <- "DN"
data_fc[which(data_fc[, "AMPK_pval"]>0.05), "grp"] <- "NS"
color_space_tmp <- c("blue", "lightgrey", "red")
names(color_space_tmp) <- c("DN", "NS", "UP")

outfile_tmp <- paste0(res_home, "fig6a.diff_gene_mouse_AMPK/Figure/data_exprs_sub.barplot.pdf")
res <- barplot_ggplot(inputData_ggplot=data_fc, outfile=outfile_tmp, xcol="symbol", ycol="AMPK_FC",
  fill_col="grp", color_col=NULL,
  scale_fill_manual=color_space_tmp, scale_fill_gradient = NULL, mid_value=0, scale_colour_manual = NULL, logyaxis=NULL,
  title=NULL, stack_barplot=FALSE, order=TRUE,
  add_text=FALSE, horizontal=FALSE, width = 7, height = 7)


# data_diff2 <- data_diff
colnames(data_fc)[c(1:3)] <- c("SYMBOL", "logFC", "adj.P.Val")
outfile_tmp <- paste0(res_home, "fig6a.diff_gene_mouse_AMPK/Figure/diff_volcano_better.pdf")
res <- diff_volcano_plot(input=data_fc, output=outfile_tmp, log2fc_cutoff=1 , padj_cutoff=0.05, 
  genes_highlight=ex_genes, width=7, height=7)
  


data_exprs_sub <- data_exprs[which(data_exprs[,1] %in% ex_genes),]
rownames(data_exprs_sub) <- data_exprs_sub[,1]
data_exprs_sub2 <- data_exprs_sub[order(data_exprs_sub[, "AMPK_FC"], decreasing=TRUE), c(2:7)]

InAcMarker_TCGAoverlap <- unique(InAcMarker_TCGAoverlap[, c("Symbol", "Role.with.Immunity")])
InAcMarker_TCGAoverlap_sup <- data.frame(Symbol=c("C10ORF54", "IL2RB"), Role.with.Immunity="Inhibit")
InAcMarker_TCGAoverlap <- rbind(InAcMarker_TCGAoverlap, InAcMarker_TCGAoverlap_sup)
rownames(InAcMarker_TCGAoverlap) <- InAcMarker_TCGAoverlap$Symbol
InAcMarker_TCGAoverlap <- InAcMarker_TCGAoverlap[rownames(data_exprs_sub2), ]
InAcMarker_TCGAoverlap <- InAcMarker_TCGAoverlap[, -1, drop=FALSE]
colnames(InAcMarker_TCGAoverlap) <- "Role"


 
color_role <- getColor_palette(pal_name="Set1", col_num=3)
names(color_role) <- unique(InAcMarker_TCGAoverlap[, "Role"])
my_colour = list(Role = color_role)

# anno_col <- data.frame(grp=c(rep("AICAR", 3), rep("CTRL", 3)))
# anno_col[, "grp"] <- factor(anno_col[, "grp"], levels=c("CTRL", "AICAR"))
# rownames(anno_col) <- colnames(data_exprs_sub2)
# data_heatmap <- getData_ComplexHeatmap(name="FPKM", inputArr=data_exprs_sub2,
#   heatmap_cols=c("blue", "white", "red"), anno_dis=FALSE, anno_cols_pal=NULL, annoType=NULL, pch=NA)
#   res_list <- list(name=name, arr=inputArr, annoType=annoType)

# data_anno <- getData_ComplexHeatmap(name="FPKM", inputArr=anno_col,
#   heatmap_cols=NULL, anno_dis=TRUE, anno_cols_pal="Set1", annoType="column", pch=NA)
#   res_list <- list(name=name, arr=inputArr, annoType=annoType)

# InAcMarker_TCGAoverlap

# # 输出
# getComplexHeatmap(outfile=outfile, heatmap_data_list=heatmap_data_list,
#   anno_dis_list=anno_dis_list, anno_barplot_list=NULL, anno_barplot_stack_list=NULL, anno_heatmap_list=NULL,
#   column_order=NULL, column_split=NULL, cluster_rows=TRUE, cluster_columns=TRUE, clust_method=distance, 
#   width=11, height=7) 

#   #' @author ZGX
# getComplexHeatmap <- function(outfile=NULL, heatmap_data_list=NULL,
#   anno_dis_list=NULL, anno_barplot_list=NULL, anno_barplot_stack_list=NULL, anno_right_list=NULL, anno_heatmap_list=NULL,
#   column_order=NULL, row_order=NULL, column_split=NULL, row_split=NULL, cluster_rows=TRUE, cluster_columns=FALSE, clust_method="euclidean",
#   vertical=TRUE, show_row_names=TRUE, show_column_names=FALSE, return_heatmap_obj=FALSE, title=NULL,
#   width=7, height=7){


library(pheatmap)
outfile_tmp <- paste0(res_home, "Figure/fig6a.diff_gene_mouse_AMPK/data_exprs_sub.pheatmap.pdf")
pdf(outfile_tmp)
pheatmap(data_exprs_sub2, cluster_rows = FALSE, annotation_row=InAcMarker_TCGAoverlap, annotation_colors = my_colour)
dev.off()




diff_genes <- data_diff[, "log2FoldChange"]
names(diff_genes) <- data_diff[, "gene_id"]
diff_genes <- sort(na.omit(diff_genes), decreasing=TRUE)
# 功能富集
gsea_res_immunesigdb <- clusterProfiler_enrich(inputgenes=diff_genes, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c7.immunesigdb.v7.5.1.symbols.gmt")
gsea_res_GO <- clusterProfiler_enrich(inputgenes=diff_genes, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c5.go.bp.v7.5.1.symbols.gmt")
gsea_res_KEGG <- clusterProfiler_enrich(inputgenes=diff_genes, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c2.cp.kegg.v7.5.1.symbols.gmt")
gsea_res_Reactome <- clusterProfiler_enrich(inputgenes=diff_genes, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c2.cp.reactome.v7.5.1.symbols.gmt")


outfile_tmp <- paste0(res_home, "3.clustering_annotation_sub/Figure/cell_anno_T/fun_res_Reactome/fun_res.", cluster_tmp, ".pdf")
res <- clusterProfiler_plot(enrich_obj=fun_res_Reactome_tmp, outfile=outfile_tmp, 
  showCategory=50, gseaplot=FALSE, geneSetID=1, child_width=7)


# 合并全部结果进行输出 -- GO富集
fun_res_GO <- merge2df(fun_res_GO_list)
fun_res_GO <- fun_res_GO[, c(".id", "Description", "pvalue", "p.adjust", "Count")]
fun_res_GO <- unique(fun_res_GO)
outfile_tmp <- paste0(res_home, "3.clustering_annotation_sub/Figure/cell_anno_T/fun_res_GO/fun_res_GO_merged.pdf")
res <- dot_matrix_ggplot2(input = fun_res_GO, outfile = outfile_tmp, 
  x_col = ".id", y_col = "Description", color_col = "Count", pval_col = "p.adjust", title = NULL, 
  color_space = c("white", "red"), color_mid = 0, order = FALSE, black_mark = TRUE,
  point_shape = "circle", better_theme=TRUE, width = 21, height = 95)
# 合并全部结果进行输出 -- KEGG富集
fun_res_KEGG <- merge2df(fun_res_KEGG_list)
fun_res_KEGG <- fun_res_KEGG[, c(".id", "Description", "pvalue", "p.adjust", "Count")]
fun_res_KEGG <- unique(fun_res_KEGG)
outfile_tmp <- paste0(res_home, "3.clustering_annotation_sub/Figure/cell_anno_T/fun_res_KEGG/fun_res_KEGG_merged.pdf")
res <- dot_matrix_ggplot2(input = fun_res_KEGG, outfile = outfile_tmp, 
  x_col = ".id", y_col = "Description", color_col = "Count", pval_col = "p.adjust", title = NULL, 
  color_space = c("white", "red"), color_mid = 0, order = FALSE, black_mark = TRUE,
  point_shape = "circle", better_theme=TRUE, width = 17, height = 15)
