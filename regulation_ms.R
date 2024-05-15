# ========================================================================================================
# Fig. 3 B-D
# Fig. S3 D
# ========================================================================================================
# 导入函数
libSources <- list.files("~/RCodes/CodeLib", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设置目录
res_home <- "~/Projects/Melanoma_AMPK/Res/"
data_home <- "~/Projects/Melanoma_AMPK/Data/"

# 设定感兴趣的基因，注意首先检查基因名是否能够在TCGA数据中找到
interest_genes <- c("TNFSF9", "PRKAA1", "HLTF")
key_gene <- "TNFSF9"
interest_cancer_type <- "SKCM"

# 确保目录存在
checkDir(paste0(res_home, "fig3.regulation_ms/Figure/"))
checkDir(paste0(res_home, "fig3.regulation_ms/Table/"))
checkDir(paste0(res_home, "fig3.regulation_ms/RData/"))

inputURL <- paste0(data_home, "protein.html")
# 读取数据
library(rvest)
file <- read_html(inputURL)
htmlNodes <- html_nodes(file, "table")
proteins_ms <- html_table(htmlNodes[1], fill=TRUE)[[1]]
proteins_ms <- as.data.frame(proteins_ms)
write.csv(proteins_ms, paste0(res_home, "fig3.regulation_ms/Table/proteins_ms.csv"), quote=FALSE)

gene_ms <- gsub(" PE.*", "", gsub(".* GN=", "", proteins_ms[, "Description"]))
proteins_ms[, "gene"] <- gene_ms

gene_ms_df <- proteins_ms
gene_ms_df[gene_ms_df[, "gene"] %in% c("GTF2IRD2B", "ADNP", "HLTF"), ]


# 读取转录因子列表
human_tfs <- read.table("~/Data/BioAnno/Human_TFs/TF_names_v_1.01.txt")[,1]

# 导入数据 
load(file=paste0(data_home_TCGA, "exprs_score/data_exprs_list.RData"))
# 提取表达谱中存在的基因
com_genes <- intersect(union("TNFSF9", gene_ms), rownames(data_exprs_list[[1]]))
length(com_genes)/length(union("TNFSF9", gene_ms))
# > length(com_genes)/length(union("TNFSF9", gene_ms))
# [1] 0.9333333
data_exprs_sub_list <- lapply(data_exprs_list, function(x){
    t(x[com_genes, ])
})

# 计算基因间相关性
cor_res_list <- list()
for(cancer_type in names(data_exprs_sub_list)){
    data_exprs_sub <- data_exprs_sub_list[[cancer_type]]
    cor_res_list[[cancer_type]] <- cor_asymmetry(input1=data_exprs_sub[, "TNFSF9", drop=FALSE], input2=data_exprs_sub, cor_method="spearman", merge_df=TRUE)
}
save(cor_res_list, file=paste0(res_home, "fig3.regulation_ms/RData/cor_res_list.RData"))


# 提取转录因子的结果
cor_res_list_merge <- lapply(cor_res_list, function(input){
    # 提取合并矩阵
    cor_res_df <- input[["cor_pval_merge"]]
    cor_res_df[, "padj"] <- p.adjust(cor_res_df[, "pval"], method="BH")
    # 提取转录因子的结果
    cor_res_df_tf <- cor_res_df[which(cor_res_df[, "var2"] %in% c("HLTF", human_tfs)), ]
    # 添加注释列
    cor_res_df_tf[, "color_col"] <- "NEG"
    cor_res_df_tf[which(cor_res_df_tf[, "corr"]>0), "color_col"] <- "POS"
    cor_res_df_tf[, "sig"] <- "NS"
    cor_res_df_tf[which(cor_res_df_tf[, "padj"]<0.05), "sig"] <- "SIG"
    return(cor_res_df_tf)
})

# names(cor_res_list)
# > names(cor_res_list)
# [1] "correlation"    "pvalue"         "cor_pval_merge"



# 输出dot matrix stat
cor_res <- merge2df(cor_res_list_merge)
colnames(cor_res)[1] <- "cancer_type"
cor_res[, "var2"] <- as.character(cor_res[, "var2"])
cor_res[, "cancer_type"] <- gsub("TCGA-", "", cor_res[, "cancer_type"])
# 输出点阵图
outfile_tmp <- paste0(res_home, "fig3.regulation_ms/Figure/corr_TNFSF9.dot_matrix_stat.pdf")
res_list <- dot_matrix_stat(input = cor_res, outfile = outfile_tmp, x_col = "var2", y_col = "cancer_type", color_col = "corr", pval_col = "padj", title = NULL,
    color_space = c("blue", "white", "red"), color_mid = 0, order = TRUE, black_mark = TRUE, text_padding_x=1, text_padding_y=1)


# 输出散点图
cor_res_sub <- subset(cor_res, var2=="HLTF")
cor_res_sub[, "log10padj"] <- -log10(cor_res_sub[, "padj"])
cor_res_sub[, "abs_corr"] <- abs(cor_res_sub[, "corr"])
cor_res_sub[which(cor_res_sub[, "padj"]>0.05), "color_col"] <- "NS"
color_space_tmp <- c("blue", "lightgrey", "red")
names(color_space_tmp) <- c("NEG", "NS", "POS")
outfile_tmp <- paste0(res_home, "fig3.regulation_ms/Figure/corr_TNFSF9_HLTF.scatter.pdf")
ggobj_scatter <- scatter_ggplot_basic(input=cor_res_sub, outfile=outfile_tmp, xcol="corr", ycol="log10padj", color_col="color_col", color_pal=color_space_tmp, 
    text_col="cancer_type", size_col="abs_corr", title = NULL, ellipse = FALSE, mean.point = FALSE, vertical_line=NULL, horizontal_line=(-log10(0.05)))


cor_res_sub2 <- subset(cor_res_sub, cancer_type!="THCA")
outfile_tmp <- paste0(res_home, "fig3.regulation_ms/Figure/corr_TNFSF9_HLTF.scatter.sub.pdf")
ggobj_scatter <- scatter_ggplot_basic(input=cor_res_sub2, outfile=outfile_tmp, xcol="corr", ycol="log10padj", color_col="color_col", color_pal=color_space_tmp, 
    text_col="cancer_type", size_col="abs_corr", title = NULL, ellipse = FALSE, mean.point = FALSE, vertical_line=NULL, horizontal_line=(-log10(0.05)))



# 设定barplot颜色
fill_space_tmp <- c("blue", "red")
names(fill_space_tmp) <- c("NEG", "POS")
color_space_tmp <- c("black", "white")
names(color_space_tmp) <- c("SIG", "NS")
# 输出barplot
for(cancer_type in names(cor_res_list_merge)){
    cor_res_tmp <- cor_res_list_merge[[cancer_type]]
    outfile_tmp <- paste0(res_home, "fig3.regulation_ms/Figure/corr_TNFSF9.", cancer_type, ".barplot.pdf")
    res <- barplot_ggplot(inputData_ggplot=cor_res_tmp, outfile=outfile_tmp, xcol="var2", ycol="corr", fill_col="color_col", color_col="sig",
        scale_fill_manual=fill_space_tmp, scale_fill_gradient = NULL, mid_value=NULL, scale_colour_manual = color_space_tmp,
        title=NULL, stack_barplot=FALSE, order=TRUE,
        add_text=FALSE, horizontal=TRUE, width = 7, height = 5)
}

cancer_type <- "TCGA-SKCM"
cor_res_tmp <- cor_res_list_merge[[cancer_type]]
cor_res_tmp[, "log10padj"] <- (-log10(cor_res_tmp[, "padj"]))
cor_res_tmp[, "abs_corr"] <- abs(cor_res_tmp[, "corr"])
cor_res_tmp[which(cor_res_tmp[, "sig"]=="NS"), "color_col"] <- "NS"
color_space_tmp <- c("blue", "lightgrey", "red")
names(color_space_tmp) <- c("NEG", "NS", "POS")
outfile_tmp <- paste0(res_home, "fig3.regulation_ms/Figure/corr_TNFSF9.", cancer_type, ".scatter.pdf")
res <- scatter_ggplot_basic(input=cor_res_tmp, outfile=outfile_tmp, xcol="corr", ycol="log10padj",
    color_col="color_col", color_pal=color_space_tmp,
    text_col="var2", size_col="abs_corr", 
    vertical_line=NULL, horizontal_line=(-log10(0.05)))




