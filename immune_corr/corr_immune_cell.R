# ========================================================================================================
# Fig. 2 F
# ========================================================================================================
# 导入函数
libSources <- list.files("~/RCodes/CodeLib", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设置目录
res_home <- "~/Projects/Melanoma_AMPK/Res/"
data_home_TCGA <- "~/Data/TCGA_secondary/"
data_home_spec <- "~/Projects/Melanoma_AMPK/Data/"


# 确保目录存在
checkDir(paste0(res_home, "Figure/phase3.corr_immune/"))
checkDir(paste0(res_home, "Table/phase3.corr_immune/"))
checkDir(paste0(res_home, "RData/phase3.corr_immune/"))

# 导入数据 
load(file=paste0(data_home_TCGA, "exprs_score/data_exprs_list.RData"))
load(file=paste0(data_home_TCGA, "exprs_score/data_ImmuCellAI_list.RData"))

ex_cols <- c("ID", "CD4.T", "CD8.T", "Tc", "Tex", "nTreg", "Th1", "Th2", "Th17", "Tfh", 
    "DC", "NKT", "B.cell", "Monocyte", "Macrophage", "NK", "Neutrophil")
data_ImmuCellAI_list <- lapply(data_ImmuCellAI_list, function(input){
    input <- input[, ex_cols]
})

com_cancer_types <- intersect(names(data_exprs_list), names(data_ImmuCellAI_list))
exprs_immune_merge_list <- list()
for(cancer_type in com_cancer_types){
    data_exprs_tmp <- t(data_exprs_list[[cancer_type]])
    data_exprs_tmp <- data.frame(ID=rownames(data_exprs_tmp), data_exprs_tmp)
    data_exprs_tmp <- data_exprs_tmp[, c("ID", interest_genes)]
    data_ImmuCellAI_tmp <- data_ImmuCellAI_list[[cancer_type]]
    exprs_immune_merge_list[[cancer_type]] <- merge(data_exprs_tmp, data_ImmuCellAI_tmp)
}

immune_feat_cols <- colnames(data_ImmuCellAI_tmp)[-1]

# 计算基因间相关性
cor_res_list <- list()
for(cancer_type in names(exprs_immune_merge_list)){
    exprs_immune_merge <- exprs_immune_merge_list[[cancer_type]]
    cor_res_list[[cancer_type]] <- cor_asymmetry(input1=exprs_immune_merge[, interest_genes], input2=exprs_immune_merge[, immune_feat_cols], cor_method="spearman", merge_df=TRUE)
}
# 提取转录因子的结果
cor_res_list_merge <- lapply(cor_res_list, function(input){
    # 提取合并矩阵
    cor_res_df <- input[["cor_pval_merge"]]
    cor_res_df[, "padj"] <- p.adjust(cor_res_df[, "pval"], method="BH")
    # 添加注释列
    cor_res_df[, "color_col"] <- "NEG"
    cor_res_df[which(cor_res_df[, "corr"]>0), "color_col"] <- "POS"
    cor_res_df[, "sig"] <- "NS"
    cor_res_df[which(cor_res_df[, "padj"]<0.05), "sig"] <- "SIG"
    return(cor_res_df)
})
save(cor_res_list_merge, file=paste0(res_home, "RData/fig4.corr_immune/cor_res_list_merge.RData"))

# 输出dot matrix stat
cor_res <- merge2df(cor_res_list_merge)
colnames(cor_res)[1] <- "cancer_type"
cor_res[, "var2"] <- as.character(cor_res[, "var2"])
cor_res_list_gene <- split(cor_res, f=cor_res[, "var1"])

# 输出点阵图
for(gene in names(cor_res_list_gene)){
    cor_res_tmp <- cor_res_list_gene[[gene]]
    outfile_tmp <- paste0(res_home, "Figure/phase3.corr_immune/corr_immuneCellAI.", gene, ".dot_matrix_stat.pdf")
    res_list <- dot_matrix_stat(input = cor_res_tmp, outfile = outfile_tmp, x_col = "var2", y_col = "cancer_type", color_col = "corr", pval_col = "padj", title = gene,
        color_space = c("blue", "white", "red"), color_mid = 0, order = TRUE, black_mark = TRUE, text_padding_x=1, text_padding_y=1)
}




# 输出dot matrix stat
cor_res <- merge2df(cor_res_list_merge)
colnames(cor_res)[1] <- "cancer_type"
cor_res[, "var2"] <- as.character(cor_res[, "var2"])

# 输出dot matrix stat
cor_res <- merge2df(cor_res_list_merge)
colnames(cor_res)[1] <- "cancer_type"
# 输出散点图
cor_res_sub <- subset(cor_res, cancer_type=="TCGA-SKCM")
cor_res_sub[, "log10padj"] <- -log10(cor_res_sub[, "padj"])
cor_res_sub[, "abs_corr"] <- abs(cor_res_sub[, "corr"])
cor_res_sub[which(cor_res_sub[, "padj"]>0.05), "color_col"] <- "NS"
color_space_tmp <- c("blue", "lightgrey", "red")
names(color_space_tmp) <- c("NEG", "NS", "POS")

cor_res_sub_list <- split(cor_res_sub, f=cor_res_sub[, "var1"])
ggobj_list <- list()
# for(gene in names(cor_res_sub_list)){
for(gene in interest_genes){
ggobj_list[[gene]] <- scatter_ggplot_basic(input=cor_res_sub_list[[gene]], outfile=NULL, 
    xcol="corr", ycol="log10padj", color_col="color_col", color_pal=color_space_tmp, 
    text_col="var2", size_col="abs_corr", title = gene, 
    ellipse = FALSE, mean.point = FALSE, vertical_line=NULL, horizontal_line=(-log10(0.05)))
}

outfile_tmp <- paste0(res_home, "Figure/phase3.corr_immune/corr_immune_genes.scatter.pdf")
multiplot_ggarrange(ggobj_list=ggobj_list, outfile=outfile_tmp, labels="AUTO", ncol=3, nrow=1, 
    legend="bottom", common.legend=FALSE, width=15, height=5)


outfile_tmp <- paste0(res_home, "Figure/fig4.corr_immune/corr_immune_genes_HLTF.barplot.pdf")
cor_res_sub_tmp <- cor_res_sub_list[["HLTF"]]
cor_res_sub_tmp <- subset(cor_res_sub_tmp, var2!="Tex")


library(ggpubr)
ggobj <- ggbarplot(cor_res_sub_tmp, x="var2", y="corr",
fill="color_col",#changefillcolorbympg_level
color="sig",#Setbarbordercolorstowhite
palette="jco",#jcojournalcolorpalett.see?ggpar
sort.val="desc",#Sortthevalueindescendingorder
sort.by.groups=FALSE,#Don'tsortinsideeachgroup
x.text.angle=90,#Rotateverticallyxaxistexts
rotate=TRUE,
ggtheme=theme_light()
)

ggsave(plot=ggobj, filename=outfile_tmp)








