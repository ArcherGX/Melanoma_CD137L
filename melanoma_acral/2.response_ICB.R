# ========================================================================================================
# Fig. 1 F
# ========================================================================================================
# 设置目录
checkDir(paste0(res_home, "Figure/8.melanoma_acral_ICB/"))
checkDir(paste0(res_home, "Table/8.melanoma_acral_ICB/"))
checkDir(paste0(res_home, "RData/8.melanoma_acral_ICB/"))

# 读取数据
Melanoma_inhouse_clinical_RNA_seq <- readRDS("~/Data/Melanoma_bulk_private/Melanoma_Patients_Processed_LiXin/Melanoma_inhouse_clinical_RNA_seq.rds")
data_exprs <- Melanoma_inhouse_clinical_RNA_seq$Expr_FPKM_log2_Combat[[2]]
data_clinical <- Melanoma_inhouse_clinical_RNA_seq$ClinicalData[[2]]

ex_cols <- c("Sample_use", "Response_use", "PFS", "PFST.Month", "OS", "OST.Month", "Subtypes")
data_clinical <- data_clinical[, ex_cols]
data_clinical <- subset(data_clinical, Subtypes=="Acral")
 

# colnames(data_clinical) <- c("pid", "date_diagnosis", "alive_dead", "os.date", "OS.time", "OS.status")
data_exprs <- data_exprs[, grep("CA$", colnames(data_exprs))]
data_exprs <- t(data_exprs[interest_genes, ,drop=FALSE])

data_exprs <- data.frame(Sample_use=rownames(data_exprs), data_exprs)
data_exprs[,"Sample_use"] <- gsub(".CA", "", data_exprs[,"Sample_use"])


data_merge <- merge(data_clinical, data_exprs, by="Sample_use")


outfile_tmp <- paste0(res_home, "Figure/8.melanoma_acral_ICB/diff_response.pdf")
res <- boxplot_ggplot2_ggarrange(input=data_merge, outfile=outfile_tmp, features=interest_genes, pval_group="Response_use", 
  paired = FALSE, test_method="wilcox.test", label="p.signif", plot_add="jitter", col_palette="lancet", 
  ncol=3, nrow=1, width=10, height=7, width_child=4, height_child=4, type_child="boxplot")
    # 对每一个特征进行计算


data_response_merge <- data_merge
data_response_merge[,"Response_use"] <- as.numeric(gsub("R", "0", gsub("NR", "1", data_response_merge[,"Response_use"])))
reg_res_list <- list()
for(interest_gene in interest_genes){
  reg_res_list[[interest_gene]] <- reg_core(inputArr=data_response_merge, x_cols=interest_gene, y_cols="Response_use", reg_family="binomial")
}
reg_res_multi <- reg_core(inputArr=data_response_merge, x_cols=interest_genes, y_cols="Response_use", reg_family="binomial")


feat_coef2 <- reg_res_multi$feat_coef_long
feat_coef2 <- feat_coef2[-1,]
outfile_tmp <- paste0(res_home, "Figure/8.melanoma_acral_ICB/feat_coef_response.pdf")
res <- barplot_ggplot(inputData_ggplot=feat_coef2, outfile=outfile_tmp, xcol="obj", ycol="coef", 
  fill_col=NULL, color_col=NULL,
  scale_fill_manual=NULL, scale_fill_gradient = NULL, mid_value=NULL, scale_colour_manual = NULL, logyaxis=NULL,
  title=NULL, stack_barplot=FALSE, order=TRUE,
  add_text=FALSE, horizontal=FALSE, width = 4, height = 7)
