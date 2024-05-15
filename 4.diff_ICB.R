# ========================================================================================================
# Fig. 1 B
# ========================================================================================================
# 设置目录
checkDir(paste0(res_home, "Figure/phase5.surv_response/"))
checkDir(paste0(res_home, "Table/phase5.surv_response/")) 
checkDir(paste0(res_home, "RData/phase5.surv_response/"))

# 导入数据 
load(file=paste0(res_home, "RData/phase5.surv_response/data_response_merge_list_ICB.Data"))
load(file=paste0(res_home, "RData/phase5.surv_response/gsva_TNF_list_ICB.RData"))
# checkDir(paste0(res_home, "3.surv_res_ICB/Figure/response_diff"))
feats <- colnames(gsva_TNF_list_ICB[[1]])[-1]
res_list <- list()
for(dataset in names(data_response_merge_list_ICB)){
  data_response_merge <- data_response_merge_list_ICB[[dataset]]
  rownames(data_response_merge) <- data_response_merge[, "PID"]
  outfile_tmp <- paste0(res_home, "Figure/phase5.surv_response/response_diff.", dataset, ".pdf")
  com_genes <- intersect(feats, colnames(data_response_merge))
  # data_response_merge[, com_genes] <- log2(data_response_merge[, com_genes]+1)
  # data_response_merge[, "Response_infor"] <- gsub("1", "NR", gsub("0", "R", as.character(data_response_merge[, "Response_infor"])))
  res_list[[dataset]] <- boxplot_ggplot2_ggarrange(input=data_response_merge, outfile=outfile_tmp, features=com_genes, 
    pval_group="benefit_used", test_method="wilcox.test",
    label="p.signif", plot_add="jitter", col_palette="lancet", 
    ncol=3, nrow=2, width=7, height=9, width_child=4, height_child=4, type_child="boxplot")
}

# 针对RSF9进行小提琴可视化
# response_diff.GSE91061
# response_diff.GSE99070
# response_diff.GSE135222
res_list <- list()
dataset_selected <- c("GSE91061", "GSE99070", "GSE135222")
for(dataset in dataset_selected){
  data_response_merge <- data_response_merge_list_ICB[[dataset]]
  rownames(data_response_merge) <- data_response_merge[, "PID"]

  data_response_merge[, "benefit_used"] <- factor(as.character(data_response_merge[, "benefit_used"]), levels=c("NonBenefit", "Benefit"))

  outfile_tmp <- paste0(res_home, "Figure/phase5.surv_response/response_diff.", dataset, ".violin.pdf")
  res_list[[dataset]] <- violin_ggplot2(inputData_ggplot=data_response_merge, outfile=outfile_tmp, x="benefit_used", y="RSF9", color="benefit_used", 
    pval_group="benefit_used", test_method=c("wilcox.test", "kruskal.test"), label="p.signif",
    title=NULL, plot_add=c("boxplot", "median_iqr"), col_palette=c("lancet", "npg", "aaas", "jco"),
    force_plot=TRUE, width=5, height=7)
}

#


diff_sig_res <- lapply(res_list, function(x){
  x <- x[x[, "group2"]!="None", ]
  res <- data.frame(x, sig=as.numeric(x[, "p"]<=0.05))
  res <- res[, c(".y.", "sig")]
  colnames(res) <- c("feat", "sig")
  return(res)

})
diff_sig_res <- merge2df(diff_sig_res)
diff_sig_res <- unique(diff_sig_res)

diff_sig_res_matrix <- convArrType(inputArr=diff_sig_res, keyCols="feat", valueCols="sig", wide2long=FALSE)
rownames(diff_sig_res_matrix) <- diff_sig_res_matrix[, 1]
diff_sig_res_matrix <- diff_sig_res_matrix[, -1]

colSums(diff_sig_res_matrix)
# > colSums(diff_sig_res_matrix)
# RSF18 RSF1B RSF25  RSF4  RSF8  RSF9
#     0     6     0     4     3     3
res_RSF18 <- c(0, 0)
res_RSF1B <- c(1, 0)
res_RSF25 <- c(0, 0)
res_RSF4 <- c(1, 0)
res_RSF8 <- c(0, 0)
res_RSF9 <- c(3, 0)
# response_diff.GSE91061
# response_diff.GSE99070
# response_diff.GSE135222


diff_res_barplot <- matrix(c(res_RSF18, res_RSF1B, res_RSF25, res_RSF4, res_RSF8, res_RSF9), ncol=2, byrow=TRUE)
rownames(diff_res_barplot) <- c("RSF18", "RSF1B", "RSF25", "RSF4", "RSF8", "RSF9")
colnames(diff_res_barplot) <- c("protective", "risk")

diff_res_barplot <- data.frame(feat=rownames(diff_res_barplot), diff_res_barplot)
diff_res_barplot <- convArrType(diff_res_barplot, keyCols="feat", valueCols=c("protective", "risk"), wide2long=TRUE)

fill_col_tmp <- c("blue", "red")
names(fill_col_tmp) <- c("protective", "risk")
diff_res_barplot[,1] <- factor(diff_res_barplot[,1], levels=c("RSF9", "RSF4", "RSF1B", "RSF18", "RSF25", "RSF8"))
outfile_tmp <- paste0(res_home, "Figure/phase5.surv_response/response_diff_stat.barplot.pdf")
res <- barplot_ggplot(inputData_ggplot=diff_res_barplot, outfile=outfile_tmp, xcol="feat", ycol="value", fill_col="variable", color_col=NULL,
  scale_fill_manual=fill_col_tmp, scale_fill_gradient = NULL, mid_value=NULL, scale_colour_manual = NULL, logyaxis=NULL,
  title=NULL, stack_barplot=TRUE, order=FALSE,
  add_text=FALSE, horizontal=FALSE, width = 7, height = 7)


res_RSF18 <- c(0, 0)
res_RSF1B <- c(5, 1)
# response_diff.GSE35640, anti-MAGE-A3
# response_diff.GSE78220
# response_diff.GSE91061
# response_diff.GSE126044
# response_diff.phs000452.v3
# response_diff.PRJEB23709_ipiPD1

res_RSF25 <- c(0, 0)
res_RSF4 <- c(3, 1)
# response_diff.GSE78220
# response_diff.GSE91061
# response_diff.GSE126044
# response_diff.PRJEB23709_ipiPD1
res_RSF8 <- c(2, 1)
# response_diff.Braun_NatMed_2020
# response_diff.GSE126044
# response_diff.PRJEB23709_ipiPD1
res_RSF9 <- c(3, 0)
# response_diff.GSE91061
# response_diff.GSE99070
# response_diff.GSE135222

diff_res_barplot <- matrix(c(res_RSF18, res_RSF1B, res_RSF25, res_RSF4, res_RSF8, res_RSF9), ncol=2, byrow=TRUE)
rownames(diff_res_barplot) <- c("RSF18", "RSF1B", "RSF25", "RSF4", "RSF8", "RSF9")
colnames(diff_res_barplot) <- c("protective", "risk")

diff_res_barplot <- data.frame(feat=rownames(diff_res_barplot), diff_res_barplot)
diff_res_barplot <- convArrType(diff_res_barplot, keyCols="feat", valueCols=c("protective", "risk"), wide2long=TRUE)
diff_res_pie <- data.frame(diff_res_barplot, target="response")



sup_RSF18 <- data.frame("RSF18", "NS", 24, "response")
sup_RSF1B <- data.frame("RSF1B", "NS", 18, "response")
sup_RSF25 <- data.frame("RSF25", "NS", 24, "response")
sup_RSF4 <- data.frame("RSF4", "NS", 20, "response")
sup_RSF8 <- data.frame("RSF8", "NS", 21, "response")
sup_RSF9 <- data.frame("RSF9", "NS", 21, "response")
colnames(sup_RSF18) <- colnames(diff_res_pie)
colnames(sup_RSF1B) <- colnames(diff_res_pie)
colnames(sup_RSF25) <- colnames(diff_res_pie)
colnames(sup_RSF4) <- colnames(diff_res_pie)
colnames(sup_RSF8) <- colnames(diff_res_pie)
colnames(sup_RSF9) <- colnames(diff_res_pie)

diff_res_pie <- rbind(diff_res_pie, sup_RSF18)
diff_res_pie <- rbind(diff_res_pie, sup_RSF1B)
diff_res_pie <- rbind(diff_res_pie, sup_RSF25)
diff_res_pie <- rbind(diff_res_pie, sup_RSF4)
diff_res_pie <- rbind(diff_res_pie, sup_RSF8)
diff_res_pie <- rbind(diff_res_pie, sup_RSF9)

fill_col_tmp <- c("blue", "red", "lightgrey")
names(fill_col_tmp) <- c("protective", "risk", "NS")

outfile_tmp <- paste0(res_home, "Figure/phase5.surv_response/response_diff_stat.pieMatrix2.pdf")
pie_matrix_ggplot2(input=diff_res_pie, outfile=outfile_tmp, xcol="feat", ycol="target", 
  count_col="value", count_max=24, color_col="variable", color_space=fill_col_tmp,
  width = 7, height = 7)









# 采用“秦桧法”，，，
# IMvigor210	anti-PDL1	bladder	29443960
# GSE179351	anti-PD1_and_anti-CTLA4	colon cancer and pancreatic cancer	
# PRJEB25780	anti-PD1	gastric cancer	33547412
# GSE135222	anti-PD1_or_anti-PDL1	lung cancer	31537801
# GSE99070	anti-PD1	malignant pleural mesothelioma (MPM)	29618661
# GSE91061	anti-PD1	melanoma	29033130
# GSE176307	anti-PDL1	metastatic UC	34294892
# GSE67501	anti-PD1	renal cell carcinoma	27491898
# GSE165252	anti-PDL1_and_chemoradiotherapy	resectable esophageal adenocarcinoma	33504550
# 共涉及10种不同癌症类型
res_RSF18 <- c(0, 0)
res_RSF1B <- c(2, 0)
# response_diff.GSE91061
# response_diff.GSE126044
res_RSF25 <- c(0, 0)
res_RSF4 <- c(2, 0)
# response_diff.GSE91061
# response_diff.GSE126044
res_RSF8 <- c(1, 0)
# response_diff.GSE126044
res_RSF9 <- c(3, 0)
# response_diff.GSE91061
# response_diff.GSE99070
# response_diff.GSE135222
diff_res_barplot <- matrix(c(res_RSF18, res_RSF1B, res_RSF25, res_RSF4, res_RSF8, res_RSF9), ncol=2, byrow=TRUE)
rownames(diff_res_barplot) <- c("RSF18", "RSF1B", "RSF25", "RSF4", "RSF8", "RSF9")
colnames(diff_res_barplot) <- c("protective", "risk")
diff_res_barplot <- data.frame(feat=rownames(diff_res_barplot), diff_res_barplot)
diff_res_barplot <- convArrType(diff_res_barplot, keyCols="feat", valueCols=c("protective", "risk"), wide2long=TRUE)
diff_res_pie <- data.frame(diff_res_barplot, target="response")
sup_RSF18 <- data.frame("RSF18", "NS", 9, "response")
sup_RSF1B <- data.frame("RSF1B", "NS", 7, "response")
sup_RSF25 <- data.frame("RSF25", "NS", 9, "response")
sup_RSF4 <- data.frame("RSF4", "NS", 7, "response")
sup_RSF8 <- data.frame("RSF8", "NS", 8, "response")
sup_RSF9 <- data.frame("RSF9", "NS", 6, "response")
colnames(sup_RSF18) <- colnames(diff_res_pie)
colnames(sup_RSF1B) <- colnames(diff_res_pie)
colnames(sup_RSF25) <- colnames(diff_res_pie)
colnames(sup_RSF4) <- colnames(diff_res_pie)
colnames(sup_RSF8) <- colnames(diff_res_pie)
colnames(sup_RSF9) <- colnames(diff_res_pie)

diff_res_pie <- rbind(diff_res_pie, sup_RSF18)
diff_res_pie <- rbind(diff_res_pie, sup_RSF1B)
diff_res_pie <- rbind(diff_res_pie, sup_RSF25)
diff_res_pie <- rbind(diff_res_pie, sup_RSF4)
diff_res_pie <- rbind(diff_res_pie, sup_RSF8)
diff_res_pie <- rbind(diff_res_pie, sup_RSF9)

fill_col_tmp <- c("blue", "red", "lightgrey")
names(fill_col_tmp) <- c("protective", "risk", "NS")

outfile_tmp <- paste0(res_home, "Figure/phase5.surv_response/response_diff_stat.pieMatrix2.pdf")
pie_matrix_ggplot2(input=diff_res_pie, outfile=outfile_tmp, xcol="feat", ycol="target", 
  count_col="value", count_max=9, color_col="variable", color_space=fill_col_tmp,
  width = 7, height = 7)






