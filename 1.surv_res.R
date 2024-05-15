# ========================================================================================================
# Fig. 1 A
# ========================================================================================================
# 设置目录
checkDir(paste0(res_home, "Figure/phase5.surv_response/"))
checkDir(paste0(res_home, "Table/phase5.surv_response/")) 
checkDir(paste0(res_home, "RData/phase5.surv_response/"))

# 导入数据 
load(file=paste0(res_home, "RData/phase5.hallmark_immune/gsva_TNF_list.RData"))
load(file=paste0(data_home_TCGA, "exprs_score/data_clinical.RData"))

# 提取数据
data_exprs_sub_list <- list()
for(cancer_type in names(gsva_TNF_list)){
	# 提取感兴趣基因的表达谱
	data_exprs_tmp <- gsva_TNF_list[[cancer_type]]
	data_exprs_tmp <- t(data_exprs_tmp)
	data_exprs_tmp <- data.frame(bcr_patient_barcode=rownames(data_exprs_tmp), data_exprs_tmp)
	if(cancer_type=="TCGA-SKCM"){
		data_exprs_tmp <- data_exprs_tmp[grep("06A$", data_exprs_tmp[, "bcr_patient_barcode"]), ]
	}
	# pat_id <- substring(rownames(data_exprs_tmp), first=1, last=12)
	# dup_pat <- pat_id[which(duplicated(pat_id))]
	# ex_idx <- which(!(pat_id %in% dup_pat))
	# data_exprs_tmp <- data_exprs_tmp[ex_idx, ]
	data_exprs_tmp[, "bcr_patient_barcode"] <- substring(data_exprs_tmp[, "bcr_patient_barcode"], first=1, last=12)
	data_exprs_sub_list[[cancer_type]] <- data_exprs_tmp
}
# 设定临床终点
outcome_list_tmp <- list(OS=c("OS", "OS.time"),
	PFI=c("PFI", "PFI.time"))

# 循环每一个基因进行计算
feats <- colnames(data_exprs_sub_list[[cancer_type]])[-1]
surv_res_list_all <- list()
for(gene in feats){
	surv_res_list_all[[gene]] <- survTCGA_single_factor(data_exprs_list=data_exprs_sub_list, data_clinical=data_clinical, interest_gene=gene, 
		outcome_list=outcome_list_tmp, quantile_grp3=NULL)
}
save(surv_res_list_all, file=paste0(res_home, "RData/phase5.surv_response/surv_res_list_all.RData"))

# 可视化输出全部KM结果
surv_res_table_list <- list()
for(gene in interest_genes){
	outdir_tmp <- paste0(res_home, "3.surv_res_TCGA/Figure/surv_res_TCGA.", gene, "_grp2")
	surv_res_table_list2 <- survTCGA_output_res(surv_res_list=surv_res_list_all[[gene]], outdir=outdir_tmp, outcome_list=c("OS", "DSS", "DFI", "PFI"))
	surv_res_table_list[[gene]] <- surv_res_table_list2
}
save(surv_res_table_list, file=paste0(res_home, "3.surv_res_TCGA/RData/surv_res_table_list.RData"))


hr_res <- merge2df(lapply(surv_res_table_list, function(x){
	x$hr_res
}))
midtime_res <- merge2df(lapply(surv_res_table_list, function(x){
	x$midtime_res
}))
colnames(hr_res)[1] <- "gene"
colnames(midtime_res)[1] <- "gene"
# 输出表格 -- grp2
writeTable(midtime_res, paste0(res_home, "3.surv_res_TCGA/Table/midtime_res.grp2.txt"))
writeTable(hr_res, paste0(res_home, "3.surv_res_TCGA/Table/hr_res.grp2.txt"))


# 可视化输出HR点阵 -- grp2
hr_res_list <- split(hr_res, f=hr_res[, "outcome"])
for(outcome in names(hr_res_list)){
	outfile_tmp <- paste0(res_home, "3.surv_res_TCGA/Figure/interest_genes.grp2.dot_matrix.", outcome, ".pdf")
	hr_res_tmp <- hr_res_list[[outcome]]
	ggobj <- dot_matrix_ggplot2(input = hr_res_tmp, outfile = outfile_tmp, 
		x_col = "cancer_type", y_col = "gene", color_col = "HR", pval_col = "Likelihood_pvalue", title = NULL, 
		color_space = c("blue", "white", "red"), color_mid = 1, order = TRUE, black_mark = TRUE,
		point_shape = "circle", better_theme=TRUE, width = 7, height = 3)
	# 设定颜色空间
	color_space_tmp <- c("blue", "red")
	names(color_space_tmp) <- c("protect", "risk")
	# 设定分组
	risk_factor <- rep("risk", dim(hr_res_tmp)[1])
	risk_factor[which(hr_res_tmp[, "HR"]<1)] <- "protect"
	hr_res_tmp[, "risk_factor"] <- risk_factor
	hr_res_tmp[, "log10Pvalue"] <- -log10(hr_res_tmp[, "Likelihood_pvalue"])
	outfile_tmp <- paste0(res_home, "3.surv_res_TCGA/Figure/interest_genes.grp2.dot_matrix_disc.", outcome, ".pdf")
	ggobj <- point_disc_ggplot(inputData_ggplot = hr_res_tmp, outfile = outfile_tmp, 
		x_col = "cancer_type", y_col = "gene", pval_col = "Likelihood_pvalue", color_col = "risk_factor", 
		color_space = color_space_tmp, clustering = TRUE, 
		sig_theme = "fading", point_shape = "circle", legend = FALSE, width = 9, height = 2.5)
}





