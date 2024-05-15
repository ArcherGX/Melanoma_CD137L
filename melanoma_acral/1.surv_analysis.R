# ========================================================================================================
# Fig. 1 F
# ========================================================================================================
# 设置目录
checkDir(paste0(res_home, "Figure/7.melanoma_acral/"))
checkDir(paste0(res_home, "Table/7.melanoma_acral/"))
checkDir(paste0(res_home, "RData/7.melanoma_acral/"))

# 读取数据
Melanoma_inhouse_clinical_RNA_seq <- readRDS("~/Data/Melanoma_bulk_private/Melanoma_Patients_Processed_LiXin/Melanoma_inhouse_clinical_RNA_seq.rds")
data_exprs <- Melanoma_inhouse_clinical_RNA_seq$Expr_FPKM_log2_Combat[[1]]
data_clinical <- Melanoma_inhouse_clinical_RNA_seq$ClinicalData[[1]]

ex_cols <- c("Sample_use", "OS", "OST.Month", "Subtypes")
data_clinical <- data_clinical[, ex_cols]
data_clinical <- subset(data_clinical, Subtypes=="Acral")
data_clinica[, "OST.Month"]

# colnames(data_clinical) <- c("pid", "date_diagnosis", "alive_dead", "os.date", "OS.time", "OS.status")
data_exprs <- data_exprs[, grep("CA$", colnames(data_exprs))]
data_exprs <- t(data_exprs[interest_genes, ,drop=FALSE])

data_exprs <- data.frame(Sample_use=rownames(data_exprs), data_exprs)
data_exprs[,"Sample_use"] <- gsub(".CA", "", data_exprs[,"Sample_use"])


data_merge <- merge(data_clinical, data_exprs, by="Sample_use")

# as.Date(data_clinical[, "Date.of.Diagnosis"])

# # ==================================================================================================
# # 计算预后 -- 基于离散分组
# # ==================================================================================================
# 设定临床终点
outcome_list <- list(OS=c("OS", "OST.Month"))

surv_res_list <- list()
for(outcome_type in names(outcome_list)){
  time_tmp <- outcome_list[[outcome_type]][2]
  event_tmp <- outcome_list[[outcome_type]][1]
  surv_res_list_tmp <- list()
  for(interest_gene in interest_genes){
    pat_cat_list <- separatePatient(input=data_merge, var_cols=interest_gene, time=time_tmp, event=event_tmp, 
      labels=c("Low", "High"), minprop=0.25, maxstat=TRUE, quantile_cutoff=0.5)
    data_clinical_merged2 <- pat_cat_list$Categorize_Arr
    data_clinical_merged2[, interest_gene] <- factor(as.character(data_clinical_merged2[, interest_gene]), levels=c("High", "Low"))
    # out_prefix_tmp <- paste0(res_home, "Figure/7.melanoma_acral/", interest_gene)
    surv_res_list_tmp[[interest_gene]] <- SurvAnalysis_univar(inputArr=data_clinical_merged2, out_prefix=NULL, 
      variable_col=interest_gene, time_col=time_tmp, status_col=event_tmp, force_plot=TRUE, cutoff=0.05, 
      conf.int=FALSE, risk.table=TRUE, title=interest_gene, legend_position="none", break.time.by=NULL, 
      width=7, height=9)
  }
  surv_res_list[[outcome_type]] <- surv_res_list_tmp
}

# # 提取surv plot
surv_list_os <- lapply(surv_res_list[["OS"]], function(x){x$plot})
# 输出两个基因的结果到一张图上
outfile_tmp <- paste0(res_home, "Figure/7.melanoma_acral/surv_os.KM_plot.maxstat.pdf")
arrange_surv_plots(ggsurvplot_list=surv_list_os, outfile=outfile_tmp, ncol =3, nrow = 1, width=18, height=8)

