# ========================================================================================================
# Fig. 1 J
# ========================================================================================================
# 导入函数
libSources <- list.files("~/RCodes/CodeLib", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设置目录
res_home <- "~/Projects/Melanoma_AMPK/Res/"
data_home <- "~/Projects/Melanoma_AMPK/Data/"
# 设置目录
checkDir(paste0(res_home, "Figure/9.melanoma_acral_mIHC/"))
checkDir(paste0(res_home, "Table/9.melanoma_acral_mIHC/"))
checkDir(paste0(res_home, "RData/9.melanoma_acral_mIHC/"))

# 读取数据
load(file=paste0(res_home, "RData/9.melanoma_acral_mIHC/data_clinical_exprs.RData"))

# ==========================================================================================================
# 基于maxstat分组
# ==========================================================================================================
sep_CD137L <- separatePatient(input=data_clinical_exprs, var_cols="CD137L", time="PFS_month", event="PFS", labels=c("Low", "High"), minprop=0.35, maxstat=TRUE, quantile_cutoff=0.5)
sep_CD137L$quantile_loc
# > sep_CD137L$quantile_loc
#      ecdf      sort 
# 0.5909091 0.5652174
sep_CD137L <- sep_CD137L$Categorize_Arr
sep_CD137L <- data.frame(PID=rownames(sep_CD137L), sep_CD137L)
sep_CD137L <- sep_CD137L[, c("PID", "CD137L")]
sep_PDL1 <- separatePatient(input=data_clinical_exprs, var_cols="PDL1", time="PFS_month", event="PFS", labels=c("Low", "High"), minprop=0.35, maxstat=TRUE, quantile_cutoff=0.5)
sep_PDL1$quantile_loc
# > sep_PDL1$quantile_loc
#      ecdf      sort
# 0.5454545 0.5217391
sep_PDL1 <- sep_PDL1$Categorize_Arr
sep_PDL1 <- data.frame(PID=rownames(sep_PDL1), sep_PDL1)
sep_PDL1 <- sep_PDL1[, c("PID", "PDL1")]
# 合并分组
sep_merge <- merge(sep_CD137L, sep_PDL1, by="PID")
colnames(sep_merge) <- c("PID", "CD137L_grp", "PDL1_grp")
data_merge <- merge(data_clinical_exprs, sep_merge, by="PID")
merge_grp <- paste(paste0("CD137L.", data_merge[, "CD137L_grp"]),
  paste0("PDL1.", data_merge[, "PDL1_grp"]), sep="_")
data_merge[, "merge_grp"] <- as.factor(merge_grp)
data_merge[, "merge_grp2"] <- as.factor(gsub("CD137L.Low_PDL1.High", "Others", gsub("CD137L.High_PDL1.Low", "Others", merge_grp)))

table(data_merge[, "merge_grp"])
# > table(data_merge[, "merge_grp"])
# CD137L.High_PDL1.High  CD137L.High_PDL1.Low  CD137L.Low_PDL1.High
#                    10                    1                     5
#   CD137L.Low_PDL1.Low
#                    13
table(data_merge[, "merge_grp2"])
# > table(data_merge[, "merge_grp2"])
# CD137L.High_PDL1.High   CD137L.Low_PDL1.Low                Others
#                    10                    13                     6
surv_res_list <- list()
data_merge[, "CD137L_grp"] <- factor(as.character(data_merge[, "CD137L_grp"]), levels=c("High", "Low"))
data_merge[, "PDL1_grp"] <- factor(as.character(data_merge[, "PDL1_grp"]), levels=c("High", "Low"))
for(grp_col in c("CD137L_grp", "PDL1_grp", "merge_grp", "merge_grp2")){
  surv_res_list[[grp_col]] <- SurvAnalysis_univar(inputArr=data_merge, out_prefix=NULL, variable_col=grp_col, 
    time_col="PFS_month", status_col="PFS", force_plot=TRUE, cutoff=0.05, conf.int=FALSE, risk.table=TRUE,
    title = grp_col, legend_position="none", break.time.by=NULL, width=7, height=9)
}
surv_plot_list <- lapply(surv_res_list, function(x){
  x$plot
})
outfile_tmp <- paste0(res_home, "Figure/9.melanoma_acral_mIHC/surv.KM.maxstat.part1.pdf")
res <- arrange_surv_plots(ggsurvplot_list=surv_plot_list[1:2], outfile=outfile_tmp, ncol = 2, nrow = 1, width=14, height=9)
outfile_tmp <- paste0(res_home, "Figure/9.melanoma_acral_mIHC/surv.KM.maxstat.part2.pdf")
res <- arrange_surv_plots(ggsurvplot_list=surv_plot_list[4], outfile=outfile_tmp, ncol = 1, nrow = 1, width=8.5, height=9)

outfile_tmp <- paste0(res_home, "Figure/9.melanoma_acral_mIHC/surv.KM.maxstat.part2.2.pdf")
res <- arrange_surv_plots(ggsurvplot_list=surv_plot_list[3], outfile=outfile_tmp, ncol = 1, nrow = 1, width=8.5, height=9)

# ggsave(plot=surv_plot_list[[3]], filename=outfile_tmp, width=9, height=9)



grp_col = "PDL1_grp"
response_col = "Response"


chisq_table <- table(data_merge[, c(response_col, grp_col)])
chisq_res <- getChisqTest(testTable=chisq_table)
  
data_barplot <- as.data.frame(chisq_table)
grp_total <- tapply(X=data_barplot[, "Freq"], INDEX=data_barplot[, grp_col], FUN=sum)
data_barplot[, "ratio"] <- data_barplot[,"Freq"]/grp_total[data_barplot[, grp_col]]
outfile_tmp <- paste0(res_home, "Figure/9.melanoma_acral_mIHC/PDL1_grp.response.barplot.stack.pdf")

fill_colors <- getColor_palette(pal_custom=NULL, pal_name="Set1", col_num=2, pal_set="brewer", outfile=NULL, force=TRUE, print_info=FALSE)
names(fill_colors) <- c("R", "NR")
res <- barplot_ggplot(inputData_ggplot=data_barplot, outfile=outfile_tmp, xcol="PDL1_grp", ycol='ratio', fill_col="Response",
  color_col=NULL, scale_fill_manual=fill_colors, scale_fill_gradient = NULL, mid_value=NULL, scale_colour_manual = NULL, logyaxis=NULL,
  title=NULL, stack_barplot=TRUE, order=FALSE,
  add_text=FALSE, horizontal=FALSE, width = 4, height = 7)


