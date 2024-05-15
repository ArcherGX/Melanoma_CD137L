# ========================================================================================================
# Fig. S3 A
# ========================================================================================================
# 导入函数
libSources <- list.files("/home/pub/RCodes/CodeLib", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设置目录
res_home <- "/home/pub/Projects/Melanoma_AMPK/Res/"
checkDir(paste0(res_home, "Figure/4.DNAme_AM/"))
checkDir(paste0(res_home, "Table/4.DNAme_AM/"))
checkDir(paste0(res_home, "RData/4.DNAme_AM/"))

# ===========================================================================================================
# 注意这个是CA vs PBMC的，没啥意义
# ===========================================================================================================
# 探针信息
# /u02/Data/Melanoma_Patients/贺毅憬数据/晶能YP186分析结果/116个样品_20170609/results
table_path <- "/u02/Data/Melanoma_Patients/贺毅憬数据/晶能YP186分析结果/116个样品_20170609/results/mainAnno.All.CpG.xls"
mainAnno <- read.table(table_path, sep="\t", header=TRUE, quote="")

mainAnno_sub <- mainAnno[grep("TNFSF9", mainAnno[, "GencodeCompV12_NAME"]), ]
tss_idx <- which(startsWith(x=mainAnno_sub[, "GencodeCompV12_Group"], prefix="TSS"))
mainAnno_sub <- mainAnno_sub[tss_idx, ]
rownames(mainAnno_sub) <- mainAnno_sub[, "probeID"]
mainAnno_sub <- mainAnno_sub[, -c(1, 118:124)]
# mainAnno_sub_avg <- as.data.frame(colMeans(mainAnno_sub))
# mainAnno_sub_avg <- data.frame(sample=rownames(mainAnno_sub_avg), mainAnno_sub_avg)
# mainAnno_sub_avg[, "grp"] <- gsub(".*[0-9]", "", mainAnno_sub_avg[,"sample"])
# colnames(mainAnno_sub_avg)[2] <- "avgMethyl"
# # csv_path <- "/u02/Data/Melanoma_Patients/贺毅憬数据/晶能YP186分析结果/分析结果/合同号_GB20160805WP18-NG_客户名_贺毅憬_甲基化芯片_20170307/Report/1.Raw_Data_and_Pre-Process/methylated.Signal.csv"
# # data_exprs2 <- read.csv(csv_path, header=TRUE)
# outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/avgMethyl.TNFSF9.pdf")
# res <- boxplot_ggplot2(inputData_ggplot=mainAnno_sub_avg, outfile=outfile_tmp, x="grp", y="avgMethyl", color="grp", 
#     pval_group="grp", label="p.signif", force_plot=TRUE, width=7, height=7)


# # [1] "ProteinGroups_combat.rds.gz" "ProteinGroups_infor.rds.gz" 
# length(intersect(colnames(data_exprs), data_info$ID))
# sum(data_exprs$Gene.names=="CD137L")
# ===========================================================================================================
# 注意这个是CA vs CJ的
# ===========================================================================================================
table_path <- "/u02/Data/Melanoma_Patients/DNAme/彭聪/项目部3/甲基化数据/TableControl__hyj_81s_bg_20170303/TableControl__hyj_81s_bg_20170303.txt"
data_exprs <- read.table(table_path, sep="\t", header=TRUE)
data_info <- readRDS("/home/pub/Data/Melanoma_bulk_private/Melanoma_ICB_cohort/ProteinGroups_infor.rds.gz")

# 下载 Infinium HumanMethylation450K v1.2 Product Files
# HumanMethylation450 v1.2 Manifest File (CSV Format)
# https://support.illumina.com/array/array_kits/infinium_humanmethylation450_beadchip_kit/downloads.html
# 
data_exprs_sub <- data_exprs[which(data_exprs[, "TargetID"] %in% rownames(mainAnno_sub)), ]
data_exprs_sub <- data_exprs_sub[, grep("Pval", colnames(data_exprs_sub), invert=TRUE)]
data_exprs_sub <- data_exprs_sub[, -c(1:4)]


data_exprs_sub_sub_avg <- as.data.frame(colMeans(data_exprs_sub))
data_exprs_sub_sub_avg <- data.frame(sample=rownames(data_exprs_sub_sub_avg), data_exprs_sub_sub_avg)

grp <- gsub(".*_", "", gsub("\\.AVG_Beta", "", data_exprs_sub_sub_avg[,"sample"]))
grp <- toupper(gsub("[0-9]*", "", grp))
grp <- gsub("XY", "", gsub("\\.", "", grp))
grp <- gsub("CP", "CJ", grp)
data_exprs_sub_sub_avg[, "grp"] <- factor(grp, levels=c("CJ", "CA"))


colnames(data_exprs_sub_sub_avg)[2] <- "avgMethyl"
# csv_path <- "/u02/Data/Melanoma_Patients/贺毅憬数据/晶能YP186分析结果/分析结果/合同号_GB20160805WP18-NG_客户名_贺毅憬_甲基化芯片_20170307/Report/1.Raw_Data_and_Pre-Process/methylated.Signal.csv"
# data_exprs2 <- read.csv(csv_path, header=TRUE)
outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/avgMethyl_CACJ.TNFSF9.pdf")
res <- boxplot_ggplot2(inputData_ggplot=data_exprs_sub_sub_avg, outfile=outfile_tmp, x="grp", y="avgMethyl", color="grp", 
    pval_group="grp", label="p.signif", force_plot=TRUE, width=7, height=7)

outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/avgMethyl_CACJ.TNFSF9.ylim.pdf")
ggobj <- res$ggobj + ylim(0, 1)
ggsave(plot=ggobj, filename=outfile_tmp, width=3)

outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/avgMethyl_CACJ.TNFSF9.density_mean.pdf")
res <- densityPlot_ggplot(inputData=data_exprs_sub_sub_avg, outfile=outfile_tmp,
  x_col='avgMethyl', group_col='grp', add="mean", rug=TRUE, title=NULL, col_palette="lancet", width=7, height=7)
outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/avgMethyl_CACJ.TNFSF9.density_median.pdf")
res <- densityPlot_ggplot(inputData=data_exprs_sub_sub_avg, outfile=outfile_tmp,
  x_col='avgMethyl', group_col='grp', add="median", rug=TRUE, title=NULL, col_palette="lancet", width=7, height=7)
ggobj <- res + xlim(0, 1)
ggsave(plot=ggobj, filename=outfile_tmp)

    library(ggpubr)



methyl_grp_avg <- tapply(X=data_exprs_sub_sub_avg[, "avgMethyl"], INDEX=data_exprs_sub_sub_avg[, "grp"], FUN=mean, na.rm=TRUE)
methyl_grp_avg <- as.data.frame(methyl_grp_avg)
methyl_grp_avg <- data.frame(grp=rownames(methyl_grp_avg), methyl_grp_avg)

color_space <- getColor_palette(col_num=3)[1:2]
names(color_space) <- c("CA", "CJ")


outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/methyl_grp_avg_CACJ.TNFSF9.pdf")
ggobj <- barplot_ggplot(inputData_ggplot=methyl_grp_avg, outfile=NULL, xcol="grp", 
  ycol="methyl_grp_avg", fill_col="grp", color_col=NULL,
  scale_fill_manual=color_space, scale_fill_gradient = NULL, mid_value=NULL, scale_colour_manual = NULL, logyaxis=NULL,
  title=NULL, stack_barplot=FALSE, order=TRUE,
  add_text=FALSE, horizontal=FALSE, width = 7, height = 7)
ggsave(plot=ggobj, filename=outfile_tmp, width=4)
# 设置ylim
outfile_tmp <- paste0(res_home, "Figure/4.DNAme_AM/methyl_grp_avg_CACJ.TNFSF9.ylim.pdf")
ggobj <- ggobj+ylim(0, 1)
ggsave(plot=ggobj, filename=outfile_tmp, width=3)

methyl_grp_avg[1,2]-methyl_grp_avg[2,2]
# > methyl_grp_avg[1,2]-methyl_grp_avg[2,2]
# [1] 0.02303845