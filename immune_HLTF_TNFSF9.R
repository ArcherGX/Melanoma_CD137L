# ========================================================================================================
# Fig. 3 M
# Fig. S3 F
# ========================================================================================================
# 导入函数
libSources <- list.files("~/RCodes/CodeLib", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设置目录
res_home <- "~/Projects/Melanoma_AMPK/Res/"
data_home_spec <- "~/Projects/Melanoma_AMPK/Data/"
data_home_independent <- "~/Data/Melanoma_bulk_public/Melanoma_GEO/"
data_home_TCGA <- "~/Data/TCGA_secondary/"
data_home_ICB_protein <- "~/Data/Melanoma_bulk_public/Melanoma_Immunotherapy_Protein_2021/"
data_home_ICB <- "~/Data/data_ICB_secondary/"

# 确保目录存在
checkDir(paste0(res_home, "Figure/phase3.immune_HLTF_TNFSF9/"))
checkDir(paste0(res_home, "Table/phase3.immune_HLTF_TNFSF9/"))
checkDir(paste0(res_home, "RData/phase3.immune_HLTF_TNFSF9/"))

# 导入数据 
load(file=paste0(data_home_TCGA, "exprs_score/data_exprs_list.RData"))
load(file=paste0(data_home_TCGA, "exprs_score/data_ImmuCellAI_list.RData"))

data_exprs_tmp <- t(data_exprs_list[["TCGA-SKCM"]])
data_immune_tmp <- data_ImmuCellAI_list[["TCGA-SKCM"]]


TNFSF9_mid <- median(data_exprs_tmp[, "TNFSF9"])
HLTF_mid <- median(data_exprs_tmp[, "HLTF"])

TNFSF9_grp <- rep("High", length(data_exprs_tmp[, "TNFSF9"]))
TNFSF9_grp[data_exprs_tmp[, "TNFSF9"]<TNFSF9_mid] <- "Low"

HLTF_grp <- rep("High", length(data_exprs_tmp[, "HLTF"]))
HLTF_grp[data_exprs_tmp[, "HLTF"]<HLTF_mid] <- "Low"

grp_merge <- paste0(paste0("TNFSF9_", TNFSF9_grp), "_&_", paste0("HLTF_", HLTF_grp))
grp_merge <- data.frame(ID=rownames(data_exprs_tmp), grp_merge=grp_merge)
grp_immune_merge <- merge(grp_merge, data_immune_tmp, by="ID")
grp_immune_merge[, "grp_merge2"] <- grp_immune_merge[, "grp_merge"]

grp_immune_merge[, "grp_merge2"] <- gsub("TNFSF9_High_&_HLTF_High", "Mix", grp_immune_merge[, "grp_merge2"])
grp_immune_merge[, "grp_merge2"] <- gsub("TNFSF9_Low_&_HLTF_Low", "Mix", grp_immune_merge[, "grp_merge2"])


level_order <- c("TNFSF9_High_&_HLTF_Low", "TNFSF9_High_&_HLTF_High", "TNFSF9_Low_&_HLTF_Low", "TNFSF9_Low_&_HLTF_High")
grp_immune_merge[, "grp_merge"] <- factor(grp_immune_merge[, "grp_merge"], levels=level_order)

level_order <- c("TNFSF9_High_&_HLTF_Low", "Mix", "TNFSF9_Low_&_HLTF_High")
grp_immune_merge[, "grp_merge2"] <- factor(grp_immune_merge[, "grp_merge2"], levels=level_order)

table(grp_immune_merge[, "grp_merge"])
# > table(grp_immune_merge[, "grp_merge"])
# TNFSF9_High_&_HLTF_High  TNFSF9_High_&_HLTF_Low  TNFSF9_Low_&_HLTF_High
#                     103                     131                     131
#   TNFSF9_Low_&_HLTF_Low
#                     103

feats_tmp <- c("CD8.T", "NK", "Macrophage", "Tc", "DC", "Tex", "Tem", "CD4.T", "B.cell", "nTreg")
outfile_tmp <- paste0(res_home, "Figure/phase3.immune_HLTF_TNFSF9/grp_immune_merge.pdf")
res <- boxplot_ggplot2_ggarrange(input=grp_immune_merge, outfile=outfile_tmp, features=feats_tmp, pval_group="grp_merge", test_method="wilcox.test",
    label="p.signif", plot_add="jitter", col_palette="Set3", ncol=5, nrow=2, width=16, height=12, width_child=4, height_child=5, type_child="boxplot")

outfile_tmp <- paste0(res_home, "Figure/phase3.immune_HLTF_TNFSF9/grp_immune_merge2.pdf")
res <- boxplot_ggplot2_ggarrange(input=grp_immune_merge, outfile=outfile_tmp, features=feats_tmp, pval_group="grp_merge2", test_method="wilcox.test",
    label="p.signif", plot_add="jitter", col_palette="Set3", ncol=5, nrow=2, width=16, height=12, width_child=4, height_child=5, type_child="boxplot")
 
color

grp_immune_merge_sub <- grp_immune_merge[grep("06A$", grp_immune_merge[, "ID"]), ]
outfile_tmp <- paste0(res_home, "Figure/phase3.immune_HLTF_TNFSF9/grp_immune_merge_sub.pdf")
res <- boxplot_ggplot2_ggarrange(input=grp_immune_merge_sub, outfile=outfile_tmp, features=feats_tmp, pval_group="grp_merge", test_method="wilcox.test",
    label="p.signif", plot_add="jitter", col_palette="Set3", ncol=5, nrow=2, width=16, height=12, width_child=4, height_child=5, type_child="boxplot")


outfile_tmp <- paste0(res_home, "Figure/phase3.immune_HLTF_TNFSF9/grp_immune_merge_sub2.pdf")
res <- boxplot_ggplot2_ggarrange(input=grp_immune_merge_sub, outfile=outfile_tmp, features=feats_tmp, pval_group="grp_merge2", test_method="wilcox.test",
    label="p.signif", plot_add="jitter", col_palette="lancet", ncol=3, nrow=2, width=16, height=12, width_child=4, height_child=5, type_child="boxplot")









 <- surv_res_list_gene <- list()
for(dataset in names(data_merge_list_ICB)){
  data_clinical_merged <- data_merge_list_ICB[[dataset]]
  # 分组 -- TNFSF9
  pat_cat_list <- separatePatient(input=data_exprs_tmp, var_cols="TNFSF9", time="PFS_Time_use", event="PFS_use", 
    labels=c("Low", "High"), minprop=0.1, maxstat=FALSE, quantile_cutoff=0.5)
  TNFSF9_grp <- pat_cat_list$Categorize_Arr
  # 分组 -- CD274
  pat_cat_list <- separatePatient(input=data_clinical_merged, var_cols="CD274", time="PFS_Time_use", event="PFS_use", 
    labels=c("Low", "High"), minprop=0.1, maxstat=FALSE, quantile_cutoff=0.5)
  CD274_grp <- pat_cat_list$Categorize_Arr
  # 合并分组
  data_surv_grp <- TNFSF9_grp
  data_surv_grp[, "CD274"] <- CD274_grp[, "CD274"]
  data_surv_grp[, "merge_grp"] <- paste0(paste0("TNFSF9_", data_surv_grp[, "TNFSF9"]), "_AND_", paste0("CD274_", data_surv_grp[, "CD274"]))
  data_surv_grp[, "merge_grp"] <- as.factor(data_surv_grp[, "merge_grp"])
  surv_res_list_gene[[dataset]] <- SurvAnalysis_univar(inputArr=data_surv_grp, out_prefix=NULL, variable_col="merge_grp", time_col="PFS_Time_use", status_col="PFS_use", 
    force_plot=TRUE, getHR=FALSE, plotHR=FALSE, cutoff=0.05, conf.int=FALSE, risk.table=TRUE, title = dataset, legend_position="none", break.time.by=NULL, width=7, height=9)
}






surv_res_list_gene <- list()
for(dataset in names(data_merge_list_ICB)){
  data_clinical_merged <- data_merge_list_ICB[[dataset]]
  # 分组 -- TNFSF9
  pat_cat_list <- separatePatient(input=data_clinical_merged, var_cols="TNFSF9", time="PFS_Time_use", event="PFS_use", 
    labels=c("Low", "High"), minprop=0.1, maxstat=FALSE, quantile_cutoff=0.5)
  TNFSF9_grp <- pat_cat_list$Categorize_Arr
  # 分组 -- CD274
  pat_cat_list <- separatePatient(input=data_clinical_merged, var_cols="CD274", time="PFS_Time_use", event="PFS_use", 
    labels=c("Low", "High"), minprop=0.1, maxstat=FALSE, quantile_cutoff=0.5)
  CD274_grp <- pat_cat_list$Categorize_Arr
  # 合并分组
  data_surv_grp <- TNFSF9_grp
  data_surv_grp[, "CD274"] <- CD274_grp[, "CD274"]
  data_surv_grp[, "merge_grp"] <- paste0(paste0("TNFSF9_", data_surv_grp[, "TNFSF9"]), "_AND_", paste0("CD274_", data_surv_grp[, "CD274"]))
  data_surv_grp[, "merge_grp"] <- as.factor(data_surv_grp[, "merge_grp"])
  surv_res_list_gene[[dataset]] <- SurvAnalysis_univar(inputArr=data_surv_grp, out_prefix=NULL, variable_col="merge_grp", time_col="PFS_Time_use", status_col="PFS_use", 
    force_plot=TRUE, getHR=FALSE, plotHR=FALSE, cutoff=0.05, conf.int=FALSE, risk.table=TRUE, title = dataset, legend_position="none", break.time.by=NULL, width=7, height=9)
}



# 计算相关性
cor_res_list <- list()
for(dataset in names(melanoma_exprs_list)){
  data_matrix_tmp <- t(melanoma_exprs_list[[dataset]])
  if(sum(colnames(data_matrix_tmp) %in% c("TNFSF9", "HLTF"))==2){
    cor_data_list <- cor_asymmetry(input1 = data_matrix_tmp[, "TNFSF9", drop=FALSE], input2 = data_matrix_tmp[, "HLTF", drop=FALSE], cor_method = "spearman", merge_df = TRUE)
    cor_res_list[[dataset]] <- cor_data_list[["cor_pval_merge"]]
  }
}
cor_res_merge <- merge2df(cor_res_list)
cor_res_merge
# > cor_res_merge
#         .id   var1 var2        corr        pval
# 1  GSE65904 TNFSF9 HLTF -0.04865507 0.478631878
# 2 GSE131521 TNFSF9 HLTF -0.37009804 0.144268640
# 3  GSE19234 TNFSF9 HLTF  0.39393939 0.008553275
# 4  GSE59455 TNFSF9 HLTF -0.02037758 0.810235697
# 5 GSE100797 TNFSF9 HLTF  0.30307692 0.140735763
# 6  GSE78220 TNFSF9 HLTF -0.03448276 0.861844337
# 7 GSE133713 TNFSF9 HLTF  0.07707510 0.726170241
# 8  GSE22153 TNFSF9 HLTF  0.09184599 0.495783465
# 9  GSE22154 TNFSF9 HLTF -0.37210615 0.088899978
cor_res_merge[, "log10pval"] <- -log10(cor_res_merge[, "pval"])
outfile_tmp <- paste0(res_home, "Figure/phase3.corr_HLTF/public_melanoma_cohort.pdf")
ggobj <- scatter_ggplot_basic(input=cor_res_merge, outfile=outfile_tmp, xcol="corr", ycol="log10pval", color_col="black", color_pal=NULL, 
  text_col=".id", size_col=NULL, title = NULL, ellipse = FALSE, mean.point = FALSE, vertical_line=NULL, horizontal_line=-log10(0.05))


# =============================================================================================================================
# TCGA
# =============================================================================================================================
# 导入数据
load(file=paste0(data_home_TCGA, "exprs_score/data_exprs_list.RData"))
# 计算相关性
cor_res_list <- list()
for(dataset in names(data_exprs_list)){
  data_matrix_tmp <- t(data_exprs_list[[dataset]])
  if(sum(colnames(data_matrix_tmp) %in% c("TNFSF9", "HLTF"))==2){
    cor_data_list <- cor_asymmetry(input1 = data_matrix_tmp[, "TNFSF9", drop=FALSE], input2 = data_matrix_tmp[, "HLTF", drop=FALSE], cor_method = "spearman", merge_df = TRUE)
    cor_res_list[[dataset]] <- cor_data_list[["cor_pval_merge"]]
  }
}
cor_res_merge <- merge2df(cor_res_list)
cor_res_merge
# > cor_res_merge
#          .id   var1 var2        corr         pval
# 1   TCGA-ACC TNFSF9 HLTF -0.02468354 8.287368e-01
# 2  TCGA-BLCA TNFSF9 HLTF -0.11665592 1.889820e-02
# 3  TCGA-BRCA TNFSF9 HLTF -0.20865981 5.862209e-12
# 4  TCGA-CESC TNFSF9 HLTF -0.18092468 1.801884e-03
# 5  TCGA-CHOL TNFSF9 HLTF -0.20411840 2.315169e-01
# 6  TCGA-COAD TNFSF9 HLTF -0.21186683 5.678889e-06
# 7  TCGA-DLBC TNFSF9 HLTF -0.26295097 7.434513e-02
# 8  TCGA-ESCA TNFSF9 HLTF -0.23616892 3.471243e-03
# 9   TCGA-GBM TNFSF9 HLTF -0.19886665 1.700071e-02
# 10 TCGA-HNSC TNFSF9 HLTF -0.15833314 4.123822e-04
# 11 TCGA-KICH TNFSF9 HLTF -0.23831117 5.591621e-02
# 12 TCGA-KIRC TNFSF9 HLTF -0.21509764 6.779463e-07
# 13 TCGA-KIRP TNFSF9 HLTF -0.30175617 2.149940e-07
# 14 TCGA-LAML TNFSF9 HLTF -0.24692524 3.840665e-03
# 15  TCGA-LGG TNFSF9 HLTF -0.25599261 7.216208e-09
# 16 TCGA-LIHC TNFSF9 HLTF -0.02717244 6.028611e-01
# 17 TCGA-LUAD TNFSF9 HLTF -0.10457851 1.819259e-02
# 18 TCGA-LUSC TNFSF9 HLTF -0.09567119 3.318752e-02
# 19 TCGA-MESO TNFSF9 HLTF -0.06551491 5.604510e-01
# 20   TCGA-OV TNFSF9 HLTF -0.15872409 2.746058e-03
# 21 TCGA-PAAD TNFSF9 HLTF -0.24223308 1.195463e-03
# 22 TCGA-PCPG TNFSF9 HLTF -0.17445404 2.094322e-02
# 23 TCGA-PRAD TNFSF9 HLTF -0.21530752 1.881204e-06
# 24 TCGA-READ TNFSF9 HLTF -0.09609512 2.221221e-01
# 25 TCGA-SARC TNFSF9 HLTF -0.01639957 7.930674e-01
# 26 TCGA-SKCM TNFSF9 HLTF -0.19302890 2.711941e-05
# 27 TCGA-STAD TNFSF9 HLTF -0.18903232 2.468366e-04
# 28 TCGA-TGCT TNFSF9 HLTF -0.28809722 3.868643e-04
# 29 TCGA-THCA TNFSF9 HLTF -0.49496206 4.475608e-32
# 30 TCGA-THYM TNFSF9 HLTF -0.01061815 9.086442e-01
# 31 TCGA-UCEC TNFSF9 HLTF -0.05827447 1.774702e-01
# 32  TCGA-UCS TNFSF9 HLTF -0.17286398 2.020863e-01
# 33  TCGA-UVM TNFSF9 HLTF -0.22963878 4.474099e-02
cor_res_merge[, "log10pval"] <- -log10(cor_res_merge[, "pval"])

cor_res_merge[, "label"] <- "NS"
cor_res_merge[cor_res_merge[, "pval"]<0.05, "label"] <- "NEG"
cor_res_merge[, "text_col"] <- cor_res_merge[, ".id"]
cor_res_merge[cor_res_merge[, "pval"]>0.05, "text_col"] <- ""
color_space_tmp <- c("blue", "grey")
names(color_space_tmp) <- c("NEG", "NS")

outfile_tmp <- paste0(res_home, "Figure/phase3.corr_HLTF/TCGA_cohorts.pdf")
ggobj <- scatter_ggplot_basic(input=cor_res_merge, outfile=outfile_tmp, xcol="corr", ycol="log10pval", color_col="label", color_pal=color_space_tmp, 
  text_col="text_col", size_col=NULL, title = NULL, ellipse = FALSE, mean.point = FALSE, vertical_line=0, horizontal_line=-log10(0.05))



# =============================================================================================================================
# in house
# =============================================================================================================================
# 读取数据
Melanoma_inhouse_clinical_RNA_seq <- readRDS("~/Data/Melanoma_bulk_private/Melanoma_Patients_Processed_LiXin/Melanoma_inhouse_clinical_RNA_seq.rds")
data_exprs_hfc <- get(load("~/Data/Melanoma_bulk_private/Melanoma_HFCbatch_cohort/RNA_seq_FPKM140.RData"))
rownames(data_exprs_hfc) <- data_exprs_hfc[, 1]
data_exprs_hfc <- data_exprs_hfc[, -1]
# 去除全是0的行
data_exprs_hfc <- data_exprs_hfc[rowSums(data_exprs_hfc)>0, ]
data_exprs_list <- Melanoma_inhouse_clinical_RNA_seq[["Expr_FPKM_log2_Combat"]]
names(data_exprs_list) <- Melanoma_inhouse_clinical_RNA_seq$Cohort
data_exprs_list[["HFC_merge"]] <- log2(data_exprs_hfc+1)
sapply(data_exprs_list, dim)
# > sapply(data_exprs_list, dim)
#      treatment_naive immune_treatment merge_cohort HFC_merge
# [1,]           31994            31994        31994     56250
# [2,]             176              145          321       141
data_exprs_list <- lapply(data_exprs_list, function(x){
  x[, grep("CA$", colnames(x))]
})
sapply(data_exprs_list, dim)
# >  sapply(data_exprs_list, dim)
#      treatment_naive immune_treatment merge_cohort HFC_merge
# [1,]           31994            31994        31994     56250
# [2,]              99               69          168        50
# 计算相关性
cor_res_list <- list()
for(dataset in names(data_exprs_list)){
  data_matrix_tmp <- t(data_exprs_list[[dataset]])
  if(sum(colnames(data_matrix_tmp) %in% c("TNFSF9", "HLTF"))==2){
    cor_data_list <- cor_asymmetry(input1 = data_matrix_tmp[, "TNFSF9", drop=FALSE], input2 = data_matrix_tmp[, "HLTF", drop=FALSE], cor_method = "pearson", merge_df = TRUE)
    cor_res_list[[dataset]] <- cor_data_list[["cor_pval_merge"]]
  }
}
cor_res_merge <- merge2df(cor_res_list)
cor_res_merge
# > cor_res_merge
#                .id   var1 var2       corr         pval
# 1  treatment_naive TNFSF9 HLTF -0.1883109 0.0619559650
# 2 immune_treatment TNFSF9 HLTF -0.2156274 0.0751724329
# 3     merge_cohort TNFSF9 HLTF -0.1998584 0.0093938213
# 4        HFC_merge TNFSF9 HLTF -0.5022369 0.0002021085
cor_res_merge[, "log10pval"] <- -log10(cor_res_merge[, "pval"])
cor_res_merge[, "label"] <- "NS"
cor_res_merge[cor_res_merge[, "pval"]<0.05, "label"] <- "NEG"
color_space_tmp <- c("blue", "grey")
names(color_space_tmp) <- c("NEG", "NS")
# 输出
outfile_tmp <- paste0(res_home, "Figure/phase3.corr_HLTF/inhouse_cohorts.pdf")
ggobj <- scatter_ggplot_basic(input=cor_res_merge, outfile=outfile_tmp, xcol="corr", ycol="log10pval", color_col="label", color_pal=color_space_tmp, 
  text_col=".id", size_col=NULL, title = NULL, ellipse = FALSE, mean.point = FALSE, vertical_line=0, horizontal_line=-log10(0.05))
# 输出显著数据集的散点图
datasets_sig <- cor_res_merge[cor_res_merge[, "pval"]<0.05, ".id"]
for(dataset in datasets_sig){
  data_matrix_tmp <- t(data_exprs_list[[dataset]])[, c("TNFSF9", "HLTF")]
  outfile_tmp <- paste0(res_home, "Figure/phase3.corr_HLTF/inhouse_cohorts.", dataset, ".pdf")
  ggobj <- scatter_ggplot_adv(inputArr=data_matrix_tmp, outfile=outfile_tmp, xcol="HLTF", ycol="TNFSF9", title = dataset, 
    color="black", color_manual=NULL, col_pal_name="lancet", add="reg.line", facet.by=NULL, Marginal_Density=FALSE)
}



# ==============================================================================================
# 黑素瘤新辅助 -- 单细胞数据
# ==============================================================================================
# 导入数据
seurat_obj <- readRDS("/lustre/home/liuq/AM_Project/new_Res/0.Maintypes/RData/seurat_obj.rds.gz")
# 输出全部细胞
exprSet <- t(as.matrix(seurat_obj@assays[["RNA"]]@data))
outfile_tmp <- paste0(res_home, "Figure/phase3.corr_HLTF/neoadjuvant_single_cell_data.all.pdf")
ggobj <- scatter_ggplot_adv(inputArr=exprSet, outfile=outfile_tmp, xcol="HLTF", ycol="TNFSF9", title = NULL, 
  color="black", color_manual=NULL, col_pal_name="lancet", add="reg.line", facet.by=NULL, Marginal_Density=FALSE)

# 只输出肿瘤细胞
seurat_obj_sub <- subset(seurat_obj, MainTypes=="MelanomaCells")
exprSet <- t(as.matrix(seurat_obj_sub@assays[["RNA"]]@data))
outfile_tmp <- paste0(res_home, "Figure/phase3.corr_HLTF/neoadjuvant_single_cell_data.mel.pdf")
ggobj <- scatter_ggplot_adv(inputArr=exprSet, outfile=outfile_tmp, xcol="HLTF", ycol="TNFSF9", title = NULL, 
  color="black", color_manual=NULL, col_pal_name="lancet", add="reg.line", facet.by=NULL, Marginal_Density=FALSE)




