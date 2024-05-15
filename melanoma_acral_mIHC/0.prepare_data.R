# ========================================================================================================
# Fig. 1 J, prepare data
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
data_clinical <- read_xls2(paste0(data_home, "免疫治疗 号子+预后 kxw new-1.1.xlsx"))
colnames(data_clinical)[1:5] <- c("PID", "Gender", "Age", "Tumor_Type", "TNM_Stage")
data_clinical[, "PID"] <- trimws(data_clinical[, "PID"])
data_clinical[, "TNM_Stage"] <- gsub("III[A-C]", "III", gsub("Ⅳ", "IV", data_clinical[, "TNM_Stage"]))
data_clinical[, "TNM_Stage"] <- factor(data_clinical[, "TNM_Stage"], levels=c("IIC", "III", "IV", "N.A."))
data_clinical_sub <- data_clinical
data_clinical_sub[, "TNM_Stage"] <- as.character(data_clinical_sub[, "TNM_Stage"])
data_clinical_sub[which(data_clinical_sub[, "TNM_Stage"]=="N.A."), "TNM_Stage"] <- "III/IV"

# data_clinical_sub <- subset(data_clinical, TNM_Stage!="N.A.")
# data_clinical_sub <- subset(data_clinical_sub, TNM_Stage!="IIC")
# data_clinical_sub[, "TNM_Stage"] <- factor(as.character(data_clinical_sub[, "TNM_Stage"]), levels=c("III", "IV"))
data_clinical_sub <- subset(data_clinical_sub, PID!="81122")
write.csv(data_clinical_sub, paste0(res_home, "Table/9.melanoma_acral_mIHC/data_clinical_sub.csv"), quote=FALSE)

library(table1)
res <- table1(~Gender+ Age+Tumor_Type+TNM_Stage|Response, data=data_clinical_sub)
write(res, paste0(res_home, "Figure/9.melanoma_acral_mIHC/test2.html"))


# dat <- data_clinical_sub

col_space <- getColor_palette(pal_custom="yye1", pal_name=NULL, col_num=30)
data_clinical_sub <- data_clinical_sub[order(data_clinical_sub[,"Response"]), ]
# data_clinical_sub[,"PID"] <- factor(data_clinical_sub[,"PID"], levels=data_clinical_sub[,"PID"])
library(ggplot2)
library(ggnewscale)
pHeat <- ggplot(data = data_clinical_sub, aes(x=PID)) +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=6, fill=Response)) +
  scale_fill_manual(limits=c("NR", "R"),values = c("#E69F00","#56B4E9")) +
  new_scale_fill() +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=5, fill=Gender)) +
  scale_fill_manual(limits=c("Male","Female"), values = col_space[c(14,18)], labels=c("Male","Female"), name="Gender") +
  new_scale_fill() +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=4, fill=Age)) +
  scale_fill_distiller(palette = "Blues", direction = 1, name="Age") + 
  new_scale_fill() +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=3, fill=Tumor_Type))+
  scale_fill_manual(limits=c("Acral", "Nodular"),values = col_space[c(2,4)], name="Subtypes") +
  new_scale_fill() +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=2, fill=TNM_Stage))+
  scale_fill_manual(limits=c("IIC", "III", "IV", "III/IV"),values = col_space[c(8:11)], name="Stage") +
  new_scale_fill() +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=1, fill=PFS_month))+
  scale_fill_gradient2(high = "purple", mid="white", low = "#56B1F7", space = "Lab", midpoint = 5.5, guide = "colourbar", aesthetics = "fill") +
  new_scale_fill() +
  geom_tile(data = data_clinical_sub, aes(x=PID, y=0, fill=OS_month))+
  scale_fill_gradient2(high = "red", mid="white", low = "#56B1F7", space = "Lab", midpoint = 16, guide = "colourbar", aesthetics = "fill") +
  theme(panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key = element_blank(),
    title = element_text(size = 8),
    legend.position = "bottom",
    legend.key.size = unit(5,"pt"), legend.text.align = 0, legend.title.align = 0, legend.direction = "vertical") +
  scale_y_discrete(limits=c(0,1,2,3,4,5,6), labels=c("OS.time (Month)","PFS.time (Month)","Stage","Subtype","Age","Gender","Response")) +
  labs(x="", subtitle = "Clinical characteristics of in house cohort")

outfile_tmp <- paste0(res_home, "Figure/9.melanoma_acral_mIHC/characteristics.patients.pdf")
ggsave(plot=pHeat, filename=outfile_tmp, width=7, height=4)


# data_clinical



data_exprs <- read_xls2(paste0(data_home, "免疫治疗 号子+预后 分析-1.1.xlsx"), sheet=2)
data_exprs <- t(data_exprs)
colnames(data_exprs) <- data_exprs[1, ]
data_exprs <- data_exprs[-1, ]
data_exprs <- data.frame(PID=rownames(data_exprs), data_exprs)

sort(setdiff(data_clinical_sub[, "PID"], data_exprs[,"PID"]))
sort(setdiff(data_exprs[, "PID"], data_clinical_sub[,"PID"]))
# > sort(setdiff(data_clinical_sub[, "PID"], data_exprs[,"PID"]))
# [1] "1362455"  "140376-3" "70213-6"  "72610"    "76645-1"  "79063"    "80406-1"
# [8] "81122"
# > sort(setdiff(data_exprs[, "PID"], data_clinical_sub[,"PID"]))
# [1] "136245"   "143076-3" "70231"    "72601"    "74475"    "76443"    "76645"
# [8] "79603"    "80406"
data_exprs[which(data_exprs[, "PID"]=="136245"), "PID"] <- "1362455"
data_exprs[which(data_exprs[, "PID"]=="143076-3"), "PID"] <- "140376-3"
data_exprs[which(data_exprs[, "PID"]=="70231"), "PID"] <- "70213-6"
data_exprs[which(data_exprs[, "PID"]=="72601"), "PID"] <- "72610"
# data_exprs[which(data_exprs[, "PID"]=="74475"), "PID"] <- "76645-1"
data_exprs[which(data_exprs[, "PID"]=="79603"), "PID"] <- "79063"
data_exprs[which(data_exprs[, "PID"]=="80406"), "PID"] <- "80406-1"
data_exprs[which(data_exprs[, "PID"]=="76645"), "PID"] <- "76645-1"





data_clinical_exprs <- merge(data_clinical_sub, data_exprs, by="PID")
data_clinical_exprs[, "CD137L"] <- as.numeric(trimws(data_clinical_exprs[, "CD137L"]))
data_clinical_exprs[, "PDL1"] <- as.numeric(trimws(data_clinical_exprs[, "PDL1"]))

rownames(data_clinical_exprs) <- data_clinical_exprs[,"PID"]
save(data_clinical_exprs, file=paste0(res_home, "RData/9.melanoma_acral_mIHC/data_clinical_exprs.RData"))
table(data_clinical_exprs[, "Response"])
table(data_clinical_exprs[,"TNM_Stage"])
# > table(data_clinical_exprs[, "Response"])
# NR  R
# 15 14
# > table(data_clinical_exprs[,"TNM_Stage"])
# IIC III  IV
#   1  13  14

