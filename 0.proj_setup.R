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
# data_home_independent <- "~/Data/TCGA_secondary/"
# data_home_singlecell <- "~/Data/TCGA_secondary/"

# 设定感兴趣的基因，注意首先检查基因名是否能够在TCGA数据中找到
interest_genes <- c("TNFSF9", "PRKAA1", "HLTF")
key_gene <- "TNFSF9"
interest_cancer_type <- "SKCM"

InAcMarker_TCGAoverlap <- read.csv("~/Data/Gene_Set/InAcMarker_extend.csv")
InAcMarker_TCGAoverlap[, "Role.with.Immunity"] <- gsub("Active", "Activate", InAcMarker_TCGAoverlap[, "Role.with.Immunity"])
InAcMarker_TCGAoverlap[which(InAcMarker_TCGAoverlap[,"Symbol"]=="CD274"), "Role.with.Immunity"] <- "Inhibit"