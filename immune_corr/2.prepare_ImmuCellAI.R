# ========================================================================================================
# Fig. 2 F, prepare data
# ========================================================================================================
checkDir(paste0(res_home, "Figure/phase4.Immune_Infiltration_independent/"))
checkDir(paste0(res_home, "Table/phase4.Immune_Infiltration_independent/")) 
checkDir(paste0(res_home, "RData/phase4.Immune_Infiltration_independent/"))

# # 导入数据 
# load(file=paste0(res_home, "RData/0.prepare_data/data_exprs_list_independent.RData"))
# # https://github.com/lydiaMyr/ImmuCellAI
# # 导入数据 
# load(file="/lustre/pub/Data/Melanoma_bulk_public/Melanoma_GEO/exprs_score/data_immune_list.RData")

# 导入数据 
# load(file=paste0(res_home, "RData/3.validate_independent/data_clinical_list_independent.RData"))
data_exprs_list_independent <- base::get(load(file=paste0(data_home_independent, "melanoma_exprs_list.Rdata")))



library(ImmuCellAI)
res_immune_list <- list()
for(dataset in names(data_exprs_list_independent)){
	data_tmp <- data_exprs_list_independent[[dataset]]
	if(dim(data_tmp)[1]>25000){
		res <- ImmuCellAI(sample=data_tmp, data_type="rnaseq", group_tag=0, response_tag=0, customer=0)
	}else{
		res <- ImmuCellAI(sample=data_tmp, data_type="microarray", group_tag=0, response_tag=0, customer=0)
	}
	
	res_immune_list[[dataset]] <- res[["Sample_abundance"]]
}
save(res_immune_list, file=paste0(res_home, "RData/phase4.Immune_Infiltration_independent/res_immune_list.RData"))

# dataset="GSE22154"


# data_immune_list
#  dataset=="GSE22154"
