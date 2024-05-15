# ========================================================================================================
# Fig. S5 A and B
# ========================================================================================================
# 设置目录
checkDir(paste0(res_home, "RData/mouse_AICAR/"))
checkDir(paste0(res_home, "Figure/mouse_AICAR/"))
checkDir(paste0(res_home, "Table/mouse_AICAR/"))


# 导入数据 
data_diff <- read.csv(paste0(data_home_spec, "AICAR_ShExp_DEseq2_Result.csv"))
# ===================================================================================================================
# 功能富集 -- GSEA
# ===================================================================================================================
# 提取差异基因排秩
rank_gene_tmp <- data_diff[, "log2FoldChange"]
names(rank_gene_tmp) <- data_diff[, "gene_id"]
rank_gene_tmp <- sort(na.omit(rank_gene_tmp), decreasing=TRUE)
# 富集分析
gsea_res_immunesigdb <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c7.immunesigdb.v7.5.1.symbols.gmt")
gsea_res_biocarta <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c2.cp.biocarta.v7.5.1.symbols.gmt")
gsea_res_GO <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c5.go.bp.v7.5.1.symbols.gmt")
gsea_res_KEGG <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c2.cp.kegg.v7.5.1.symbols.gmt")
gsea_res_Reactome <- clusterProfiler_enrich(inputgenes=rank_gene_tmp, min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GSEA", GSEA_gmtfile="/home/pub/Data/Gene_Set/MsigDB/version7.5/c2.cp.reactome.v7.5.1.symbols.gmt")
save(gsea_res_biocarta, file=paste0(res_home, "RData/mouse_AICAR/gsea_res_biocarta.RData"))
save(gsea_res_GO, file=paste0(res_home, "RData/mouse_AICAR/gsea_res_GO.RData"))
save(gsea_res_KEGG, file=paste0(res_home, "RData/mouse_AICAR/gsea_res_KEGG.RData"))
save(gsea_res_Reactome, file=paste0(res_home, "RData/mouse_AICAR/gsea_res_Reactome.RData"))
save(gsea_res_immunesigdb, file=paste0(res_home, "RData/mouse_AICAR/gsea_res_immunesigdb.RData"))

# 可视化
gsea_res_list <- list(biocarta=gsea_res_biocarta, GO=gsea_res_GO, KEGG=gsea_res_KEGG,
  Reactome=gsea_res_Reactome, immunesigdb=gsea_res_immunesigdb)
for(fun_type in names(gsea_res_list)){
  gsea_res_tmp <- gsea_res_list[[fun_type]]
  # 输出富集结果
  if(dim(gsea_res_tmp)[1]>0){
    outfile_tmp <- paste0(res_home, "Figure/mouse_AICAR/gsea_res_", fun_type, ".pdf")
    res <- clusterProfiler_plot(enrich_obj=gsea_res_tmp, outfile=outfile_tmp, 
      showCategory=50, gseaplot=TRUE, geneSetID=1, child_width=7)
  }
}
 
# ===================================================================================================================
# 功能富集 -- 超几何
# ===================================================================================================================
# 提取显著差异的基因
data_diff_sig <- subset(data_diff, padj<0.05)
data_diff_sig <- subset(data_diff_sig, abs(log2FoldChange)>1)
diff_genes <- unique(data_diff_sig$gene_id)
diff_genes <- id_convert(input=diff_genes, from="SYMBOL", to="ENTREZID")
# 功能富集
fun_res_GO <- clusterProfiler_enrich(inputgenes=na.omit(diff_genes$ENTREZID), min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="GO")
fun_res_KEGG <- clusterProfiler_enrich(inputgenes=na.omit(diff_genes$ENTREZID), min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="KEGG")
fun_res_Reactome <- clusterProfiler_enrich(inputgenes=na.omit(diff_genes$ENTREZID), min_size=5, cutoff=0.05,
  adjustMethod="BH", enrich_type="Reactome")
save(fun_res_GO, file=paste0(res_home, "RData/mouse_AICAR/fun_res_GO.RData"))
save(fun_res_KEGG, file=paste0(res_home, "RData/mouse_AICAR/fun_res_KEGG.RData"))
save(fun_res_Reactome, file=paste0(res_home, "RData/mouse_AICAR/fun_res_Reactome.RData"))

fun_res_list <- list(GO=fun_res_GO, KEGG=fun_res_KEGG, Reactome=fun_res_Reactome)
# 可视化
for(fun_type in names(fun_res_list)){
  fun_res_tmp <- fun_res_list[[fun_type]]
  # 输出富集结果
  if(dim(fun_res_tmp)[1]>0){
    outfile_tmp <- paste0(res_home, "Figure/mouse_AICAR/fun_res_", fun_type, ".pdf")
    res <- clusterProfiler_plot(enrich_obj=fun_res_tmp, outfile=outfile_tmp, 
      showCategory=50, gseaplot=FALSE, geneSetID=1, child_width=7)
  }
}
 
fun_res_KEGG_sub <- fun_res_KEGG[grep("cancer", fun_res_KEGG$Description, invert=TRUE), , asis=TRUE]

rm_idx <- c(5:7, 12, 15, 17, 22:23, 26:29, 32, 36, 38, 42:43, 45, 47:50, 52)
fun_res_KEGG_sub <- fun_res_KEGG_sub[-rm_idx, , asis=TRUE]
outfile_tmp <- paste0(res_home, "Figure/mouse_AICAR/fun_res_KEGG_sub.pdf")
ggobj_list <- clusterProfiler_plot(enrich_obj=fun_res_KEGG_sub, outfile=outfile_tmp, 
  showCategory=50, gseaplot=FALSE, geneSetID=1, child_width=7)

# 输出结果
outfile_tmp <- paste0(res_home, "Figure/mouse_AICAR/fun_res_KEGG_sub.part1.pdf")
ggobj_barplot <- barplot(fun_res_KEGG_sub, showCategory=30) #+ theme(axis.text.y=element_text(size=10))
ggsave(plot=ggobj_barplot, filename=outfile_tmp, width=7)

