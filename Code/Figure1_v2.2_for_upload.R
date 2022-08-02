# 2022.03.01
# heatmap without MUT/VLM & VOUS annotation. Suitable for AI. 

# Figure 1: column annotation, color amendment & TMB update
rm(list = ls())

library(dplyr)
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(openxlsx)
library(eoffice)

# pass in working directory
args <- commandArgs(T)
data_path <- paste0(args[1], "Data/Figure1/")
res_path <- paste0(args[1], "Results/Figure1/")
dir.create(res_path)

# #################################### Heatmap data ##################################
# ############### copy from <V4_heatmap_final_snv_cnv_table_20220209.R> ##############
# ###### Codes before embellishment. Thanks to pervious work done by colleagues ######
# ####################################################################################
# get_mut <- function(path, res_path) {
# 	# 2022-02-09: 显示TMB level: low mid high, 而不是TMB数值
# 	# 补充了pathway基因DDR2
	
# 	library(dplyr)
# 	library(maftools)
# 	library(ComplexHeatmap)
# 	library(circlize)
# 	library(RColorBrewer)
	
# 	setwd(path)

# 	# 确定的最终版本SNV
# 	dat <- read.csv("snv_id_new_new_add_Variant_Classification_namechanged.csv", stringsAsFactors = F,header = T)
	
# 	# 删除冗余的列和淋巴结的变异，删除OK检测的结果
# 	dat2 <- dat %>% dplyr::select(-c("cp", "start", "end", "ref", "alt", "sample_id")) %>% filter(Group != "lymph" & testing_type != "OK" & hxid != "NA") 
	
# 	dat3 <- dat2 %>% mutate(new_var = 
# 		ifelse(dat2$Variant_Classification == "Frameshift" & (dat2$Reference_Allele == "-" | (nchar(dat2$Reference_Allele) < nchar(dat2$Tumor_Seq_Allele2))) , "Frame_Shift_Ins",
# 		ifelse(dat2$Variant_Classification == "Frameshift" & (dat2$Tumor_Seq_Allele2 == "-" | (nchar(dat2$Reference_Allele) > nchar(dat2$Tumor_Seq_Allele2))), "Frame_Shift_Del",
# 		ifelse(dat2$Variant_Classification == "In-frame indel" & (dat2$Reference_Allele == "-" | (nchar(dat2$Reference_Allele) < nchar(dat2$Tumor_Seq_Allele2))), "In_Frame_Ins",
# 		ifelse(dat2$Variant_Classification == "In-frame indel" & (dat2$Tumor_Seq_Allele2 == "-" | (nchar(dat2$Reference_Allele) > nchar(dat2$Tumor_Seq_Allele2))), "In_Frame_Del",
# 		ifelse(dat2$Variant_Classification == "Missense", "Missense_Mutation", 
# 		ifelse(dat2$Variant_Classification == "Splice", "Splice_Site",
# 		ifelse(dat2$Variant_Classification == "Stopgain", "Nonsense_Mutation", "5'UTR"))))))))
# 	dat3$Variant_Classification <- dat3$new_var
# 	table(dat3$Variant_Classification)
	
# 	samples = read.csv("idcp.samples.with.HXid.update.csv",stringsAsFactors = F)
# 	samples$hxid_short = gsub("HXIDCP-","",samples$hxid)
# 	samples$group[samples$group == "adenocarcinoma"] = "PCA"
# 	samples$new_id = paste(samples$hxid_short, samples$group, sep="_")
	
# 	pos = match(dat3$tumor, samples$Tumor)
# 	dat3 = cbind(dat3, samples[pos,])
# 	dat3$Tumor_Sample_Barcode <- dat3$new_id
# 	dat3 <- dat3[which(!is.na(dat3$new_id)), ]
	
# 	sample_order = unique(sort(dat3$Tumor_Sample_Barcode))
	
# 	# 确定基因
# 	PathwayGenes <- read.csv("PathwayGene_20220209.csv")
	
# 	#########preprocess using maftools
# 	maf=read.maf(maf = dat3)
	
# 	SNV_gene <- PathwayGenes[PathwayGenes$gene %in% dat3$Hugo_Symbol, "gene"]
	
# 	pdf(file.path("./","44.select_PathwayGene_oncoplot.pdf"))
# 	oncoplot(maf = maf,
# 		#top = 30,
# 		showTumorSampleBarcodes = T,
# 		sampleOrder=sample_order,
# 		genes = SNV_gene,
# 		barcode_mar = 4,
# 		gene_mar = 6,	
# 		writeMatrix = T, # 输出矩阵，Complexheatmap Input
# 		keepGeneOrder = T,# GeneOrderSort = F, # 保持基因顺序不要变
# 		removeNonMutated = F, # 不要删除没有突变的样本
# 		annotationFontSize = 0.8)
# 	dev.off()

# 	################################################################################
# 	# 增加TMB
# 	# 53582S10 06_PCA手动计算 
# 	TMB_MSI <- read.csv("tmb_wes.csv", header = T)
# 	pos = match(TMB_MSI$sampleID, samples$Tumor)
# 	TMB_MSI = cbind(TMB_MSI, samples[pos,])
# 	TMB_MSI <- na.omit(TMB_MSI)
# 	TMB_MSI <- TMB_MSI[order(TMB_MSI$new_id), ] 
# 	TMB_MSI[TMB_MSI$Level=="low", "Level"] <- "Low"
# 	TMB_MSI$Level <- factor(TMB_MSI$Level, levels = c("Low", "Mid", "High"))
# 	write.csv(TMB_MSI, "TMB_MSI.csv", quote = F)
	
# 	################################################################################
# 	# 增加CNV burden
# 	CNV_burden <- read.csv("cnv_burden_20211214.csv")
# 	CNV_burden$cnv_burden <- CNV_burden$Total_length/3137161264
	
# 	pos <- match(CNV_burden$Tumor, samples$Tumor)
# 	CNV_burden <- cbind(CNV_burden, samples[pos,])
# 	CNV_burden <- CNV_burden[order(CNV_burden$new_id), ]
# 	write.csv(CNV_burden, "CNV_burden.csv", quote = F)
	
# 	################################################################################
# 	# Mutation signature 
# 	Mut_Sig <- read.csv("Percent_weights_6_signatures.csv") 
# 	Mut_Sig$X <- sub("HXIDCP-", "", Mut_Sig$X) 
# 	Mut_Sig$X <- sub("-", "_", Mut_Sig$X)
# 	rownames(Mut_Sig) <- Mut_Sig[,1]
# 	Mut_Sig <- Mut_Sig[,-1]
# 	write.csv(Mut_Sig, "Mut_Sig.csv", quote = F)
	
# 	################################################################################
# 	# Tumor purity
# 	TF = read.table("80.samples.list",stringsAsFactors = F,header = T)
# 	pos <- match(TF$tumor, samples$Tumor)
# 	TF <- cbind(TF, samples[pos,])
# 	TF <- TF[order(TF$new_id), ]
# 	TF <- na.omit(TF)
# 	write.csv(TF, "TF.csv", quote = F)

# 	# 2021-12-24
# 	# 增加了提供的临床信息，但是目前信息不完整
# 	# 2022-01-12 修改了PSA的结果 患者基本信息-20220107.xlsx
# 	clin_info <- read.delim("20220114_clinical_info.txt", header = T)
# 	clin_info$Baseline.PSA <- as.numeric(clin_info$Baseline_PSA)
# 	clin_info$Neoadjuvant.ADT <- as.factor(clin_info$NeoadjuvantADT)
# 	clin_info$T.stage <- as.factor(clin_info$Tstage)
# 	clin_info$N.stage <- as.factor(clin_info$Nstage)
# 	pos <- match(clin_info$new_id, samples$new_id)
# 	clin_info <- cbind(clin_info, samples[pos,])
# 	clin_info <- clin_info[order(clin_info$new_id), ]
	
# 	# 增加CNV信息
# 	cnv <- read.csv("cnv_id_new_new_20211224.csv")
# 	cnv_filter <- cnv %>% dplyr::filter(Group != "lymph" & testing_type == "WES" & hxid != "NA") 
# 	cnv_filter$Sub_group[cnv_filter$Sub_group == "adenomatous"] <- "PCA"
# 	cnv_filter$new_id <- paste(unlist(lapply(cnv_filter$hxid, function(x){strsplit(x, "-")[[1]][2]})), cnv_filter$Sub_group, sep = "_")
	
# 	# 构建CNV gain loss突变矩阵
# 	sampleId <- unique(dat3$new_id)
# 	gene <- PathwayGenes$gene
	
# 	mutationMatrix <- matrix(data = "", ncol = length(sampleId), nrow = length(gene))
# 	colnames(mutationMatrix) <- sampleId
# 	row.names(mutationMatrix) <- gene
# 	loc1 <- match(cnv_filter$gene, gene)
# 	loc2 <- match(cnv_filter$new_id, sampleId)
# 	# 赋值gain loss
# 	for(i in 1:nrow(cnv_filter)){
# 		mutationMatrix[loc1[i], loc2[i]] <- cnv_filter[i, "type"]
# 	}
# 	df <- as.data.frame(mutationMatrix)
	
# 	# 在onco_matrix.csv基础上修改增加CNV的信息
# 	onco_table <- read.csv("onco_matrix.csv", header = T)
# 	colnames(onco_table) <- sub("X", "", colnames(onco_table))
# 	onco_table[is.na(onco_table)] = ""
# 	rownames(onco_table) <- onco_table[,1]
# 	onco_table = onco_table[, -1]
	
# 	# SNV非空 & CNV空，则将SNV的结果覆盖掉CNV的位置
# 	# 实现合并SNV和CNV
# 	for (i in rownames(onco_table)) {
# 		for (j in colnames(onco_table)) {
# 			if ((df[i,j] == "") & (onco_table[i,j] != "")) {
# 				df[i,j] <- onco_table[i,j] 
# 			}#else if(!is.na(df[i,j]) & !is.na(onco_table[i,j])){
# 			#	df[i,j] <- paste(df[i,j], onco_table[i,j], sep = ";")
# 			#}
# 		}
# 	}
	
# 	mut <- as.matrix(df)

# 	col = c("In_Frame_Del" = "#80b1d3", "In_Frame_Ins" = "#fb8072", "Missense_Mutation" = "#8dd3c7", "Nonsense_Mutation" = "#ffffb3", "Splice_Site" = "#bebada", "Frame_Shift_Ins" = "#fdb462", "Frame_Shift_Del" = "#b3de69", "Multi_Hit" = "#bc80bd", "gain" = "#fccde5", "loss" = "#d9d9d9")
	
# 	alter_fun = list(
# 		background = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"), 
# 			gp = gpar(fill = "#DDDDDD", col = "white"))},
# 		In_Frame_Del = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["In_Frame_Del"], col = "white"))},
# 		In_Frame_Ins = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["In_Frame_Ins"], col = "white"))},
# 		Missense_Mutation = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Missense_Mutation"], col = "white"))},
# 		Nonsense_Mutation = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Nonsense_Mutation"], col = "white"))},
# 		Splice_Site = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Splice_Site"], col = "white"))},
# 		Frame_Shift_Ins = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Frame_Shift_Ins"], col = "white"))},
# 		Frame_Shift_Del = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Frame_Shift_Del"], col = "white"))},
# 		Frame_Shift_Del = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Frame_Shift_Del"], col = "white"))},
# 		Multi_Hit = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["Multi_Hit"], col = "white"))},
# 		gain = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["gain"], col = "white"))},
# 		loss = function(x, y, w, h) {
# 			grid.rect(x, y, w-unit(0.1, "mm"), h-unit(0.1, "mm"),
# 			gp = gpar(fill = col["loss"], col = "white"))}
# 	)
	
# 	heatmap_legend_param = list(
# 		title = "SNVs and CNVs", at = c("In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "Frame_Shift_Del", "Multi_Hit", "gain", "loss"), 
# 		labels = c("In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "Frame_Shift_Del", "Multi_Hit", "CNV gain", "CNV loss"))
	
# 	set.seed(1)
	
# 	Mut_Sig_col <- c("Signature.1" = "#248DCA", "Signature.3" = "#FC8412", "Signature.5" = "#BFD595", "Signature.6" = "#BC9AC5", "Signature.15" = "#BC4E97", "Others" = "#878787")
# 	ha_bottom <- HeatmapAnnotation(
# 		'Mutation Signature' = anno_barplot(Mut_Sig, bar_width = 1, gp = gpar(fill = Mut_Sig_col), height = unit(2, "cm")),
# 		Group = TF$group,
# 		col = list('Mutation Signature' = c("Signature.1" = "#B91934", "Signature.3" = "#F16A45", "Signature.5" = "#FAD16C", "Signature.6" = "#BC9AC5", "Signature.15" = "#BC4E97", "Others" = "#878787"), 'Group' = c("IDCP" = "#ADCF34", "PCA" = "#78BEC2")),
# 		show_annotation_name = TRUE,
# 		annotation_name_gp = gpar(fontsize = 10, color = "black")
# 	)
	
# 	lgd_list = list(
# 	Legend(labels = c("Signature.1", "Signature.3", "Signature.5", "Signature.6", "Signature.15", "Others"), title = "Mutation Signature", type = "points", pch = 16, legend_gp = gpar(col = c("#248DCA", "#FC8412", "#BFD595", "#4199B8", "#2D72A1", "#878787")), background = "white"))
	
# 	col_psa = colorRamp2(c(0, 100), c("#B2D9D8", "#0B5C51"))
# 	col_CNV = colorRamp2(c(0, 1), c("white", "red"))
# 	col_Pathological = colorRamp2(c(0.1, 0.6), c("#FFC0CB", "#FF1493"))
# 	col_age = colorRamp2(c(49, 80), c("white", "#373D61"))
	
# 	ha_top <- HeatmapAnnotation(Age = clin_info$Age,
# 		'Baseline PSA' = clin_info$Baseline.PSA,
# 		'ISUP grading for PSA' = clin_info$ISUP_grading_for_PAC,
# 		'T stage' = clin_info$T.stage,
# 		'N stage' = clin_info$N.stage,
# 		'Neoadjuvant ADT' = clin_info$Neoadjuvant.ADT,
# 		'CNV burden' = CNV_burden$cnv_burden,
# 		'TMB level' = TMB_MSI$Level,
# 		'Tumor fraction' = TF$purity,
# 		'Pathological IDC-P proportion' = clin_info$Pathological_IDCP_proportion,
# 		cbar = anno_oncoprint_barplot(),
# 		#指定颜色
# 		col = list('Gleason score' = c("3+4=7" = "#E6E2D8", "4+3=7" = "#D7EBC5", "4+4=8" = "#BDD99A", "4+5=9" = "#73A86F", "5+4=9" = "#4B7C99"),
# 		'Age' = col_age,
# 		'ISUP grading for PSA' = c("Group 2" = "#4199B8", "Group 3" = "#2D72A1", "Group 4" = "#255E81", "Group 5" = "#195963"), 
# 		'T stage' = c("T2a" = "#DC8C75","T2c" = "#BF534B", "T3a" = "#A4292B", "T3b"="#8F2748"),
# 		'N stage' = c("0" = "#BC9AC5", "1" = "#BC4E97", "2" = "#2A2A40", "x"="#F53183"),
# 		'Neoadjuvant ADT' = c("0" = "#65BBA4", "1" = "#489D43"),
# 		'Baseline PSA' = col_psa,
# 		'CNV burden' = col_CNV,
# 		'TMB level' = c("Low" = "#E2ED63", "Mid" = "#BCD451", "High" = "#78962D"),
# 		'Pathological IDC-P proportion'= col_Pathological),
# 		show_annotation_name = TRUE,
# 		annotation_name_gp = gpar(fontsize = 10, color = "black")
# 	)
	
# 	# 保存t(mut)
# 	a <- as.data.frame(mut)
# 	write.table(as.data.frame(mut), "landscape_mutation.txt", quote = F, row.names = T, sep = "\t")
	
# 	oncoplot_combined <- oncoPrint(mut,
# 		alter_fun = alter_fun, 
# 		remove_empty_columns = FALSE, remove_empty_rows = FALSE,
# 		alter_fun_is_vectorized = FALSE,
# 		col = col,
# 		pct_side = "right", 
# 		row_names_side = "left",
# 		show_column_names = T, 
# 		row_order = PathwayGenes$gene, 
# 		column_order = colnames(mut),
# 		# top_annotation = ha_top,
# 		# bottom_annotation = ha_bottom,
# 		heatmap_legend_param = heatmap_legend_param,
# 		height=unit(114.114,"mm"))
# 	# draw(oncoplot, annotation_legend_list = lgd_list, annotation_legend_side = "right")
# 	# draw(oncoplot, annotation_legend_side = "right")
# 	# PDF 9 14
	
	
# 	# 2022-01-25 
# 	# 为了获取指定基因的变异数barplot
# 	# 单独绘制IDCP的landscape
# 	mut_IDCP <- read.delim("landscape_mutation_IDCP_20220209.txt", header = T, sep = "\t")
# 	mut_IDCP[is.na(mut_IDCP)] <- ""
# 	rownames(mut_IDCP) <- mut_IDCP[,1]
# 	mut_IDCP = mut_IDCP[, -1]
# 	colnames(mut_IDCP) <- sub("X","", colnames(mut_IDCP))
	
# 	clin_info_IDCP <- clin_info[which(clin_info$group=="IDCP"),]
# 	CNV_burden_IDCP <- CNV_burden[which(CNV_burden$group=="IDCP"),]
# 	TMB_MSI_IDCP <- TMB_MSI[which(TMB_MSI$group=="IDCP"),]
# 	TF_IDCP <- TF[which(TF$group=="IDCP"),]
	
# 	ha_top_IDCP <- HeatmapAnnotation(Age = clin_info_IDCP$Age,
# 		'Baseline PSA' = clin_info_IDCP$Baseline.PSA,
# 		'ISUP grading for PSA' = clin_info_IDCP$ISUP_grading_for_PAC,
# 		'T stage' = clin_info_IDCP$T.stage,
# 		'N stage' = clin_info_IDCP$N.stage,
# 		'Neoadjuvant ADT' = clin_info_IDCP$Neoadjuvant.ADT,
# 		'CNV burden' = CNV_burden_IDCP$cnv_burden,
# 		'TMB level' = TMB_MSI_IDCP$Level,
# 		'Tumor fraction' = TF_IDCP$purity,
# 		'Pathological IDC-P proportion' = clin_info_IDCP$Pathological_IDCP_proportion,
# 		cbar = anno_oncoprint_barplot(),
# 		#指定颜色
		
		
# 		col = list('Gleason score' = c("3+4=7" = "#E6E2D8", "4+3=7" = "#D7EBC5", "4+4=8" = "#BDD99A", "4+5=9" = "#73A86F", "5+4=9" = "#4B7C99"),
# 		'Age' = col_age,
# 		'ISUP grading for PSA' = c("Group 2" = "#4199B8", "Group 3" = "#2D72A1", "Group 4" = "#255E81", "Group 5" = "#195963"), 
# 		'T stage' = c("T2a" = "#DC8C75","T2c" = "#BF534B", "T3a" = "#A4292B", "T3b"="#8F2748"),
# 		'N stage' = c("0" = "#BC9AC5", "1" = "#BC4E97", "2" = "#2A2A40", "x"="#F53183"),
# 		'Neoadjuvant ADT' = c("0" = "#65BBA4", "1" = "#489D43"),
# 		'Baseline PSA' = col_psa,
# 		'CNV burden' = col_CNV,
# 		'TMB level' =c("Low" = "#E2ED63", "Mid" = "#BCD451", "High" = "#78962D"),
# 		'Pathological IDC-P proportion'= col_Pathological),
# 		show_annotation_name = TRUE,
# 		annotation_name_gp = gpar(fontsize = 10, color = "black")
# 	)
	
	
# 	Mut_Sig_IDCP <- Mut_Sig[grep("IDCP", rownames(Mut_Sig)), ]
# 	ha_bottom_IDCP <- HeatmapAnnotation('Mutation Signature' = anno_barplot(Mut_Sig_IDCP, bar_width = 1, gp = gpar(fill = Mut_Sig_col), height = unit(2, "cm")),
# 		Group = TF_IDCP$group,
# 		col = list('Mutation Signature' = c("Signature.1" = "#B91934", "Signature.3" = "#F16A45", "Signature.5" = "#FAD16C", "Signature.6" = "#DDF680", "Signature.15" = "#7AC470", "Others" = "#878787"),
# 		 'Group' = c("IDCP" = "#ADCF34", "PCA" = "#78BEC2")
# 		),
# 		show_annotation_name = TRUE,
# 		annotation_name_gp = gpar(fontsize = 10, color = "black")
# 	)
	
# 	oncoplot_IDCP <- oncoPrint(mut_IDCP,
# 		alter_fun = alter_fun, 
# 		remove_empty_columns = FALSE, remove_empty_rows = FALSE,
# 		alter_fun_is_vectorized = FALSE,
# 		col = col,
# 		pct_side = "right", 
# 		row_names_side = "left",
# 		show_column_names = T, 
# 		row_order = PathwayGenes$gene, 
# 		column_order = colnames(mut_IDCP),
# 		# top_annotation = ha_top_IDCP,
# 		# bottom_annotation = ha_bottom_IDCP,
# 		heatmap_legend_param = heatmap_legend_param,
# 		height=unit(114.114,"mm"))
# 	# draw(oncoplot, annotation_legend_list = lgd_list, annotation_legend_side = "right")
# 	# 9 12	


# 	# 2.单独绘制PCA的landscape
# 	mut_PCA <- read.delim("landscape_mutation_PCA_20220209.txt", header = T, sep = "\t")
# 	mut_PCA[is.na(mut_PCA)] <- ""
# 	rownames(mut_PCA) <- mut_PCA[,1]
# 	mut_PCA = mut_PCA[, -1]
# 	colnames(mut_PCA) <- sub("X","", colnames(mut_PCA))
	
# 	clin_info_PCA <- clin_info[which(clin_info$group=="PCA"),]
# 	CNV_burden_PCA <- CNV_burden[which(CNV_burden$group=="PCA"),]
# 	TMB_MSI_PCA <- TMB_MSI[which(TMB_MSI$group=="PCA"),]
# 	TF_PCA <- TF[which(TF$group=="PCA"),]
	
# 	ha_top_PCA <- HeatmapAnnotation(Age = clin_info_PCA$Age,
# 		'Baseline PSA' = clin_info_PCA$Baseline.PSA,
# 		'ISUP grading for PSA' = clin_info_PCA$ISUP_grading_for_PAC,
# 		'T stage' = clin_info_PCA$T.stage,
# 		'N stage' = clin_info_PCA$N.stage,
# 		'Neoadjuvant ADT' = clin_info_PCA$Neoadjuvant.ADT,
# 		'CNV burden' = CNV_burden_PCA$cnv_burden,
# 		'TMB level' = TMB_MSI_PCA$Level,
# 		'Tumor fraction' = TF_PCA$purity,
# 		'Pathological IDC-P proportion' = clin_info_PCA$Pathological_IDCP_proportion,
# 		cbar = anno_oncoprint_barplot(),
# 		#指定颜色
	
# 		 col = list('Gleason score' = c("3+4=7" = "#E6E2D8", "4+3=7" = "#D7EBC5", "4+4=8" = "#BDD99A", "4+5=9" = "#73A86F", "5+4=9" = "#4B7C99"),
# 		'Age' = col_age,
# 		'ISUP grading for PSA' = c("Group 2" = "#4199B8", "Group 3" = "#2D72A1", "Group 4" = "#255E81", "Group 5" = "#195963"), 
# 		'T stage' = c("T2a" = "#DC8C75","T2c" = "#BF534B", "T3a" = "#A4292B", "T3b"="#8F2748"),
# 		'N stage' = c("0" = "#BC9AC5", "1" = "#BC4E97", "2" = "#2A2A40", "x"="#F53183"),
# 		'Neoadjuvant ADT' = c("0" = "#65BBA4", "1" = "#489D43"),
# 		'Baseline PSA' = col_psa,
# 		'CNV burden' = col_CNV,
# 		'TMB level' = c("Low" = "#E2ED63", "Mid" = "#BCD451", "High" = "#78962D"),
# 		'Pathological IDC-P proportion'= col_Pathological),
# 		show_annotation_name = TRUE,
# 		annotation_name_gp = gpar(fontsize = 10, color = "black")
# 	)
	
# 	Mut_Sig_PCA <- Mut_Sig[grep("PCA", rownames(Mut_Sig)), ]
# 	ha_bottom_PCA <- HeatmapAnnotation('Mutation Signature' = anno_barplot(Mut_Sig_PCA, bar_width = 1, gp = gpar(fill = Mut_Sig_col), height = unit(2, "cm")),
# 		Group = TF_PCA$group,
# 		col = list('Mutation Signature' = c("Signature.1" = "#B91934", "Signature.3" = "#F16A45", "Signature.5" = "#FAD16C", "Signature.6" = "#DDF680", "Signature.15" = "#7AC470", "Others" = "#878787"), 'Group' = c("IDCP" = "#ADCF34", "PCA" = "#78BEC2")),
# 		show_annotation_name = TRUE,
# 		annotation_name_gp = gpar(fontsize = 10, color = "black")
# 	)
	
# 	oncoplot_PCA <- oncoPrint(mut_PCA,
# 		alter_fun = alter_fun, 
# 		remove_empty_columns = FALSE, remove_empty_rows = FALSE,
# 		alter_fun_is_vectorized = FALSE,
# 		col = col,
# 		pct_side = "right", 
# 		row_names_side = "left",
# 		show_column_names = T, 
# 		row_order = PathwayGenes$gene, 
# 		column_order = colnames(mut_PCA),
# 		# top_annotation = ha_top_PCA,
# 		# bottom_annotation = ha_bottom_PCA,
# 		heatmap_legend_param = heatmap_legend_param,
# 		height=unit(114.114,"mm"))
# 	# draw(oncoplot_PCA, annotation_legend_list = lgd_list, annotation_legend_side = "right")
	

# 	############## retuen usefull results	
# 	pdf(paste0(res_path, "v2.1_0_oncoplot_ratio.pdf"))
# 	print(oncoplot_combined)
# 	print(oncoplot_IDCP)
# 	print(oncoplot_PCA)
# 	dev.off()


# 	return(mut)
# }

# mut <- get_mut(paste0(args[1], "Data/Figure1/pre_embellishment/"), res_path)
# write.xlsx(as.data.frame(mut), paste0(res_path, "mut.xlsx"), rowNames=TRUE, overwrite=TRUE)

mut <- as.matrix(read.xlsx(paste0(res_path, "mut_manual.xlsx"), rowNames=TRUE))
mut_patho <- as.matrix(read.xlsx(paste0(res_path, "mut_patho.xlsx"), rowNames=TRUE))


################################## top annotation ##################################
# load data
clin_data_ori <- read.xlsx(paste0(data_path, "Oncoprint_clinical.xlsx"), rowNames=TRUE, sep.names= " ")
tmb <- read.xlsx(paste0(data_path, "TMB.xlsx"))
rownames(tmb) <- paste0(tmb$CASE, "_", tmb$group)

# integrate data to one clinical data matrix
clin_data_pre <- clin_data_ori
for(row in 1:nrow(clin_data_ori)){
	clin_data_pre[(2*row-1):(2*row),] <- clin_data_ori[row,]
	rownames(clin_data_pre)[(2*row-1):(2*row)] <- paste0(rownames(clin_data_ori[row,]), "_", c("IDC-P", "PAC"))
}
clin_data <- merge(clin_data_pre, tmb, by=0, all.x=TRUE)
clin_data$`Pathological pattern for IDC-P` <- factor(clin_data$`Pathological pattern for IDC-P`, levels=c("Solid", "Dense cribriform", "Loose cribriform"))
clin_data$`NCCN risk group` <- factor(clin_data$`NCCN risk group`, levels=c("Intermediate risk", "High risk", "Very high risk"))

# heatmap top annotation
column_color <- list(`Age (Year)`=c(">70"="#d53e4f","≤70"="#eeb1b8"),
	`Baseline PSA (ng/mL)`=c(">20"="#f46d43","≤20"="#fac4b3"),
	`Gleason pattern for PAC`=c("3"="#ffffbf","4"="#fee08b","5"="#fdae61"),
	`Pathological pattern for IDC-P`=c("Solid"="#e6f598","Dense cribriform"="#abdda4","Loose cribriform"="#66c2a5"),
	`NCCN risk group`=c("Intermediate risk"="#adcfe4","High risk"="#6fabd0","Very high risk"="#3288bd"),
	`Neoadjuvant ADT`=c("No"="#beb8d9","Yes"="#5e4fa2"))
column_ha_top = HeatmapAnnotation(
	`Age (Year)`=clin_data$`Age (Year)`,
	`Baseline PSA (ng/mL)`=clin_data$`Baseline PSA (ng/mL)`,
	`Gleason pattern for PAC`=clin_data$`Gleason pattern for PAC`,
	`Pathological pattern for IDC-P`=clin_data$`Pathological pattern for IDC-P`,
	`NCCN risk group`=clin_data$`NCCN risk group`,
	`Neoadjuvant ADT`=clin_data$`Neoadjuvant ADT`,
	TMB=anno_barplot(as.numeric(clin_data$TMB), gp=gpar(fill="#8e83bd", border=NA, lty="blank"), ylim=range(0,5.2), border=FALSE, axis=TRUE, axis_param=list(side="left", at=c(0,1,2,3,4,5)), height=unit(15, "mm")),
	col=column_color,
	gp=gpar(col="white"),
	show_legend=FALSE,
	show_annotation_name=TRUE,
	gap=unit(0.5,"mm"),
	# height=unit(35, "mm"),
	# width=unit(155.155, "mm"),
	annotation_name_gp=gpar(fontsize=10))
column_ha_top_legend = HeatmapAnnotation(
  `Age (Year)`=clin_data$`Age (Year)`,
  `Baseline PSA (ng/mL)`=clin_data$`Baseline PSA (ng/mL)`,
  `Gleason pattern for PAC`=clin_data$`Gleason pattern for PAC`,
  `Pathological pattern for IDC-P`=clin_data$`Pathological pattern for IDC-P`,
  `NCCN risk group`=clin_data$`NCCN risk group`,
  `Neoadjuvant ADT`=clin_data$`Neoadjuvant ADT`,
  TMB=anno_barplot(as.numeric(clin_data$TMB), gp=gpar(fill="#8e83bd", border=NA, lty="blank"), ylim=range(0,5.2), border=FALSE, axis=TRUE, axis_param=list(side="left", at=c(0,1,2,3,4,5)), height=unit(15, "mm")),
  col=column_color,
  gp=gpar(col="white"),
  show_legend=TRUE,
  gap=unit(0.5,"mm"),
  # height=unit(35, "mm"),
  # width=unit(155.155, "mm"),
  annotation_name_gp=gpar(fontsize=10))


################################## bottom annotation ##################################
# load data
mut_sig <- read.xlsx(paste0(data_path, "Mut_Sig.xlsx"))
rownames(mut_sig) <- paste0(mut_sig$CASE, "_", mut_sig$Pathology)

# heatmap bottom annotation
column_color <- list(`IDC-P/PAC`=c("IDC-P"="#cf0101","PAC"="#0072bb"))
column_ha_bottom = HeatmapAnnotation(
	`Mutation Signature`=anno_barplot(mut_sig[,-c(1,2)], gp=gpar(fill=c("#3288bd", "#99d594", "#e6f598", "#fee08b", "#fc8d59", "#f2f2f2"), border=TRUE), bar_width=1, border=FALSE, axis=TRUE, height=unit(18, "mm")),
	`IDC-P/PAC`=mut_sig$Pathology,
	col=column_color,
	gp=gpar(col="white"),
	show_legend=TRUE,
	show_annotation_name=TRUE,
	gap=unit(1,"mm"),
	# height=unit(23, "mm"),
	# width=unit(155.155, "mm"),
	annotation_name_gp=gpar(fontsize=10))

################################## draw annotation ##################################
heatmap_blank <- data.frame(matrix(runif(31*44, 0, 1), nrow=31, ncol=44))
names(heatmap_blank) <- rownames(clin_data)
rownames(heatmap_blank) <- c("AR","FOXA1","NCOR1","NCOR2","APC","FAT1","LRP1B","SOX9","BRAF","HRAS","KRAS","PIK3CA","PIK3C2B","TSC2","CDK12","ATM","BAP1","BRCA2","PMS2","FANCC","MSH6","CHD1","BCOR","KMT2C","ZFHX3","PTEN","RB1","DDR2","SPOP","STAG2","TP53")
legends <- Heatmap(as.matrix(heatmap_blank),
  rect_gp=gpar(type="none"), 
  # heatmap_height=unit(8,"cm"), 
  top_annotation=column_ha_top_legend, 
  bottom_annotation=column_ha_bottom,
  show_column_dend=FALSE, 
  show_column_names=FALSE, 
  show_row_dend=FALSE,
  show_row_names=FALSE,
  row_order=rownames(heatmap_blank),
  column_order=names(heatmap_blank),
  column_title=NULL)

mut_sig_legend <- Legend(title="Mutation Signature", at=c(names(mut_sig[,-c(1,2)]), "MUT/VLM", "VOUS"), legend_gp=gpar(fill=c("#3288bd", "#99d594", "#e6f598", "#fee08b", "#fc8d59", "#f2f2f2", "white", "white")))

pdf(paste0(res_path, "v2.2_0_anno_legend.pdf"), paper="USr", width=unit(279,"mm"), height=unit(216,"mm"))
print(legends)
draw(mut_sig_legend, just=c("right", "bottom"))
dev.off()


################################### Heatmap body ##################################
# locate MUT/VLM & VOUS cells
mut_patho_rev <- apply(mut_patho, 2, rev)
cell_mut <- data.frame(which(mut_patho_rev=="MUT/VLM", arr.ind=TRUE))
cell_vous <- data.frame(which(mut_patho_rev=="VOUS", arr.ind=TRUE))

# heatmap body
heatmap <- Heatmap(mut,
	show_heatmap_legend = TRUE,
	heatmap_legend_param = list(title="SNVs and CNVs"),
	height=unit(114.114, "mm"),
	width=unit(195.155, "mm"),
	na_col="#f2f2f2",
	rect_gp=gpar(col="white", lwd=1.5),
	col=structure(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd"), names=c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", "In_Frame_Ins", "In_Frame_Del", "Frame_Shift_Ins", "Frame_Shift_Del", "CNV gain", "CNV loss", "Multi_Hit")),
	# cell_fun = function(j, i, x, y, width, height, fill){
		# grid.rect(x=unit((cell_mut$col-0.5)/ncol(mut), "npc"), y=unit((cell_mut$row-0.5)/nrow(mut), "npc"), width=width, height=height, gp=gpar(col="black", fill=NA, lwd=1.5))
	# },
	row_names_side="left",
	row_names_gp = gpar(fontsize=8),
	top_annotation=column_ha_top,
	bottom_annotation=column_ha_bottom,
	column_split=as.character(apply(data.frame(paste0("CASE", sprintf("%02d", seq(01,22,1))), paste0("CASE", sprintf("%02d", seq(01,22,1)))), 1, unlist)),
	column_gap=unit(2,"mm"),
	column_title="",
	column_label=as.character(apply(data.frame(rep("",22), paste0("CASE", sprintf("%02d", seq(01,22,1)))), 1, unlist)),
	column_names_rot=55)

pdf(paste0(res_path, "v2.2_0_Figure1_pre_20220302.pdf"), paper="USr", width=unit(279,"mm"), height=unit(232,"mm"))
draw(heatmap)
decorate_annotation("TMB", {
    grid.lines(unit(c(0, 195.155), "mm"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9"))
})
decorate_annotation("TMB", {
    grid.lines(unit(c(0, 195.155), "mm"), unit(c(1, 1), "native"), gp = gpar(col = "#d9d9d9"))
})
decorate_annotation("TMB", {
    grid.lines(unit(c(0, 195.155), "mm"), unit(c(2, 2), "native"), gp = gpar(col = "#d9d9d9"))
})
decorate_annotation("TMB", {
    grid.lines(unit(c(0, 195.155), "mm"), unit(c(3, 3), "native"), gp = gpar(col = "#d9d9d9"))
})
decorate_annotation("TMB", {
    grid.lines(unit(c(0, 195.155), "mm"), unit(c(4, 4), "native"), gp = gpar(col = "#d9d9d9"))
})
decorate_annotation("TMB", {
    grid.lines(unit(c(0, 195.155), "mm"), unit(c(5, 5), "native"), gp = gpar(col = "#d9d9d9"))
})
dev.off()
