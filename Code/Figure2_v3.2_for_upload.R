# Figure 1: column annotation, color amendment & TMB update
rm(list = ls())

library(dplyr)
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(openxlsx)
library(eoffice)
library(plyr)
library(ggplot2)
library(ggbreak) 
library(patchwork)
library(aplot)

# pass in working directory
args <- commandArgs(T)
data_path <- paste0(args[1], "Data/Figure2/")
res_path <- paste0(args[1], "Results/Figure2/")

dir.create(res_path)

############################### data ###############################
# load data
sample_pre <- read.xlsx(paste0(data_path, "sample.xlsx"), rowNames=TRUE, sep.names= " ")
patient <- read.xlsx(paste0(data_path, "patient.xlsx"), rowNames=TRUE, sep.names= " ")

# 2ed & 3rd layer barplot data: mutation fractions per patient
snv_frac <- patient[,1:3]/100
scna_frac <- patient[,4:6]/100


############################### plot ##############################
########################## oncoplot body ##########################
# column annotation
feature_data <- patient[,grep("Evolution", names(patient)):ncol(patient)]
feature_data$`Evolution patterns` <- factor(feature_data$`Evolution patterns`, levels=c("Indepentant origin", "Early divergent", "Late divergent"))
feature_data$`FOXA1 mutation` <- factor(feature_data$`FOXA1 mutation`, levels=c("PAC", "IDC-P", "IDC-P|PAC", "None"))
feature_data$`Spatial distance` <- factor(feature_data$`Spatial distance`, levels=c("Distant", "Closed"))
feature_data$`IDC-P pattern` <- factor(feature_data$`IDC-P pattern`, levels=c("Solid", "Dense cribriform", "Loose cribriform"))

# wanted oncoplot order
oncoplot_order <- rownames(feature_data[order(match(feature_data$`Evolution patterns`, c("Indepentant origin", "Early divergent", "Late divergent"))),])

oncoplot_data <- t(data.frame(feature_data[,-grep("Spatial distance|PAC pattern|IDC-P pattern", names(feature_data))]))
oncoplot_color <- c(
	# Evolution patterns
	"Indepentant origin"="#dbf06a", "Early divergent"="#abdda4", "Late divergent"="#66c2a5", 
	# FOXA1 mutation
	"IDC-P"="#3288bd", "PAC"="#3288bd", "None"="#e5e5e5", 
	# CDK12 mutation
	"IDC-P and PAC"="#3288bd", 
	# SCNA & DNA methylation & mRNA
	"Cluster seperate"="#7e72b4", "Cluster together"="#72b47e", "Unavailable"="#e5e5e5", 
	# Cytological features
	"Resemblance"="#fc9272", "Difference"="#ef3b2c", 
	# PFS ≤ 24 mo
	"No"="#e5e5e5", "Yes"="black")
oncoplot <- oncoPrint(as.matrix(oncoplot_data), 
	alter_fun=list(
		background = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = "#e5e5e5", col = "white"))
			grid.polygon(
				unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
				unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
				gp = gpar(fill = "#e5e5e5", col = "white"))
			},
		# Evolution patterns
		`Indepentant origin` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Indepentant origin"], col = "white"))
		},
		`Early divergent` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Early divergent"], col = "white"))
		},
		`Late divergent` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Late divergent"], col = "white"))
		},
		# FOXA1 mutation & CKD12 mutation
		None = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["None"], col = "white"))
		},
		# FOXA1 mutation
		`IDC-P` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["IDC-P"], col = "white"))
		},
		PAC = function(x, y, w, h) {
			grid.polygon(
				unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
				unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
				gp = gpar(fill = oncoplot_color["PAC"], col = "white"))
		},
		# CDK12 mutation
		`IDC-P and PAC` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["IDC-P and PAC"], col="white"))
		},
		# SCNA & DNA methylation & mRNA
		`Cluster seperate` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Cluster seperate"], col = "white"))
		},
		`Cluster together` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Cluster together"], col = "white"))
		},
		`Unavailable` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Unavailable"], col = "white"))
		},
		# Cytological features
		`Resemblance` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Resemblance"], col = "white"))
		},
		`Difference` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Difference"], col = "white"))
		},
		# PFS ≤ 24 mo
		`No` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["No"], col = "white"))
		},
		`Yes` = function(x, y, w, h) {
			grid.polygon(
				unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w, x + 0.5*w), 
				unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h, y - 0.5*h),
				gp = gpar(fill = oncoplot_color["Yes"], col = "white"))
		}
	),
	col=oncoplot_color,
	height=unit(42, "mm"),
	width=unit(112,"mm"),
	row_order=rownames(oncoplot_data),
	column_order=colnames(oncoplot_data),
	row_split=c("a", rep("b", 5), rep("c", 2)), 
	column_split=factor(oncoplot_data["Evolution.patterns",], levels=c("Indepentant origin", "Early divergent", "Late divergent")),
	show_column_names=TRUE)

pdf(paste0(res_path, "0_oncoplot.pdf"), paper="USr", width=unit(279,"mm"), height=unit(216,"mm"))
print(oncoplot)
dev.off()

########################## top annotation: single barplot ##########################
# 1st layer barplot data
sample <- sample_pre
sample$patient <- gsub("_.*", "", rownames(sample_pre))
sample$feature <- gsub(".*_", "", rownames(sample_pre))
# re-order sample sheet
sample <- sample[order(match(sample$patient, oncoplot_order)),]
sample <- sample[,1:2]

# plot
column_ha = HeatmapAnnotation(
	# `Somatic exonic mutations, count`=function(index){
	# 	n_column = length(index)
	# 	pushViewport(viewport(xscale = c(0, n_column), yscale = c(0, max(sample))))
	# 	n_bar = ncol(sample)
	# 	for(i in seq_along(index)) {  # for each column in the oncoplot
	# 		for(j in 1:n_bar) {
	# 			grid.rect(x = n_column-i + (j-0.5)/n_bar, y = 0, width = unit(2, "mm"), height = unit(sample[i, j]/max(sample)*20, "mm"), just = c("bottom"), default.unit = "native", gp = gpar(fill = j, col = "white"))
	# 		}
	# 	}
	# 	grid.xaxis()
	# 	popViewport()
	# },
	`SNV, fraction`=anno_barplot(snv_frac, gp=gpar(fill=c("#0068b3", "#f7941d", "#f15923"), col="white"), border=FALSE, height=unit(20,"mm"), bar_width=unit(5,"mm")),
	`SCNA, fraction`=anno_barplot(scna_frac, gp=gpar(fill=c("#0068b3", "#f7941d", "#f15923"), col="white"), border=FALSE, height=unit(20,"mm"), bar_width=unit(5,"mm")),
	show_legend=TRUE,
	annotation_name_gp=gpar(fontsize=10),
	width=unit(110,"mm"),
	gap=unit(2,"mm"))

# Legends
legend1 <- Legend(title="Somatic exonic mutation, count", at=c("Exclusive mutation", "Shared mutation"), legend_gp=gpar(fill=c("#f15923", "#0068b3")))
legend2 <- Legend(title="SNV, fraction", at=c("SNV, PAC only (%)", "SNV, IDC-P only (%)", "Shared SNV (%)"), legend_gp=gpar(fill=c("#f15923", "#f7941d", "#0068b3")))
legend3 <- Legend(title="SCNA, fraction", at=c("SCNA, PAC only (%)", "SCNA, IDC-P only (%)", "Shared SNV (%)"), legend_gp=gpar(fill=c("#f15923", "#f7941d", "#0068b3")))
legend4 <- Legend(title="Evolution patterns", at=c("Indepentant origin", "Early divergent", "Late divergent"), legend_gp=gpar(fill=c("#dbf06a", "#abdda4", "#66c2a5")))
legend5 <- Legend(title="FOXA1 mutation", at=c("PAC only", "IDC−P only", "IDC−P and PAC", "None"), legend_gp=gpar(fill=c("#3288bd", "#3288bd", "#3288bd", "#e5e5e5")))
legend6 <- Legend(title="SCNA", at=c("Cluster seperate", "Cluster together"), legend_gp=gpar(fill=c("#7e72b4", "#72b47e")))
legend7 <- Legend(title="DNA methylation", at=c("Cluster seperate", "Cluster together", "Unavailable"), legend_gp=gpar(fill=c("#7e72b4", "#72b47e", "#e5e5e5")))
legend8 <- Legend(title="mRNA", at=c("Cluster seperate", "Cluster together", "Unavailable"), legend_gp=gpar(fill=c("#7e72b4", "#72b47e", "#e5e5e5")))
legend9 <- Legend(title="Cytological features", at=c("Resemblance", "Difference"), legend_gp=gpar(fill=c("#fc9272", "#ef3b2c")))
legend10 <- Legend(title="PFS ≤ 24 mo", at=c("No", "Yes"), legend_gp=gpar(fill=c("#e5e5e5", "black")))
legend_pack <- packLegend(legend1, legend2, legend3, legend4, legend5, legend6, legend7, legend8, legend9, legend10, direction="vertical", max_height=unit(110,"mm"))

# save to file
pdf(paste0(res_path, "0_single_bar_v2.pdf"), paper="USr", width=unit(279,"mm"), height=unit(216,"mm"))
draw(column_ha)
draw(legend_pack)
dev.off()

########################## top annotation: multiple barplot ##########################
# 1st layer barplot data
sample <- data.frame(unlist(sample_pre))
names(sample) <- "mutation_count"
sample$feature <- sapply(rownames(sample), function(x){gsub("[0-9]*", "", x)})
sample$patient <- sapply(rownames(sample_pre), function(x){gsub("_.*", "", x)})
sample$sample <- sapply(rownames(sample_pre), function(x){gsub(".*_", "", x)})
sample$feature <- factor(sample$feature, levels=unique(sample$feature))
sample$patient <- factor(sample$patient, levels=oncoplot_order)

# plot
stack_bar <- ggplot(sample, aes(x=sample, y=mutation_count, fill=feature)) +
	geom_bar(stat="identity", width=0.9, colour="white", lwd=0.75) +
	facet_grid(. ~ patient, scales="free_x", space="free") +
	scale_fill_manual(values=c(`Exclusive mutation`="#f15923", `shared mutation`="#0068b3")) +
	theme(legend.position="none",
		panel.background=element_blank(),
		panel.grid=element_blank(),
		panel.border=element_blank(),
		axis.title=element_blank(),
		axis.text=element_blank(),
		axis.ticks.x=element_blank())

pdf(paste0(res_path, "0_stack_bar.pdf"), width=unit(110,"mm"), height=unit(16.5,"mm"))
print(stack_bar)
dev.off()
