library(ggplot2)
library(readxl)
library(openxlsx)
library(dplyr)

path="path/to/Endophyte-Fungus-Activity/Species-Indentification/scripts"
# path="/media/thong/sda/RCID/A.Phong/endophytic_fungi/Endophyte-Fungus-Activity/Species-Indentification/scripts"
setwd(path)

for (sN in c("MD-H1", "MD-L1",  "MD-L2",  "MD-L3",  "MD-L4",  "MD-L5",  "MD-R1",  "MD-T1",  "MD-T10", "MD-T2")) {

	print(sN)

	f=list.files(paste0("../Blast/", sN, "/"), pattern="*-Alignment-Descriptions.csv", full.names = T)

	dt=read.table(f, header=T, sep=",", check.names = F, colClasses=c("NULL", rep(NA, 8)))
	outDT=dt[,c(1,4,5,6,8)]

	df=data.frame(outDT %>%
	group_by(`Scientific Name`) %>%
	summarise(nIdent = n(), meanIdent= mean(`Per. ident`)))

	df$Scientific.Name=factor(df$Scientific.Name, levels=df[order(df[,2]),"Scientific.Name"])

	names(df)[3]="Identity (%)"
	ggplot(df, aes(x=Scientific.Name , y=`nIdent`, fill=`Identity (%)`)) +
	geom_bar(stat="identity", alpha=.6, width=.4) +
	coord_flip() +
	xlab("") +
	ylab("Number of BLAST hits") +
	#guides(fill=guide_legend(title="Độ tương đồng (%)")) +
	#scale_y_continuous() +
	theme_light(base_size=20)

	outp=paste0("../results/", sN)
	dir.create(outp, recursive = T)

	ggsave(paste0(outp, "/Taxonomy.png"), bg="white", height = 10, width = 10)

	write.xlsx(outDT, (paste0(outp, "/BLAST-short.xlsx")))

}