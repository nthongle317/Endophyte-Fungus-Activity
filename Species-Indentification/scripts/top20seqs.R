library(msa)
library(stringr)
library(seqinr)
library(openxlsx)
library(ggtree)
library(ape)
library(treeio)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(phangorn)

options(ignore.negative.edge=TRUE)
path="/media/thong/sda/RCID/A.Phong/endophytic_fungi/Endophyte-Fungus-Activity/Species-Indentification/"
setwd(path)

for (sN in c("MD-H1", "MD-L1",  "MD-L2",  "MD-L3",  "MD-L4",  "MD-L5",  "MD-R1",  "MD-T1",  "MD-T10", "MD-T2")) {

	print(sN)
	outp=paste0("results/", sN)

	mySequenceFile=paste0("Blast/", sN, "/", sN, "_blast.fasta")
	mySequences=readAAStringSet(mySequenceFile)
	seqNames=mySequences@ranges@NAMES
	mySequences@ranges@NAMES=str_split_fixed(str_split_fixed(seqNames, " ", 2)[,1], ":", 2)[,1]

	metaData=read.xlsx(paste0(outp, "/BLAST-short.xlsx"))
	metaData=metaData[order(metaData$Per..ident, decreasing=T),]

	sl=mySequences@ranges@NAMES%in%c(metaData[1:20,5], sN)

	mySequences=mySequences[sl]

	writeXStringSet(mySequences, paste0("Blast/", sN, "/", sN, "_top20_blast.fasta"))
}