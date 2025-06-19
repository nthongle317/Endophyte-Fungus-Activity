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

path="/media/thong/sda/RCID/A.Phong/10 vi nấm nội sinh từ cây hương phụ -Ngân/thong/"
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

	myFirstAlignment=msa(mySequences)
	myFirstAlignment

	msaPrettyPrint(myFirstAlignment, shadingMode="functional", output="pdf", shadingModeArg="accessible area", logoColors="accessible area", askForOverwrite=FALSE)

	file.rename(from = "myFirstAlignment.pdf",  to = paste0(outp, "/", sN, "_MSA.pdf"))

	myFirstAlignment1=msaConvert(myFirstAlignment, type="seqinr::alignment")
	d=dist.alignment(myFirstAlignment1, "identity")
	tree=nj(d)

	temp1=data.frame(seqN=as.phylo(tree)$tip.label)
	df=merge(temp1, metaData[1:20,], by.x="seqN", by.y="Accession", all=T)
	df[df[,1]==sN, 2]="This study"

	x=tibble(label=df[,1], Name=df[,2])
	tree=full_join(tree, x, by="label")
	df=na.omit(as_tibble(tree))
	# df$Name[[1]]=`Fungal endophyte`
	# tree@ extraInfo$Name[[22]]=`Fungal endophyte`

	p=ggplot(tree, aes(x, y), branch.length='none') + 
	geom_hilight(df, aes(node=node, fill=Name), align="both", type = "gradient", gradient.direction = 'rt', alpha = .8) +
	# geom_point(aes(x, y, color=Name), size=5) +
	geom_tree() + 
	geom_tiplab(aes(color=Name)) +
	theme_tree()
	ggsave(paste0(outp, "/", sN, "_nj_tree.pdf"), height=7, width=10)

	tree=upgma(d)

	temp1=data.frame(seqN=as.phylo(tree)$tip.label)
	df=merge(temp1, metaData[1:20,], by.x="seqN", by.y="Accession", all=T)
	df[df[,1]==sN, 2]="This study"

	x=tibble(label=df[,1], Name=df[,2])
	tree=full_join(tree, x, by="label")
	df=na.omit(as_tibble(tree))

	p=ggplot(tree, aes(x, y), branch.length='none') + 
	geom_hilight(df, aes(node=node, fill=Name), align="both", type = "gradient", gradient.direction = 'rt', alpha = .8) +
	# geom_point(aes(x, y, color=Name), size=5) +
	geom_tree() + 
	geom_tiplab(aes(color=Name)) +
	theme_tree()
	ggsave(paste0(outp, "/", sN, "_upgma_tree.pdf"), height=7, width=10)


}