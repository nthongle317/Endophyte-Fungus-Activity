library(ggtree)
library(ape)
library(treeio)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(scales)

path="path/to/Endophyte-Fungus-Activity/Species-Indentification/scripts"
# path="/media/thong/sda/RCID/A.Phong/endophytic_fungi/Endophyte-Fungus-Activity/Species-Indentification/scripts"
setwd(path)

for (sN in c("MD-H1", "MD-L1",  "MD-L2",  "MD-L3",  "MD-L4",  "MD-L5",  "MD-R1",  "MD-T1",  "MD-T10", "MD-T2")) {

tree=read.tree(paste0("../results/", sN, "/MD-H1_top20_blast.aln.treefile"))

temp1=data.frame(seqN=as.phylo(tree)$tip.label)
metaData=read.xlsx(paste0("../results/", sN, "/BLAST-short.xlsx"))
metaData=metaData[order(metaData$Per..ident, decreasing=T),]
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
ggsave(paste0(outp, "/", sN, "_tree.pdf"), height=7, width=10)

}