# 1. Multiple seqeunce alignment using `MAFFT` v7.453

cd path/to/Endophyte-Fungus-Activity/mainFigure/scripts/

fa=../fasta/merged.fasta
mafft $fa > ${fa/.fasta/.aln}

# 2. Tree estimation using `iqtree` 1.6.12
fa=../fasta/merged.fasta
iqtree -s ${fa/.fasta/.aln}

# 3. Plot tree
Rscripts -e '
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

tree=read.tree("../fasta/merged.aln.treefile")

# Annotation data
dat=read.xlsx("../Fungi Major Clade ITS list for phylogenetic tree.xlsx")
names(dat)=c("label", "Species", "Class", "Phylum", "ID")
tree=full_join(tree, dat, by="label")
tree@phylo$tip.label=c("MG543732.1", "PQ269859.1", "KM513579.1", "MH712154.1", "MK243486.1", "OP238466.1", "OP060705.1", "MW016414.1", "MH712156.1", "JN624882.1", "MW581306.1", "ON979428.1", "MD-L1", "OQ102622.1", "MH497385.1", "MH497386.1", "MH497387.1", "MW581325.1", "MW581274.1", "MW581323.1", "HQ700354.1", "NR_198859.1", "OR939713.1", "HG995564.1", "KX442659.1", "MW600265.1", "MT066189.1", "OP801679.1", "MZ436987.1", "MD-L2", "MD-L4", "ON063274.1", "JX949177.1", "PQ302434.1", "MG980397.1", "JN030357.1", "OL413441.1", "KR816861.1", "PQ285209.1", "EU002891.1", "MD-T2", "EU250374.1", "JN545760.1", "NR_172743.1", "NR_154181.1", "KY965395.1", "MT159541.1", "MW186478.1", "JF705215.1", "KF176581.1", "MW369639.1", "KF176606.1", "MG837526.1", "MG837524.1", "MG837525.1", "NR_198123.1", "NR_185411.1", "NR_160640.1", "NR_158908.1", "NR_182807.1", "NR_177586.1", "NR_158795.1", "NR_158799.1", "NR_171823.1", "NR_175084.1", "NR_154671.1", "NR_158838.1", "NR_148062.1", "NR_155564.1", "NR_111082.1", "NR_111686.1", "NR_073292.1", "NR_191070.1", "NR_191069.1", "NR_198185.1", "NR_197471.1", "NR_191180.1", "NR_197799.1", "NR_103640.1", "NR_198092.1", "NR_191184.1", "NR_191183.1", "NR_191182.1", "NR_191181.1", "NR_198363.1", "NR_197470.1", "MD-R1", "NR_191201.1", "NR_198844.1", "NR_073343.1", "NR_189987.1", "NR_174711.1", "NR_174710.1", "NR_189760.1", "NR_174931.1", "NR_165942.1", "NR_185630.1", "NR_130693.1", "NR_182836.1", "NR_132806.1", "NR_119533.1", "NR_119625.1", "NR_198529.1", "NR_198528.1", "NR_189889.1", "NR_198524.1", "NR_198527.1", "NR_198525.1", "NR_198526.1", "NR_189888.1", "NR_198523.1", "NR_189887.1", "NR_198522.1", "NR_186908.1", "NR_185353.1", "NR_119547.1", "NR_182822.1", "NR_164211.1", "NR_147650.1", "NR_119552.1", "NR_121280.1", "NR_121478.1", "NR_121477.1", "NR_121476.1", "NR_189926.1", "NR_168761.1", "NR_168424.1", "NR_185673.1", "NR_168426.1", "NR_186946.1", "NR_186945.1", "NR_121448.1", "NR_189511.1", "NR_121469.1", "NR_121468.1", "NR_121199.1", "NR_197890.1", "NR_164223.1", "NR_164222.1", "NR_177181.1", "NR_164279.1", "NR_119555.1", "NR_165200.1", "NR_154392.1", "NR_154391.1", "NR_197859.1", "NR_111240.1", "NR_111241.1", "NR_159603.1", "NR_077175.1", "NR_198113.1", "NR_172188.1", "NR_197934.1", "NR_121206.1", "NR_172369.1", "NR_172284.1", "NR_160073.1", "NR_160070.1", "NR_160071.1", "NR_160065.1", "NR_164089.1", "NR_121541.1", "NR_182530.1", "NR_182528.1", "NR_198845.1", "NR_198339.1", "NR_198342.1", "NR_198341.1", "NR_198340.1", "NR_151797.1", "NR_198597.1", "FJ224123.1", "NR_138418.1", "NR_198848.1", "NR_198847.1", "MT913002.1", "KM520010.1", "MZ157058.1", "MK123443.1", "ON746571.1", "KJ625631.1", "PP537918.1", "OP288215.1", "MT993618.1", "MT928533.1", "MD-L3", "NR_160104.1", "NR_159869.1", "MN227420.1", "MN227440.1", "MN227405.1", "MN227407.1", "MN227446.1", "MN227478.1", "MN227422.1", "MN227428.1", "MN227423.1", "MN227410.1", "MD-T1", "NR_198858.1", "PQ699345.1", "MK281561.1", "KX925279.1", "FJ037739.1", "PQ302374.1", "GU066664.1", "KU671359.1", "GU066649.1", "MD-L5", "EU054404.1", "OM281089.1", "MD-H1", "MD-T10")
df=na.omit(as_tibble(tree))

n <- 38
colors <- hue_pal(h = c(15, 375), c = 100, l = 65)(n)
names(colors)=unique(df$Class)
colors["This study"]='black'

p=ggtree(tree, layout='circular', linetype=1) + 
geom_hilight(df, aes(node=node, fill=Phylum), type = "gradient", gradient.direction = 'rt', alpha = .8, to.bottom = TRUE) +
# geom_point(aes(x, y, shape=Phylum, color=Phylum), size=5) +
scale_shape_manual(values=seq(0,10)) +
geom_tree() + 
geom_tiplab(aes(color=Class), align=TRUE) +
scale_color_manual(values = colors) +
theme_tree()
ggsave("../figure/physeq_tree-distance.pdf", width=15, height=15)

p=ggtree(tree, layout='circular', branch.length='none', linetype=1) + 
geom_hilight(df, aes(node=node, fill=Phylum), type = "gradient", gradient.direction = 'rt', alpha = .8, to.bottom = TRUE) +
# geom_point(aes(x, y, shape=Phylum, color=Phylum), size=5) +
scale_shape_manual(values=seq(0,10)) +
geom_tree() + 
geom_tiplab(aes(color=Class), align=TRUE) +
scale_color_manual(values = colors) +
theme_tree()
ggsave("../figure/physeq_tree.pdf", width=15, height=15)
'

