path="path/to/Endophyte-Fungus-Activity/Species-Indentification/scripts"
# path="/media/thong/sda/RCID/A.Phong/endophytic_fungi/Endophyte-Fungus-Activity/Species-Indentification/scripts"

cd $path 
for sN in MD-H1 MD-L1 MD-L2 MD-L3 MD-L4 MD-L5 MD-R1 MD-T1 MD-T10 MD-T2; do
	outp="../results/$sN"
	fa="../Blast/$sN/${sN}_top20_blast.fasta"

	# Sequence alignment
	mafft $fa > $outp/${sN}_top20_blast.aln
	# phylogenetic tree
	iqtree -s $outp/${sN}_top20_blast.aln

done


