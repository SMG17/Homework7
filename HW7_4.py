import gffutils
import pysam
import numpy as np
import scipy.stats as st

db = gffutils.FeatureDB('yeast.db')
bamFileBY = pysam.AlignmentFile("outputBY.sorted.bam","rb")
bamFileRM = pysam.AlignmentFile("outputRM.sorted.bam","rb")
output = open("HW7_4_data.txt","w")

GeneNames = []
k_BY = []
k_RM = []
N_BY = 0.0
N_RM = 0.0
LengthPerGene = []
output.write("gene_name\tgene_length\tBY_expression\tRM_expression\n")
for mRNA in db.features_of_type("mRNA"):
	if mRNA.chrom == "chrmt":
		continue
	countBY = bamFileBY.count(mRNA.chrom,mRNA.start,mRNA.stop) #k11 or countBY
	countRM = bamFileRM.count(mRNA.chrom,mRNA.start,mRNA.stop) #k21 or countRM
	mRNA_length = 0.0
	for CDS in db.children(mRNA, featuretype = "CDS"):
		CDS_length = CDS.stop - CDS.start + 1
		mRNA_length += CDS_length
	output.write("%s\t%d\t%d\t%d\n" % (mRNA["Name"],mRNA_length,countBY,countRM))

output.close()
