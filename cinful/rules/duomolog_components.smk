import pandas as pd
import pyhmmer
from functools import reduce


SAMPLES, = glob_wildcards("cinfulOut/01_orf_homology/{sample}_prodigal")
COMPONENTS, = glob_wildcards("{component}.verified.pep")


def load_blast(blastout, hasHeaders=False):
	if hasHeaders:
		blastpDF = pd.read_csv(blastout, sep = "\t")
	else:
		blastpDF = pd.read_csv(blastout, sep = "\t", names = ["qseqid","sseqid","pident","length","mismatch",			"gapopen","qstart","qend","sstart","send","evalue","bitscore"])
	return blastpDF

def load_hmm(hmmFile):
    with pyhmmer.plan7.HMMFile(hmmFile) as h:
        hmm = next(h)
    return hmm

def run_hmmsearch(queryFile, hmmFile):
	hmm = load_hmm(hmmFile)
	with pyhmmer.easel.SequenceFile(queryFile) as seq_file:
		sequences = [ seq.digitize(hmm.alphabet) for seq in seq_file ]
	pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
	hits = pipeline.search_hmm(hmm, sequences) # Has lots of goodies!
	
	return hits , hmm.name.decode() # [hit.name.decode() for hit in hmmerOut]





rule merged_results:
	input:
		microcins = "cinfulOut/01_orf_homology/microcins/blast_v_hmmer.csv",
		immunity_proteins = "cinfulOut/01_orf_homology/immunity_proteins/blast_v_hmmer.csv",
		unfilteredCvaB = "cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.csv",
		QC_Cvab="cinfulOut/01_orf_homology/CvaB/QC.csv"
	output:
		"cinfulOut/02_homology_results/all_merged.csv"
	run:
		# componentDFs = []
		# for componentHomologFile in [input.microcins, input.immunity_proteins, input.CvaB]:
		# 	componentDFs.append(pd.read_csv(componentHomologFile))
		
		microcinDF = pd.read_csv(input.microcins)
		immunity_proteinDF = pd.read_csv(input.immunity_proteins)
		QC_CvabDF = pd.read_csv(input.QC_Cvab)
		unfilteredCvaBDF = pd.read_csv(input.unfilteredCvaB)
		# filteredCvaBDF = pd.read_csv(input.filteredCvaB)
		bestCvaB = unfilteredCvaBDF[unfilteredCvaBDF["qseqid"].isin(QC_CvabDF["id"]) ]
		
		componentDFs = [microcinDF, immunity_proteinDF, bestCvaB]
		mergedDf = pd.concat(componentDFs)
		mergedDf.to_csv(output[0], index = None)
		# with open(output[0],"w") as out:
			# out.write("test\n")
		






# rule hmmsearch:



# rule duomolog:
# 	input:
# 		verified_microcins = "verified_microcins.pep",
# 		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa",
# 		blastout="{sample}_cinfulOut/{sample}.verified_microcins.blast.txt",
# 		hmm="verified_microcins.hmm"
# 	output:
# 		"{sample}_cinfulOut/duomolog_microcin/summary_out.txt"
# 	shell:
# 		"""duomolog blast_v_hmmer --inFile {input.verified_microcins} --queryFile {input.input_seqs} \
# 			--blastFile {input.blastout} \
# 			--intersectOnly \
# 			--hmmFile {input.hmm}	\
# 			--summaryOut {output}
# 		"""

# rule duomolog:
# 	input:
# 		verified_immunity_proteins = "verified_immunity_proteins.pep",
# 		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa",
# 		blastout="{sample}_cinfulOut/{sample}.verified_immunity_proteins.blast.txt",
# 		hmm="verified_immunity_proteins.hmm"
# 	output:
# 		"{sample}_cinfulOut/duomolog_immunity_protein/summary_out.txt"
# 	shell:
# 		"""duomolog blast_v_hmmer --inFile {input.verified_immunity_proteins} --queryFile {input.input_seqs} \
# 			--blastFile {input.blastout} \
# 			--intersectOnly \
# 			--hmmFile {input.hmm}	\
# 			--summaryOut {output}
# 		"""	


# rule duomolog:
# 	input:
# 		verified_CvaB = "verified_CvaB.pep",
# 		input_seqs = "{sample}_cinfulOut/{sample}.faa",
# 		blastout = "{sample}_cinfulOut/{sample}.verified_CvaB.blast.txt",
# 		hmm = "verified_CvaB.hmm"
# 	output:
# 		"{sample}_cinfulOut/duomolog_CvaB/summary_out.txt"
# 	shell:
# 		"""duomolog blast_v_hmmer --inFile {input.verified_CvaB} --queryFile {input.input_seqs} \
# 			--blastFile {input.blastout} \
# 			--intersectOnly \
# 			--hmmFile {input.hmm}	\
# 			--summaryOut {output}
# 		"""
