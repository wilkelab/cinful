import pandas as pd
from Bio import SeqIO

SAMPLES, = glob_wildcards("{sample}.fna")

def prodigalFa2DF(fa):
  prodigalDict = {"id":[],"seq":[],"contig":[],"start":[],"stop":[],"strand":[]}
  for seq_record in SeqIO.parse(fa,"fasta"):
    descriptionParts = seq_record.description.split("#") 

    start = descriptionParts[1].strip()
    stop = descriptionParts[2].strip()
    strand = descriptionParts[3].strip()
    contig = '_'.join(seq_record.id.split("_")[:-1])

    prodigalDict["id"].append(seq_record.id)
    prodigalDict["seq"].append(seq_record.seq)
    prodigalDict["contig"].append(contig)
    prodigalDict["stop"].append(stop)
    prodigalDict["start"].append(start)
    prodigalDict["strand"].append(strand)

  return pd.DataFrame.from_dict(prodigalDict)


def homology_withProdigal(homology_results_file, prodigalDF):
	homology_resultsDF = pd.read_csv(homology_results_file)

	return homology_resultsDF.merge(prodigalDF, left_on="qseqid",right_on="pephash")#, left_on = "qseqid", right_on = "id")

def componentDFs(homology_withProdigalDF):	
	microcinDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="microcins.verified"]
	immunity_proteinDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="immunity_proteins.verified"]
	CvaBDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="CvaB.verified"]
	return microcinDF, immunity_proteinDF, CvaBDF

def bestHits(microcinDF,immunity_proteinDF,CvaBDF):
  best_microcin_idx = microcinDF.groupby(['contig'])['bitscore'].transform(max) == microcinDF['bitscore']
  best_microcinDF = microcinDF[best_microcin_idx]

  best_immunity_protein_idx = immunity_proteinDF.groupby(['contig'])['bitscore'].transform(max) == immunity_proteinDF['bitscore']
  best_immunity_proteinDF = immunity_proteinDF[best_immunity_protein_idx]

  best_CvaB_idx = CvaBDF.groupby(['contig'])['bitscore'].transform(max) == CvaBDF['bitscore']
  best_CvaBDF = CvaBDF[best_CvaB_idx]

  return pd.concat([best_microcinDF, best_immunity_proteinDF,best_CvaBDF])




rule best_hits:
	input:
		merged_homology_results = "cinfulOut/02_homology_results/all_merged.csv", 
		nr_csv = "cinfulOut/01_orf_homology/CvaB/filtered_nr.csv"
	output:
		"cinfulOut/03_best_hits/best_hits.csv"
	run:
		print("merged_homology_results:",input.merged_homology_results)
		# for inFileIndex in range(len(input.merged_homology_results)):
		homologyFile = input.merged_homology_results #[inFileIndex]
		########### will break here ##########
		# prodigalPepFile = input.nr_csv #[inFileIndex] 
		prodigalDF =  pd.read_csv(input.nr_csv) #prodigalFa2DF(prodigalPepFile)
		######################################
		
		homology_withProdigalDF = homology_withProdigal(homologyFile, prodigalDF)
		microcinDF, immunity_proteinDF, CvaBDF = componentDFs(homology_withProdigalDF)
		if not microcinDF.empty:
			best_hitsDF = bestHits(microcinDF,immunity_proteinDF,CvaBDF)
			best_hitsDF.to_csv(output[0], index = False)

# rule best_hits_combined:
# 	input:
# 		"cinfulOut/02_homology_results/best_hits.csv", sample = SAMPLES)
# 	output:
# 		"cinfulOut/03_best_hits/best_hits.csv"
# 	run:
# 		inDFs = []
		
# 		for inFile in input:
# 			sample = inFile.replace("cinfulOut/02_homology_results/","").replace(".best_hits.csv","")
# 			# print(f"################\n{sample}\n###############")
# 			inDF = pd.read_csv(inFile)
# 			inDF.insert(0,"sample",sample)
# 			inDFs.append(inDF)
# 		outDF = pd.concat(inDFs)
# 		outDF.to_csv(output[0], index=False)

