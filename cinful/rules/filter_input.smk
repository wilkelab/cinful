from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import seqhash

SAMPLES, = glob_wildcards("{sample}.fna")


rule filter_microcin:
	input:
		"cinfulOut/01_orf_homology/prodigal_out/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/microcins/all_samples/{sample}.filtered.fa"
	shell:
		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"

rule nonredundant_microcin:
	input:
		expand("cinfulOut/01_orf_homology/microcins/all_samples/{sample}.filtered.fa", sample=SAMPLES)
	output:
		fasta="cinfulOut/01_orf_homology/microcins/filtered_nr.fa",
		csv="cinfulOut/01_orf_homology/microcins/filtered_nr.csv"
	run:
		hashDict = {}
		idDict = {}
		for file in input:
			sample = file.split("cinfulOut/01_orf_homology/microcins/all_samples/")[1].strip(".filtered.fa")
			with open(file) as handle:
				for seq_record in SeqIO.parse(handle, "fasta"):
					pephash = seqhash.seqhash(seq_record.seq,dna_type='PROTEIN')
					sequence = str(seq_record.seq)
					
					descriptionParts = seq_record.description.split("#") 
					start = descriptionParts[1].strip()
					stop = descriptionParts[2].strip()
					strand = descriptionParts[3].strip()
					contig = '_'.join(seq_record.id.split("_")[:-1])
					hashDict[pephash] = sequence
					seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
					idDict[seqID] = [pephash,sample,contig,start,stop,strand,sequence]
		idDF = pd.DataFrame.from_dict(idDict,orient="index")
		idDF.columns = ["pephash","sample","contig","start","stop","strand","seq"]
									
		idDF.to_csv(output.csv)


		with open(output.fasta,"w") as fasta_file:
			for pephash in hashDict:
				outRecord = SeqRecord(
					Seq(hashDict[pephash]),
					id=pephash,
					description=""
				)
				SeqIO.write(outRecord, fasta_file, "fasta")
			



rule filter_immunity_protein:
	input:
		"cinfulOut/01_orf_homology/prodigal_out/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/immunity_proteins/all_samples/{sample}.filtered.fa"
	shell:
		"seqkit seq -m 30 -M 250  {input} | seqkit rmdup -s > {output}"

rule nonredundant_immunity_protein:
	input:
		expand("cinfulOut/01_orf_homology/immunity_proteins/all_samples/{sample}.filtered.fa", sample=SAMPLES)
	output:
		fasta="cinfulOut/01_orf_homology/immunity_proteins/filtered_nr.fa",
		csv="cinfulOut/01_orf_homology/immunity_proteins/filtered_nr.csv"
	run:
		hashDict = {}
		idDict = {}
		for file in input:
			sample = file.split("cinfulOut/01_orf_homology/immunity_proteins/all_samples/")[1].strip(".filtered.fa")
			with open(file) as handle:
				for seq_record in SeqIO.parse(handle, "fasta"):
					pephash = seqhash.seqhash(seq_record.seq,dna_type='PROTEIN')
					sequence = str(seq_record.seq)
					
					descriptionParts = seq_record.description.split("#") 
					start = descriptionParts[1].strip()
					stop = descriptionParts[2].strip()
					strand = descriptionParts[3].strip()
					contig = '_'.join(seq_record.id.split("_")[:-1])
					hashDict[pephash] = sequence
					seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
					idDict[seqID] = [pephash,sample,contig,start,stop,strand,sequence]
		idDF = pd.DataFrame.from_dict(idDict,orient="index")
		idDF.columns = ["pephash","sample","contig","start","stop","strand","seq"]
									
		idDF.to_csv(output.csv)


		with open(output.fasta,"w") as fasta_file:
			for pephash in hashDict:
				outRecord = SeqRecord(
					Seq(hashDict[pephash]),
					id=pephash,
					description=""
				)
				SeqIO.write(outRecord, fasta_file, "fasta")


rule filter_CvaB:
	input:
		"cinfulOut/01_orf_homology/prodigal_out/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/CvaB/all_samples/{sample}.filtered.fa"
	shell:
		"seqkit rmdup -s {input} > {output}"

rule nonredundant_CvaB:
	input:
		expand("cinfulOut/01_orf_homology/CvaB/all_samples/{sample}.filtered.fa", sample=SAMPLES)
	output:
		fasta="cinfulOut/01_orf_homology/CvaB/filtered_nr.fa",
		csv="cinfulOut/01_orf_homology/CvaB/filtered_nr.csv"
	run:
		hashDict = {}
		idDict = {}
		for file in input:
			sample = file.split("cinfulOut/01_orf_homology/CvaB/all_samples/")[1].strip(".filtered.fa")
			with open(file) as handle:
				for seq_record in SeqIO.parse(handle, "fasta"):
					pephash = seqhash.seqhash(seq_record.seq,dna_type='PROTEIN')
					sequence = str(seq_record.seq)

					descriptionParts = seq_record.description.split("#") 
					start = descriptionParts[1].strip()
					stop = descriptionParts[2].strip()
					strand = descriptionParts[3].strip()
					contig = '_'.join(seq_record.id.split("_")[:-1])
					hashDict[pephash] = sequence
					seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
					idDict[seqID] = [pephash,sample,contig,start,stop,strand,sequence]
		idDF = pd.DataFrame.from_dict(idDict,orient="index")
		idDF.columns = ["pephash","sample","contig","start","stop","strand","seq"]
									
		idDF.to_csv(output.csv)


		with open(output.fasta,"w") as fasta_file:
			for pephash in hashDict:
				outRecord = SeqRecord(
					Seq(hashDict[pephash]),
					id=pephash,
					description=""
				)
				SeqIO.write(outRecord, fasta_file, "fasta")



# TODO: add a merge nonredundant step for each component