from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import seqhash
from csv import writer

SAMPLES, = glob_wildcards("{sample}.fna")
if SAMPLES == []:
	SAMPLES, = glob_wildcards("{config[outdir]}/01_orf_homology/prodigal_out/{sample}.faa")


def hasAllStandardAA(seq, alphabet="ACDEFGHIKLMNPQRSTVWY",ignore="*"):
	return (set(seq) - set(alphabet+ignore)) == set()


rule nonredundant_prodigal:
	input:
		expand(config["outdir"]+"/01_orf_homology/prodigal_out/{sample}.faa", sample=SAMPLES)
	output:
		fasta=config["outdir"] + "/01_orf_homology/prodigal_out.all.nr.faa",
		csv=config["outdir"] + "/01_orf_homology/prodigal_out.all.nr_expanded.csv"
	threads:threads_max
	message:
		"Combining Prodigal outputs into fasta and csv files."
	run:
		idDict = {}
		hashDict = {}
		csvDF = pd.DataFrame.from_dict(idDict, orient="index", columns=["pephash","sample","contig","start","stop","strand","allStandardAA","startCodon","RBSmotif","gcContent","seq"]).reset_index()
		csvDF.rename(columns={'index': 'cinful_id'}, inplace = True)
		csvDF.to_csv(output.csv, index = None)
		fastaDF = pd.DataFrame()
		fastaDF.to_csv(output.fasta, index=False, header=0)
		with open(output.csv, "a") as csv_file, open(output.fasta, "w") as fasta_file:
			for file in input:
					sample = file.split("01_orf_homology/prodigal_out/")[1].strip(".faa")
					with open(file) as handle:
						for seq_record in SeqIO.parse(handle, "fasta"):
							sequence = str(seq_record.seq)
							pephash = seqhash.seqhash(sequence.strip("*"),dna_type='PROTEIN')
							hashDict[pephash] = sequence

							descriptionParts = seq_record.description.split("#")
							gc_content = seq_record.description.split("=")[-1]
							start_codon = seq_record.description.split(";")[2].strip("start_type=")
							rbs_motif = seq_record.description.split(";")[3].strip("rbs_motif=")
							start = descriptionParts[1].strip()
							stop = descriptionParts[2].strip()
							strand = descriptionParts[3].strip()
							contig = '_'.join(seq_record.id.split("_")[:-1])
							allStandardAA = hasAllStandardAA(sequence)
							
							seqID = f"{sample}|{contig}|{start}:{stop}:{strand}"
							csv_input = [seqID, pephash, sample, contig, start, stop, strand, allStandardAA, start_codon, rbs_motif, gc_content, sequence]
							csv_writer = writer(csv_file)
							csv_writer.writerow(csv_input)
			for pephash in hashDict:
				outRecord = SeqRecord(
					Seq(hashDict[pephash]),
					id=pephash,
					description=""
				)
				SeqIO.write(outRecord, fasta_file, "fasta")


rule filter_microcin:
	input:
		config["outdir"] + "/01_orf_homology/prodigal_out.all.nr.faa"
	output:
		config["outdir"] + "/01_orf_homology/microcins/filtered_nr.fa"
	shell:
		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"

rule filter_immunity_protein:
	input:
		config["outdir"] + "/01_orf_homology/prodigal_out.all.nr.faa"
	output:
		config["outdir"] + "/01_orf_homology/immunity_proteins/filtered_nr.fa"
	shell:
		"seqkit seq -m 30 -M 250  {input} | seqkit rmdup -s > {output}"

rule filter_CvaB:
	input:
		config["outdir"] + "/01_orf_homology/prodigal_out.all.nr.faa"
	output:
		config["outdir"] + "/01_orf_homology/CvaB/filtered_nr.fa"
	shell:
		"seqkit seq -m 600 -M 800 {input} | seqkit rmdup -s > {output}"

rule filter_MFP:
	input:
		config["outdir"] + "/01_orf_homology/prodigal_out.all.nr.faa"
	output:
		config["outdir"] + "/01_orf_homology/MFP/filtered_nr.fa"
	shell:
		"seqkit seq -m 375 -M 450 {input} | seqkit rmdup -s > {output}"
