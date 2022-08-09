import pandas as pd
import pyhmmer
from functools import reduce


SAMPLES, = glob_wildcards("{config[outdir]}/01_orf_homology/{sample}_prodigal")
COMPONENTS, = glob_wildcards("{component}.verified.pep")


def load_blast(blastout, hasHeaders=False):
	if hasHeaders:
		blastpDF = pd.read_csv(blastout, sep = "\t")
	else:
		blastpDF = pd.read_csv(blastout, sep = "\t", names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
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


rule makeblastdb:
	input:
		"{component}.verified.pep"
	output:
		"{component}.verified.pep.phr"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast:
	input:
		verified_component = "{component}.verified.pep",
		blastdb = "{component}.verified.pep.phr",
		input_seqs = "{config[outdir]}/01_orf_homology/{sample}_prodigal/{component}/{sample}.filtered.fa"
	output:
		"{config[outdir]}/01_orf_homology/{sample}_prodigal/{component}/blast.txt"
	shell:
		"blastp -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

rule msa:
	input:
		"{component}.verified.pep"
	output:
		"{component}.verified.aln"
	shell:
		"mafft --auto {input} > {output}"

rule buildhmm:
	input:
		"{component}.verified.aln"
	output:
		"{component}.verified.hmm"
	shell:
		"hmmbuild {output} {input}"

rule blast_v_hmmer:
	input:
		verifiedHMM = "{component}.verified.hmm",
		input_seqs = "{config[outdir]}/01_orf_homology/{sample}_prodigal/{component}/{sample}.filtered.fa",
		blastOut = "{config[outdir]}/01_orf_homology/{sample}_prodigal/{component}/blast.txt"
	output:
		"{config[outdir]}/01_orf_homology/{sample}_prodigal/{component}/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		# print(blastDF.shape)
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)#hmmer_hitsHeaders in blastDF["qseqid"]
		blastDF.to_csv(output[0], index = False)

rule merged_results:
	input:
		blast_v_hmmer = expand("{config[outdir]}/01_orf_homology/{sample}_prodigal/{component}/blast_v_hmmer.csv",
		                       sample=SAMPLES, component=COMPONENTS),
		prodigalGFF = "{config[outdir]}/01_orf_homology/{sample}_prodigal/{sample}.gff3"
	output:
		"{config[outdir]}/02_homology_results/{sample}.all_merged.csv"
	run:
		# print(len(input.blast_v_hmmer), input.blast_v_hmmer)
		componentDFs = []
		for componentHomologFile in input.blast_v_hmmer:
			componentDFs.append(pd.read_csv(componentHomologFile))
		mergedDf = pd.concat(componentDFs)
		mergedDf.to_csv(output[0], index = None)
		# with open(output[0],"w") as out:
			# out.write("test\n")
