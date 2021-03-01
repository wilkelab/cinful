rule final:
	input:
		"putative_microcins.txt"


rule input_seqs:
	input:
		"input_seqs"
	output:
		"input_seqs.fa"
	shell:
		"touch {output}"

rule verified_microcins:
	input:
		"verified_microcins"
	output:
		"verified_microcins.fa"
	shell:	
		"touch {output}"

rule filter_input:
	input:
		"input_seqs.fa"
	output:
		"input_seqs.short.fa"
	shell:
		"seqkit fx2tab -l {input} | awk '$3 < 150' > {output}"

rule blast:
	input:
		verified_microcins = "verified_microcins.fa",
		input_seqs = "input_seqs.short.fa"
	output:
		"verified_microcins.blast.txt"
	shell:
		"""
		makeblastdb -in {input.verified_microcins}
		blastp -db {input.verified_microcins} -query {input.input_seqs} -outfmt 6 -out {output}
		"""

rule verified_microcinsMSA:
	input:
		"verified_microcins.fa"
	output:
		"verified_microcins.aln"


rule buildhmm:
	input:
		"verified_microcins.aln"
	output:
		"verified_microcins.hmm"
	shell:
		"hmmbuild {input} {output}"

rule hmmsearch:
	input:
		verified_microcinsHMM = "verified_microcins.hmm",
		input_seqs = "input_seqs.short.fa"
	output:
		hmmerAlignment = "verified_microcins.hmmerAlignment.txt",
		hmmerResults = "verified_microcins.hmmerOut.txt"
	shell:
		"hmmsearch -hmm {input.verified_microcinsHMM} -tblout {output.hmmerResults} > {output.hmmerAlignment}"

rule getBestHits:
	input:
		blast_hits = "verified_microcins.blast.txt",
		hmmer_hits = "verified_microcins.hmmerOut.txt",
		seq = "input_seqs.short.fa"
	output:
		"verified_microcins.bestHits.fa"
	shell:
		"touch {output}"

rule bestHitsMSA:
	input:
		"verified_microcins.bestHits.fa"
	output:
		"verified_microcins.bestHits.aln"
	shell:
		"touch {output}"

rule evaluateMSA:
	input:
		"verified_microcins.bestHits.aln"
	output:
		"verified_microcins.evaluateMSA.txt"
rule checkSignalSeq:
	input:
		"verified_microcins.bestHits.aln"
	output:
		"verified_microcins.signalCheck.txt"

rule neuBI:
	input:
		"verified_microcins.bestHits.fa"
	output:
		"verified_microcins.bestHits.neuBI.txt"
	shell:
		"neuBI {input} {output}"

rule putative_microcins:
	input:
		seqs = "verified_microcins.bestHits.fa",
		evaluateMSA = "verified_microcins.evaluateMSA.txt",
		signalCheck = "verified_microcins.signalCheck.txt",
		neuBI = "verified_microcins.bestHits.neuBI.txt"
	output:
		"putative_microcins.txt"





		
