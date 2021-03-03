SAMPLES, = glob_wildcards("{sample}.faa")
rule final:
	input:
		expand("{sample}_duomolog.blast_intersect_hmmer.fa", sample = SAMPLES)

rule filter_input:
	input:
		"{sample}.faa"
	output:
		"{sample}.30_150.fa"
	shell:
		"seqkit seq -M 150 -m 30 -M 150 {input} > {output}"

rule blast:
	input:
		verified_microcins = "verified_microcins.pep",
		input_seqs = "{sample}.30_150.fa"
	output:
		"{sample}.verified_microcins.blast.txt"
	shell:
		"""
		makeblastdb -dbtype prot -in {input.verified_microcins}
		blastp -db {input.verified_microcins} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1
		"""

rule verified_microcinsMSA:
	input:
		"verified_microcins.pep"
	output:
		"verified_microcins.aln"
	shell:
		"mafft --auto {input} > {output}"


rule buildhmm:
	input:
		"verified_microcins.aln"
	output:
		"verified_microcins.hmm"
	shell:
		"hmmbuild {output} {input}"


rule duomolog:
	input:
		verified_microcins = "verified_microcins.pep",
		input_seqs = "{sample}.30_150.fa",
		blastout="{sample}.verified_microcins.blast.txt",
		hmm="verified_microcins.hmm"
	output:
		"{sample}_duomolog.blast_intersect_hmmer.fa"
	shell:
		"""duomolog blast_v_hmmer -i {input.verified_microcins} -q {input.input_seqs} \
			--intersect_only \
			--blastout {input.blastout} \
			--hmm {input.hmm}	\
			-o {output}
		"""


rule hmmsearch:
	input:
		verified_microcinsHMM = "verified_microcins.hmm",
		input_seqs = "{sample}.30_150.fa"
	output:
		hmmerAlignment = "verified_microcins.hmmerAlignment.txt",
		hmmerResults = "verified_microcins.hmmerOut.txt"
	shell:
		"hmmsearch -hmm {input.verified_microcinsHMM} -tblout {output.hmmerResults} > {output.hmmerAlignment}"

rule getBestHits:
	input:
		blast_hits = "{sample}.verified_microcins.blast.txt",
		hmmer_hits = "verified_microcins.hmmerOut.txt",
		seq = "{sample}.30_150.fa"
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





		
