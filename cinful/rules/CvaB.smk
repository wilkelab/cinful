SAMPLES, = glob_wildcards("{sample}.faa")

rule final:
	input:
		expand("{sample}_duomolog_CvaB/summary_out.txt", sample = SAMPLES)




# rule filter_input:
# 	input:
# 		"input_seqs.fa"
# 	output:
# 		"input_seqs.short.fa"
# 	shell:
# 		"seqkit fx2tab -l {input} | awk '$3 < 150' > {output}"

rule makeblastdb:
	input:
		"verified_CvaB.fa"
	output:
		"verified_CvaB.fa.phr"		
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast:
	input:
		verified_CvaB = "verified_CvaB.fa",
		blastdb = "verified_CvaB.pep.phr",
		input_seqs = "{sample}.faa"
	output:
		"{sample}.verified_CvaB.blast.txt"
	shell:
		"blastp -db {input.verified_CvaB} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"
		


rule verified_CvaBMSA:
	input:
		"verified_CvaB.fa"
	output:
		"verified_CvaB.aln"
	shell:
		"mafft --auto {input} > {output}"
		
rule buildhmm:
	input:
		"verified_CvaB.aln"
	output:
		"verified_CvaB.hmm"
	shell:
		"hmmbuild {input} {output}"

rule duomolog:
	input:
		verified_CvaB = "verified_CvaB.fa",
		input_seqs = "{sample}.faa",
		blastout = "{sample}.verified_CvaB.blast.txt",
		hmm = "verified_CvaB.hmm"
	output:
		"{sample}_duomolog_CvaB/summary_out.txt"
	shell:
		"""duomolog blast_v_hmmer --inFile {input.verified_CvaB} --queryFile {input.input_seqs} \
			--blastFile {input.blastout} \
			--intersectOnly \
			--hmmFile {input.hmm}	\
			--summaryOut {output}
		"""

rule getBestHits:
	input:
		blast_hits = "verified_CvaB.blast.txt",
		hmmer_hits = "verified_CvaB.hmmerOut.txt",
		seq = "input_seqs.short.fa"
	output:
		"verified_CvaB.bestHits.fa"
	shell:
		"touch {output}"

rule bestHitsMSA:
	input:
		"verified_CvaB.bestHits.fa"
	output:
		"verified_CvaB.bestHits.aln"
	shell:
		"touch {output}"

rule evaluateMSA:
	input:
		"verified_CvaB.bestHits.aln"
	output:
		"verified_CvaB.evaluateMSA.txt"

rule catalytic_triad:
	input:
		"verified_CvaB.evaluateMSA.txt"
	output:
		"verified_CvaB.catalytic_triad.txt"


rule putative_CvaB:
	input:
		seqs = "verified_CvaB.bestHits.fa",
		evaluateMSA = "verified_CvaB.evaluateMSA.txt",
		catalytic_triad = "verified_CvaB.catalytic_triad.txt"

	output:
		"putative_CvaB.txt"