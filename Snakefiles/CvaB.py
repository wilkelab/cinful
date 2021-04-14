SAMPLES, = glob_wildcards("{sample}.faa")

rule final:
	input:
		expand("{sample}.cvaB_diamond.txt", sample = SAMPLES)

rule cvaB_diamond:
	input:
		cvaB = "cvaB.fa",
		db = "{sample}.faa"
	output:
		"{sample}.cvaB_diamond.txt"
	shell:"""
		diamond makedb --in {input.db} --out {input.db}.seq.db -d {input.db}
        diamond blastp -d {input.db}.dmnd -q {input.cvaB} -o {output} -p {threads} -e 1E-5
	"""
rule verified_CvaB:
	input:
		"verified_CvaB"
	output:
		"verified_CvaB.fa"
	shell:	
		"touch {output}"

rule verified_CvaBMSA:
	input:
		"verified_CvaB.fa"
	output:
		"verified_CvaB.aln"
	shell:
		"mafft --auto {input} > {output}"

rule filter_input:
	input:
		"input_seqs.fa"
	output:
		"input_seqs.short.fa"
	shell:
		"seqkit fx2tab -l {input} | awk '$3 < 150' > {output}"

rule blast:
	input:
		verified_CvaB = "verified_CvaB.fa",
		input_seqs = "input_seqs.short.fa"
	output:
		"verified_CvaB.blast.txt"
	shell:
		"""
		makeblastdb -in {input.verified_CvaB}
		blastp -db {input.verified_CvaB} -query {input.input_seqs} -outfmt 6 -out {output}
		"""

rule buildhmm:
	input:
		"verified_CvaB.aln"
	output:
		"verified_CvaB.hmm"
	shell:
		"hmmbuild {input} {output}"

rule hmmsearch:
	input:
		verified_CvaBHMM = "verified_CvaB.hmm",
		input_seqs = "input_seqs.short.fa"
	output:
		hmmerAlignment = "verified_CvaB.hmmerAlignment.txt",
		hmmerResults = "verified_CvaB.hmmerOut.txt"
	shell:
		"hmmsearch -hmm {input.verified_CvaBHMM} -tblout {output.hmmerResults} > {output.hmmerAlignment}"

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