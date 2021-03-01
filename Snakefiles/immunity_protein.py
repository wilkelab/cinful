rule final:
	input:
		"putative_immunity_proteins.txt"

rule input_seqs:
	input:
		"input_seqs"
	output:
		"input_seqs.fa"
	shell:
		"touch {output}"

rule verified_immunity_proteins:
	input:
		"verified_immunity_proteins"
	output:
		"verified_immunity_proteins.fa"
	shell:	
		"touch {output}"

rule verified_immunity_proteinsMSA:
	input:
		"verified_immunity_proteins.fa"
	output:
		"verified_immunity_proteins.aln"
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
		verified_immunity_proteins = "verified_immunity_proteins.fa",
		input_seqs = "input_seqs.short.fa"
	output:
		"verified_immunity_proteins.blast.txt"
	shell:
		"""
		makeblastdb -in {input.verified_immunity_proteins}
		blastp -db {input.verified_immunity_proteins} -query {input.input_seqs} -outfmt 6 -out {output}
		"""

rule buildhmm:
	input:
		"verified_immunity_proteins.aln"
	output:
		"verified_immunity_proteins.hmm"
	shell:
		"hmmbuild {input} {output}"

rule hmmsearch:
	input:
		verified_immunity_proteinsHMM = "verified_immunity_proteins.hmm",
		input_seqs = "input_seqs.short.fa"
	output:
		hmmerAlignment = "verified_immunity_proteins.hmmerAlignment.txt",
		hmmerResults = "verified_immunity_proteins.hmmerOut.txt"
	shell:
		"hmmsearch -hmm {input.verified_immunity_proteinsHMM} -tblout {output.hmmerResults} > {output.hmmerAlignment}"

rule getBestHits:
	input:
		blast_hits = "verified_immunity_proteins.blast.txt",
		hmmer_hits = "verified_immunity_proteins.hmmerOut.txt",
		seq = "input_seqs.short.fa"
	output:
		"verified_immunity_proteins.bestHits.fa"
	shell:
		"touch {output}"

rule bestHitsMSA:
	input:
		"verified_immunity_proteins.bestHits.fa"
	output:
		"verified_immunity_proteins.bestHits.aln"
	shell:
		"touch {output}"

rule evaluateMSA:
	input:
		"verified_immunity_proteins.bestHits.aln"
	output:
		"verified_immunity_proteins.evaluateMSA.txt"

rule subcellular_localization:
	input:
		"verified_immunity_proteins.bestHits.fa"
	output:
		"verified_immunity_proteins.bestHits.subcellular_localization.txt"
	shell:
		"touch {output}"

rule transmembrane_helix:
	input:
		"verified_immunity_proteins.bestHits.fa"
	output:
		"verified_immunity_proteins.bestHits.transmembrane_helix.txt"
	shell:
		"touch {output}"



rule putative_immunity_proteins:
	input:
		seqs = "verified_immunity_proteins.bestHits.fa",
		evaluateMSA = "verified_immunity_proteins.evaluateMSA.txt",
		subcellular_localization = "verified_immunity_proteins.bestHits.subcellular_localization.txt",
		transmembrane_helix = "verified_immunity_proteins.bestHits.transmembrane_helix.txt"

	output:
		"putative_immunity_proteins.txt"