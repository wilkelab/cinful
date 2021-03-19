SAMPLES, = glob_wildcards("{sample}.faa")
rule final:
	input:
		expand("{sample}_duomolog_immunity_protein/summary_out.txt", sample = SAMPLES)

rule filter_input:
	input:
		"{sample}.faa"
	output:
		"{sample}.30_250.fa"
	shell:
		"seqkit seq -m 30 -M 250  {input} | seqkit rmdup -s > {output}"


rule makeblastdb:
	input:
		"verified_immunity_proteins.pep"
	output:
		"verified_immunity_proteins.pep.phr"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast:
	input:
		verified_immunity_proteins = "verified_immunity_proteins.pep",
		blastdb = "verified_immunity_proteins.pep.phr",
		input_seqs = "{sample}.30_150.fa"
	output:
		"{sample}.verified_immunity_proteins.blast.txt"
	shell:
		"blastp -db {input.verified_immunity_proteins} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

rule verified_immunity_proteinsMSA:
	input:
		"verified_immunity_proteins.pep"
	output:
		"verified_immunity_proteins.aln"
	shell:
		"mafft --auto {input} > {output}"

rule buildhmm:
	input:
		"verified_immunity_proteins.aln"
	output:
		"verified_immunity_proteins.hmm"
	shell:
		"hmmbuild {output} {input}"

rule duomolog:
	input:
		verified_immunity_proteins = "verified_immunity_proteins.pep",
		input_seqs = "{sample}.30_150.fa",
		blastout="{sample}.verified_immunity_proteins.blast.txt",
		hmm="verified_immunity_proteins.hmm"
	output:
		"{sample}_duomolog_immunity_protein/summary_out.txt"
	shell:
		"""duomolog blast_v_hmmer --inFile {input.verified_immunity_proteins} --queryFile {input.input_seqs} \
			--blastFile {input.blastout} \
			--hmmFile {input.hmm}	\
			--summaryOut {output}
		"""		





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