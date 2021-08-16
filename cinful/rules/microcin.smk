from io import StringIO
from Bio import SeqIO


SAMPLES, = glob_wildcards("cinfulOut/01_orf_homology/{sample}_prodigal/")

# print(SAMPLES)
# rule final:
# 	# input:
	# 	expand("cinfulOut/01_orf_homology/{sample}_prodigal/microcins/{sample}.filtered.fa", sample = SAMPLES)

# rule filter_microcin:
# 	input:
# 		"cinfulOut/01_orf_homology/{sample}_prodigal/{sample}.faa"
# 	output:
# 		"cinfulOut/01_orf_homology/{sample}_prodigal/microcins/{sample}.filtered.fa"
# 	shell:
# 		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"


rule makeblastdb_microcin:
	input:
		"cinfulOut/00_dbs/microcins.verified.pep"
	output:
		"cinfulOut/00_dbs/microcins.verified.pep.phr"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast_microcin:
	input:
		verified_component = "cinfulOut/00_dbs/microcins.verified.pep",
		blastdb = "cinfulOut/00_dbs/microcins.verified.pep.phr",
		input_seqs = "cinfulOut/01_orf_homology/{sample}_prodigal/microcins/{sample}.filtered.fa"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/microcins/blast.txt"
	shell:
		"blastp -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

rule msa_microcin:
	input:
		"cinfulOut/00_dbs/microcins.verified.pep"
	output:
		"cinfulOut/00_dbs/microcins.verified.aln"
	shell:
		"mafft --auto {input} > {output}"

rule buildhmm_microcin:
	input:
		"cinfulOut/00_dbs/microcins.verified.aln"
	output:
		"cinfulOut/00_dbs/microcins.verified.hmm"
	shell:
		"hmmbuild {output} {input}"



rule blast_v_hmmer_microcin:
	input:
		verifiedHMM = "cinfulOut/00_dbs/microcins.verified.hmm",
		input_seqs = "cinfulOut/01_orf_homology/{sample}_prodigal/microcins/{sample}.filtered.fa",
		blastOut = "cinfulOut/01_orf_homology/{sample}_prodigal/microcins/blast.txt"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/microcins/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)#hmmer_hitsHeaders in blastDF["qseqid"]
		blastDF.to_csv(output[0], index = False)


# rule makeblastdb:
# 	input:
# 		"verified_microcins.pep"
# 	output:
# 		"verified_microcins.pep.phr"
# 	shell:
# 		"makeblastdb -dbtype prot -in {input}"

# rule blast:
# 	input:
# 		verified_microcins = "verified_microcins.pep",
# 		blastdb = "verified_microcins.pep.phr",
# 		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa"
# 	output:
# 		"{sample}_cinfulOut/{sample}.verified_microcins.blast.txt"
# 	shell:
# 		"blastp -db {input.verified_microcins} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

# rule verified_microcinsMSA:
# 	input:
# 		"verified_microcins.pep"
# 	output:
# 		"verified_microcins.aln"
# 	shell:
# 		"mafft --auto {input} > {output}"


# rule buildhmm:
# 	input:
# 		"verified_microcins.aln"
# 	output:
# 		"verified_microcins.hmm"
# 	shell:
# 		"hmmbuild {output} {input}"


# rule duomolog:
# 	input:
# 		verified_microcins = "verified_microcins.pep",
# 		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa",
# 		blastout="{sample}_cinfulOut/{sample}.verified_microcins.blast.txt",
# 		hmm="verified_microcins.hmm"
# 	output:
# 		"{sample}_cinfulOut/duomolog_microcin/summary_out.txt"
# 	shell:
# 		"""duomolog blast_v_hmmer --inFile {input.verified_microcins} --queryFile {input.input_seqs} \
# 			--blastFile {input.blastout} \
# 			--intersectOnly \
# 			--hmmFile {input.hmm}	\
# 			--summaryOut {output}
# 		"""


# rule hmmsearch:
# 	input:
# 		verified_microcinsHMM = "verified_microcins.hmm",
# 		input_seqs = "{sample}.30_150.fa"
# 	output:
# 		hmmerAlignment = "verified_microcins.hmmerAlignment.txt",
# 		hmmerResults = "verified_microcins.hmmerOut.txt"
# 	shell:
# 		"hmmsearch -hmm {input.verified_microcinsHMM} -tblout {output.hmmerResults} > {output.hmmerAlignment}"

# rule getBestHits:
# 	input:
# 		blast_hits = "{sample}.verified_microcins.blast.txt",
# 		hmmer_hits = "verified_microcins.hmmerOut.txt",
# 		seq = "{sample}.30_150.fa"
# 	output:
# 		"verified_microcins.bestHits.fa"
# 	shell:
# 		"touch {output}"

# rule bestHitsMSA:
# 	input:
# 		"verified_microcins.bestHits.fa"
# 	output:
# 		"verified_microcins.bestHits.aln"
# 	shell:
# 		"touch {output}"

# rule evaluateMSA:
# 	input:
# 		"verified_microcins.bestHits.aln"
# 	output:
# 		"verified_microcins.evaluateMSA.txt"
# rule checkSignalSeq:
# 	input:
# 		"verified_microcins.bestHits.aln"
# 	output:
# 		"verified_microcins.signalCheck.txt"

# rule neuBI:
# 	input:
# 		"verified_microcins.bestHits.fa"
# 	output:
# 		"verified_microcins.bestHits.neuBI.txt"
# 	shell:
# 		"neuBI {input} {output}"

# rule putative_microcins:
# 	input:
# 		seqs = "verified_microcins.bestHits.fa",
# 		evaluateMSA = "verified_microcins.evaluateMSA.txt",
# 		signalCheck = "verified_microcins.signalCheck.txt",
# 		neuBI = "verified_microcins.bestHits.neuBI.txt"
# 	output:
# 		"putative_microcins.txt"





		
