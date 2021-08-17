from io import StringIO
from Bio import SeqIO



SAMPLES, = glob_wildcards("cinfulOut/01_orf_homology/{sample}_prodigal/")


# rule final:
	# input:
		# expand("cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/{sample}.filtered.fa", sample = SAMPLES)





# rule filter_CvaB:
# 	input:
# 		"cinfulOut/01_orf_homology/{sample}_prodigal/{sample}.faa"
# 	output:
# 		"cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/{sample}.filtered.fa"
# 	shell:
# 		"seqkit rmdup -s {input} > {output}" 



rule makeblastdb_CvaB:
	input:
		"cinfulOut/00_dbs/CvaB.verified.pep"
	output:
		"cinfulOut/00_dbs/CvaB.verified.pep.phr"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast_CvaB:
	input:
		verified_component = "cinfulOut/00_dbs/CvaB.verified.pep",
		blastdb = "cinfulOut/00_dbs/CvaB.verified.pep.phr",
		input_seqs = "cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/{sample}.filtered.fa"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/blast.txt"
	shell:
		"blastp -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

rule msa_CvaB:
	input:
		"cinfulOut/00_dbs/CvaB.verified.pep"
	output:
		"cinfulOut/00_dbs/CvaB.verified.aln"
	shell:
		"mafft --auto {input} > {output}"

rule buildhmm_CvaB:
	input:
		"cinfulOut/00_dbs/CvaB.verified.aln"
	output:
		"cinfulOut/00_dbs/CvaB.verified.hmm"
	shell:
		"hmmbuild {output} {input}"



rule blast_v_hmmer_CvaB:
	input:
		verifiedHMM = "cinfulOut/00_dbs/CvaB.verified.hmm",
		input_seqs = "cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/{sample}.filtered.fa",
		blastOut = "cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/blast.txt"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/CvaB/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)#hmmer_hitsHeaders in blastDF["qseqid"]
		blastDF.to_csv(output[0], index = False)

# NOTE: at some point I will need to check if the CvaB results actually make sense, and aren't just latching on to
# ABC transporter domain
# perhaps just some blast filtering will suffice, based on overlap. 

# rule makeblastdb:
# 	input:
# 		"verified_CvaB.pep"
# 	output:
# 		"verified_CvaB.pep.phr"		
# 	shell:
# 		"makeblastdb -dbtype prot -in {input}"

# rule blast:
# 	input:
# 		verified_CvaB = "verified_CvaB.pep",
# 		blastdb = "verified_CvaB.pep.phr",
# 		input_seqs = "{sample}_cinfulOut/{sample}.faa"
# 	output:
# 		"{sample}_cinfulOut/{sample}.verified_CvaB.blast.txt"
# 	shell:
# 		"blastp -db {input.verified_CvaB} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"
		


# rule verified_CvaBMSA:
# 	input:
# 		"verified_CvaB.pep"
# 	output:
# 		"verified_CvaB.aln"
# 	shell:
# 		"mafft --auto {input} > {output}"
		
# rule buildhmm:
# 	input:
# 		"verified_CvaB.aln"
# 	output:
# 		"verified_CvaB.hmm"
# 	shell:
# 		"hmmbuild {output} {input}"

# rule duomolog:
# 	input:
# 		verified_CvaB = "verified_CvaB.pep",
# 		input_seqs = "{sample}_cinfulOut/{sample}.faa",
# 		blastout = "{sample}_cinfulOut/{sample}.verified_CvaB.blast.txt",
# 		hmm = "verified_CvaB.hmm"
# 	output:
# 		"{sample}_cinfulOut/duomolog_CvaB/summary_out.txt"
# 	shell:
# 		"""duomolog blast_v_hmmer --inFile {input.verified_CvaB} --queryFile {input.input_seqs} \
# 			--blastFile {input.blastout} \
# 			--intersectOnly \
# 			--hmmFile {input.hmm}	\
# 			--summaryOut {output}
# 		"""

# rule getBestHits:
# 	input:
# 		blast_hits = "verified_CvaB.blast.txt",
# 		hmmer_hits = "verified_CvaB.hmmerOut.txt",
# 		seq = "input_seqs.short.fa"
# 	output:
# 		"verified_CvaB.bestHits.fa"
# 	shell:
# 		"touch {output}"

# rule bestHitsMSA:
# 	input:
# 		"verified_CvaB.bestHits.fa"
# 	output:
# 		"verified_CvaB.bestHits.aln"
# 	shell:
# 		"touch {output}"

# rule evaluateMSA:
# 	input:
# 		"verified_CvaB.bestHits.aln"
# 	output:
# 		"verified_CvaB.evaluateMSA.txt"

# rule catalytic_triad:
# 	input:
# 		"verified_CvaB.evaluateMSA.txt"
# 	output:
# 		"verified_CvaB.catalytic_triad.txt"


# rule putative_CvaB:
# 	input:
# 		seqs = "verified_CvaB.bestHits.fa",
# 		evaluateMSA = "verified_CvaB.evaluateMSA.txt",
# 		catalytic_triad = "verified_CvaB.catalytic_triad.txt"

# 	output:
# 		"putative_CvaB.txt"