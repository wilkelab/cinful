from io import StringIO
from Bio import SeqIO




rule makeblastdb_immunity_protein:
	input:
		"cinfulOut/00_dbs/immunity_proteins.verified.pep"
	output:
		"cinfulOut/00_dbs/immunity_proteins.verified.pep.phr"
	benchmark:
		"cinfulOut/benchmark/makeblastdb_immunity_protein.txt"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast_immunity_protein:
	input:
		verified_component = "cinfulOut/00_dbs/immunity_proteins.verified.pep",
		blastdb = "cinfulOut/00_dbs/immunity_proteins.verified.pep.phr",
		input_seqs = "cinfulOut/01_orf_homology/immunity_proteins/filtered_nr.fa"
	output:
		"cinfulOut/01_orf_homology/immunity_proteins/blast.txt"
	threads:threads_max
	benchmark:
		"cinfulOut/benchmark/blast_immunity_protein.txt"
	shell:
		"blastp -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1 -num_threads {threads}"

rule msa_immunity_protein:
	input:
		"cinfulOut/00_dbs/immunity_proteins.verified.pep"
	output:
		"cinfulOut/00_dbs/immunity_proteins.verified.aln"
	benchmark:
		"cinfulOut/benchmark/msa_immunity_protein.txt"
	shell:
		"mafft --auto {input} > {output}"

rule buildhmm_immunity_protein:
	input:
		"cinfulOut/00_dbs/immunity_proteins.verified.aln"
	output:
		"cinfulOut/00_dbs/immunity_proteins.verified.hmm"
	benchmark:
		"cinfulOut/benchmark/buildhmm_immunity_protein.txt"
	shell:
		"hmmbuild {output} {input}"



rule blast_v_hmmer_immunity_protein:
	input:
		verifiedHMM = "cinfulOut/00_dbs/immunity_proteins.verified.hmm",
		input_seqs = "cinfulOut/01_orf_homology/immunity_proteins/filtered_nr.fa",
		blastOut = "cinfulOut/01_orf_homology/immunity_proteins/blast.txt"
	output:
		"cinfulOut/01_orf_homology/immunity_proteins/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)#hmmer_hitsHeaders in blastDF["qseqid"]
		blastDF.to_csv(output[0], index = False)
