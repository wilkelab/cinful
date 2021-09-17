from io import StringIO
from Bio import SeqIO






rule makeblastdb_microcin:
	input:
		"cinfulOut/00_dbs/microcins.verified.pep"
	output:
		"cinfulOut/00_dbs/microcins.verified.pep.phr"
	benchmark:
		"cinfulOut/benchmark/makeblastdb_microcin.txt"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast_microcin:
	input:
		verified_component = "cinfulOut/00_dbs/microcins.verified.pep",
		blastdb = "cinfulOut/00_dbs/microcins.verified.pep.phr",
		input_seqs = "cinfulOut/01_orf_homology/microcins/filtered_nr.fa"
	output:
		"cinfulOut/01_orf_homology/microcins/blast.txt"
	threads:threads_max
	benchmark:
		"cinfulOut/benchmark/blast_microcin.txt"
	shell:
		"blastp -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1 -num_threads {threads}"

rule msa_microcin:
	input:
		"cinfulOut/00_dbs/microcins.verified.pep"
	output:
		"cinfulOut/00_dbs/microcins.verified.aln"
	benchmark:
		"cinfulOut/benchmark/msa_microcin.txt"
	shell:
		"mafft --auto {input} > {output}"

rule buildhmm_microcin:
	input:
		"cinfulOut/00_dbs/microcins.verified.aln"
	output:
		"cinfulOut/00_dbs/microcins.verified.hmm"
	benchmark:
		"cinfulOut/benchmark/buildhmm_microcin.txt"
	shell:
		"hmmbuild {output} {input}"



rule blast_v_hmmer_microcin:
	input:
		verifiedHMM = "cinfulOut/00_dbs/microcins.verified.hmm",
		input_seqs = "cinfulOut/01_orf_homology/microcins/filtered_nr.fa",
		blastOut = "cinfulOut/01_orf_homology/microcins/blast.txt"
	output:
		"cinfulOut/01_orf_homology/microcins/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)#hmmer_hitsHeaders in blastDF["qseqid"]
		blastDF.to_csv(output[0], index = False)

