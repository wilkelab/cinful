from io import StringIO
from Bio import SeqIO

rule makeblastdb_microcin:
	input:
		config["outdir"] + "/00_dbs/microcins.verified.pep"
	output:
		config["outdir"] + "/00_dbs/microcins.verified.pep.phr"
	threads:threads_max
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast_microcin:
	input:
		verified_component = config["outdir"] + "/00_dbs/microcins.verified.pep",
		blastdb = config["outdir"] + "/00_dbs/microcins.verified.pep.phr",
		input_seqs = config["outdir"] + "/01_orf_homology/microcins/filtered_nr.fa"
	output:
		config["outdir"] + "/01_orf_homology/microcins/blast.txt"
	threads:threads_max
	shell:
		"blastp -num_threads {threads} -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

rule msa_microcin:
	input:
		config["outdir"] + "/00_dbs/microcins.verified.pep"
	output:
		config["outdir"] + "/00_dbs/microcins.verified.aln"
	threads:threads_max
	shell:
		"mafft --thread {threads} --auto {input} > {output}"

rule buildhmm_microcin:
	input:
		config["outdir"] + "/00_dbs/microcins.verified.aln"
	output:
		config["outdir"] + "/00_dbs/microcins.verified.hmm"
	threads:threads_max
	shell:
		"hmmbuild --cpu {threads} {output} {input}"

rule signalSeqHMM:
	input:
		input_seqs = config["outdir"] + "/01_orf_homology/microcins/filtered_nr.fa",
		signalSeqAln = config["outdir"] + "/00_dbs/verified_SP.aln"
	output:
		config["outdir"] + "/01_orf_homology/microcins/signalSeq.hit.csv"
	run:
		signalSeqHMM = build_hmm(input.signalSeqAln)

		signalSeqHits = hmmsearch(input.input_seqs, signalSeqHMM)
		signalSeqHitStr = [hit.name.decode('utf-8') for hit in signalSeqHits]
		matchDF = pd.DataFrame.from_dict({"signalMatch":signalSeqHitStr})
		matchDF.to_csv(output[0])

rule blast_v_hmmer_microcin:
	input:
		verifiedHMM = config["outdir"] + "/00_dbs/microcins.verified.hmm",
		input_seqs = config["outdir"] + "/01_orf_homology/microcins/filtered_nr.fa",
		blastOut = config["outdir"] + "/01_orf_homology/microcins/blast.txt"
	output:
		config["outdir"] + "/01_orf_homology/microcins/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)
		blastDF.to_csv(output[0], index = False)
