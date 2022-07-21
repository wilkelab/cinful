rule msa_MFP:
	input:
		config["outdir"] + "/00_dbs/MFP.verified.pep"
	output:
		config["outdir"] + "/00_dbs/MFP.verified.aln"
	threads:threads_max
	shell:
		"mafft --thread {threads} --auto {input} > {output}"

rule buildhmm_MFP:
	input:
		config["outdir"] + "/00_dbs/MFP.verified.aln"
	output:
		config["outdir"] + "/00_dbs/MFP.verified.hmm"
	shell:
		"hmmbuild {output} {input}"
