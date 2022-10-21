from io import StringIO
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd

def qcMFP(best_hitsPep, best_hitsAln):
	alignment = AlignIO.read(open(best_hitsAln), "fasta")
	gappyDict = {"id":[], "nonGap":[], "percentGap":[], "alignment":[]}
	for record in alignment:
		percentGap = record.seq.count("-")/len(record.seq)
		gappyDict["id"].append(record.id)
		gappyDict["nonGap"].append( len(record.seq) - record.seq.count("-"))
		gappyDict["percentGap"].append(percentGap)
		gappyDict["alignment"].append(str(record.seq))
	gappyDF = pd.DataFrame.from_dict(gappyDict)
	lenDict = {"id":[], "len":[]}
	for record in SeqIO.parse(best_hitsPep, "fasta"):
		lenDict["id"].append(record.id)
		lenDict["len"].append(len(record.seq))
	lenDF = pd.DataFrame.from_dict(lenDict)
	lenGapDF = gappyDF.merge(lenDF)
	lenGapDF["percentTrim"] = (lenGapDF["len"]-lenGapDF["nonGap"])/ lenGapDF["len"]
	return lenGapDF

rule makeblastdb_MFP:
	input:
		config["outdir"] + "/00_dbs/MFP.verified.pep"
	output:
		config["outdir"] + "/00_dbs/MFP.verified.pep.dmnd"
	threads:threads_max
	shell:
		"diamond makedb --in {input} -d {input} -p {threads}"


rule blast_MFP:
	input:
		verified_component = config["outdir"] + "/00_dbs/MFP.verified.pep.dmnd",
		input_seqs = config["outdir"] + "/01_orf_homology/MFP/filtered_nr.fa"
	output:
		config["outdir"] + "/01_orf_homology/MFP/blast.txt"
	threads:threads_max
	shell:
		"diamond blastp -d {input.verified_component} -q {input.input_seqs}   --evalue 0.001 -k 1 -o {output} -p {threads}"

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
	threads:threads_max
	shell:
		"hmmbuild --cpu {threads} {output} {input}"

rule blast_v_hmmer_MFP:
	input:
		verifiedHMM = config["outdir"] + "/00_dbs/MFP.verified.hmm",
		input_seqs = config["outdir"] + "/01_orf_homology/MFP/filtered_nr.fa",
		blastOut = config["outdir"] + "/01_orf_homology/MFP/blast.txt"
	output:
		config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)
		blastDF.to_csv(output[0], index = False)

rule best_MFP_headers:
	input:
		config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.csv"
	output:
		config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.headers"
	shell:
		"cut -d, -f1 {input} > {output}"

rule best_MFP_fasta:
	input:
		headers=config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.headers",
		input_seqs = config["outdir"] + "/01_orf_homology/MFP/filtered_nr.fa"
	output:
		config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.fa"
	threads:threads_max
	shell:
		"seqkit -j {threads} grep -f {input.headers} {input.input_seqs} > {output}"

rule align_with_verifiedMFP:
	input:
		best_MFP = config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.fa",
		verified_MFP = config["outdir"] + "/00_dbs/MFP.verified.aln"
	output:
		config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.with_verified.aln"
	threads:threads_max
	shell:
		"mafft --thread {threads} --inputorder --keeplength --add {input.best_MFP} --auto "
		"{input.verified_MFP} > {output}"

rule filter_MFP_hits:
	input:
		best_MFP=config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.fa",
		align_with_verifiedMFP=config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.with_verified.aln",
		nr_filtered_csv=config["outdir"] + "/01_orf_homology/prodigal_out.all.nr_expanded.csv"
	output:
		preQC=config["outdir"] + "/01_orf_homology/MFP/preQC.csv",
		QC=config["outdir"] + "/01_orf_homology/MFP/QC.csv"
	run:
		preQC_MFP = qcMFP(input.best_MFP, input.align_with_verifiedMFP)
		nr_filteredDF = pd.read_csv(input.nr_filtered_csv)
		nr_filteredDF.merge(preQC_MFP, left_on="pephash", right_on="id").to_csv(output.preQC)

		lowGapMFP = preQC_MFP[preQC_MFP["percentGap"] <0.1]
		lowTrimMFP = lowGapMFP[lowGapMFP["percentTrim"] <0.1]
		lowTrimMFP.to_csv(output.QC)
