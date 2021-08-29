from io import StringIO
from Bio import SeqIO
from Bio import AlignIO
import pandas as pd



def qcCvab(best_hitsPep, best_hitsAln):
	alignment = AlignIO.read(open(best_hitsAln), "fasta")
	gappyDict = {"id":[], "nonGap":[], "percentGap":[], "alignment":[]}

	for record in alignment:	
		percentGap = record.seq.count("-")/len(record.seq)
		gappyDict["id"].append(record.id)
		gappyDict["nonGap"].append( len(record.seq) - record.seq.count("-"))
		gappyDict["percentGap"].append(percentGap)
		gappyDict["alignment"].append(str(record.seq))
	gappyDF = pd.DataFrame.from_dict(gappyDict)

	# catalytic triad
	gappyDF["C34"] = gappyDF["alignment"].astype(str).str[34]
	gappyDF["H107"] = gappyDF["alignment"].astype(str).str[107]
	gappyDF["D123"] = gappyDF["alignment"].astype(str).str[123]

	lenDict = {"id":[], "len":[]}
	for record in SeqIO.parse(best_hitsPep, "fasta"):
		lenDict["id"].append(record.id)
		lenDict["len"].append(len(record.seq))
	lenDF = pd.DataFrame.from_dict(lenDict)

	lenGapDF = gappyDF.merge(lenDF)
	lenGapDF["percentTrim"] = (lenGapDF["len"]-lenGapDF["nonGap"])/ lenGapDF["len"]

	return lenGapDF




# rule final:
	# input:
		# expand("cinfulOut/01_orf_homology/CvaB/filtered_nr.fa", sample = SAMPLES)





# rule filter_CvaB:
# 	input:
# 		"cinfulOut/01_orf_homology/{sample}.faa"
# 	output:
# 		"cinfulOut/01_orf_homology/CvaB/filtered_nr.fa"
# 	shell:
# 		"seqkit rmdup -s {input} > {output}" 



rule makeblastdb_CvaB:
	input:
		"cinfulOut/00_dbs/CvaB.verified.pep"
	output:
		"cinfulOut/00_dbs/CvaB.verified.pep.dmnd"
	shell:
		"diamond makedb --in {input} -d {input}"
		# "makeblastdb -dbtype prot -in {input}"

rule blast_CvaB:
	input:
		verified_component = "cinfulOut/00_dbs/CvaB.verified.pep.dmnd",
		input_seqs = "cinfulOut/01_orf_homology/CvaB/filtered_nr.fa"
	output:
		"cinfulOut/01_orf_homology/CvaB/blast.txt"
	threads:threads_max
	shell:
		"diamond blastp -d {input.verified_component} -q {input.input_seqs}   --evalue 0.001 -k 1 -o {output} -p {threads}"
		# "blastp -db {input.verified_component} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"



rule msa_CvaB:
	input:
		"cinfulOut/00_dbs/CvaB.verified.pep"
	output:
		"cinfulOut/00_dbs/CvaB.verified.aln"
	threads:threads_max
	shell:
		"mafft --thread {threads} --auto {input} > {output}"

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
		input_seqs = "cinfulOut/01_orf_homology/CvaB/filtered_nr.fa",
		blastOut = "cinfulOut/01_orf_homology/CvaB/blast.txt"
	output:
		"cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.csv"
	run:
		blastDF = load_blast(input.blastOut)
		hmmer_hits, hmm_name = run_hmmsearch(input.input_seqs, input.verifiedHMM)
		hmmer_hitsHeaders = [hit.name.decode() for hit in hmmer_hits]
		blastDF["component"] = hmm_name
		blastDF["hmmerHit"] = blastDF["qseqid"].isin(hmmer_hitsHeaders)#hmmer_hitsHeaders in blastDF["qseqid"]
		blastDF.to_csv(output[0], index = False)


rule best_Cvab_headers:
	input:
		"cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.csv"
	output:
		"cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.headers"
	shell:
		"cut -d, -f1 {input} > {output}"

rule best_Cvab_fasta:
	input:
		headers="cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.headers",
		input_seqs = "cinfulOut/01_orf_homology/CvaB/filtered_nr.fa"
	output:
		"cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.fa"
	shell:
		"seqkit grep -f {input.headers} {input.input_seqs} > {output}"


# rule best_Cvab_fa:
# 	input:
# 		hitsCvaB="cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.csv",
# 		nrPeps="cinfulOut/01_orf_homology/CvaB/filtered_nr.fa"
# 	output:
# 		"cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.fa"
# 	run:
# 		hitsCvaBDF = pd.read_csv(input.hitsCvaB)
# 		with open(output[0],"w") as out:
# 			for record in SeqIO.parse(input.nrPeps, "fasta"):
# 				if record.id in set(hitsCvaBDF["qseqid"]):
# 					SeqIO.write(record,out,"fasta")

rule align_with_verifiedCvab:
	input:
		best_CvaB="cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.fa",
		verified_CvaB="cinfulOut/00_dbs/CvaB.verified.aln"
	output:
		"cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.with_verified.aln"
	shell:
		"mafft --inputorder --keeplength --add {input.best_CvaB} --auto "
		"{input[1]} > {output}"

rule filter_CvaB_hits:
	input:
		best_CvaB="cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.fa", 
		align_with_verifiedCvab="cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.with_verified.aln",
		nr_filtered_csv="cinfulOut/01_orf_homology/prodigal_out.all.nr_expanded.csv"
	output:
		preQC="cinfulOut/01_orf_homology/CvaB/preQC.csv",
		QC="cinfulOut/01_orf_homology/CvaB/QC.csv"
	run:
		preQC_Cvab = qcCvab(input.best_CvaB, input.align_with_verifiedCvab)
		nr_filteredDF = pd.read_csv(input.nr_filtered_csv)
		nr_filteredDF.merge(preQC_Cvab, left_on="pephash", right_on="id").to_csv(output.preQC)

		lowGapCvaB = preQC_Cvab[preQC_Cvab["percentGap"] <0.1]
		lowTrimCvaB = lowGapCvaB[lowGapCvaB["percentTrim"] <0.1]
		catalyticTriadCvab = lowTrimCvaB[(lowTrimCvaB["C34"] + lowTrimCvaB["H107"] + lowTrimCvaB["D123"]) == "CHD"]
		catalyticTriadCvab.to_csv(output.QC)



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