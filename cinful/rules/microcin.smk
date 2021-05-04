from io import StringIO
from Bio import SeqIO

verified_microcins ="""
>E492_sp_Q9Z4N4_MCEA_KLEPN Microcin E492 OS=Klebsiella pneumoniae OX=573 GN=mceA PE=1 SV=2
MREISQKDLNLAFGAGETDPNTQLLNDLGNNMAWGAALGAPGGLGSAALGAAGGALQTVG
QGLIDHGPVNVPIPVLIGPSWNGSGSGYNSATSSSGSGS
>H47_sp_P62530_MCHB_ECOLX Microcin H47 OS=Escherichia coli OX=562 GN=mchB PE=1 SV=1
MREITESQLRYISGAGGAPATSANAAGAAAIVGALAGIPGGPLGVVVGAVSAGLTTAIGS
TVGSGSASSSAGGGS
>I47_tr_Q712Q0_Q712Q0_ECOLX MchS2 protein OS=Escherichia coli OX=562 GN=mchS2 PE=4 SV=1
MREISDNMLDSVKGGMNLNGLPASTNVIDLRGKDMGTYIDANGACWAPDTPSIIMYPGGS
GPSYSMSSSTSSANSGS
>M_tr_Q83TS1_Q83TS1_ECOLX Colicin-V (Microcin-V bacteriocin) OS=Escherichia coli OX=562 GN=mcmA PE=4 SV=1
MRKLSENEIKQISGGDGNDGQAELIAIGSLAGTFISPGFGSIAGAYIGDKVHSWATTATV
SPSMSPSGIGLSSQFGSGRGTSSASSSAGSGS
>G492_tr_B4DCT5_B4DCT5_KLEPN MceL OS=Klebsiella pneumoniae OX=573 PE=4 SV=1
MRALTENDFFAVSGADRGDAAVAGAVAGGTAGAAAGGWAGAQMGATVGSLAGPVGTVVGF
VAGAAAGRYGGAFIYDSFSSPSNSSSSGS
>V_sp_P22522_CEAV_ECOLX Colicin-V OS=Escherichia coli OX=562 GN=cvaC PE=1 SV=1
MRTLTLNELDSVSGGASGRDIAMAIGTLSGQFVAGGIGAAAGGVAGGAIYDYASTHKPNP
AMSPSGLGGTIKQKPEGIPSEAWNYAAGRLCNWSPNNLSDVCL
>L_tr_Q841V4_Q841V4_ECOLX Microcin L OS=Escherichia coli OX=562 GN=mclC PE=4 SV=1
MREITLNEMNNVSGAGDVNWVDVGKTVATNGAGVIGGAFGAGLCGPVCAGAFAVGSSAAV
AALYDAAGNSNSAKQKPEGLPPEAWNYAEGRMCNWSPNNLSDVCL
>N_tr_C3VUZ5_C3VUZ5_ECOLX McnN OS=Escherichia coli OX=562 PE=4 SV=1
MRELDREELNCVGGAGDPLADPNSQIVRQIMSNAAWGAAFGARGGLGGMAVGAAGGVTQT
VLQGAAAHMPVNVPIPKVPMGPSWNGSKG
>PDI_tr_I6ZU90_I6ZU90_ECOLX McpM OS=Escherichia coli OX=562 GN=mcpM PE=4 SV=1
MANIRELTLDEITLVSGGNANSNFEGGPRNDRSSGARNSLGRNAPTHIYSDPSTVKCANA
VFSGMIGGAIKGGPIGMARGTIGGAVVGQCLSDHGSGNGSGNRGSSSSCSGNNVGGTCNR
>S_tr_H9ZMG7_H9ZMG7_ECOLX Microcin S OS=Escherichia coli OX=562 GN=mcsS PE=4 SV=1
MSNIRELSFDEIALVSGGNANSNYEGGGSRSRNTGARNSLGRNAPTHIYSDPSTVKCANA
VFSGMVGGAIKGGPVGMTRGTIGGAVIGQCLSGGGNGNGGGNRAGSSNCSGSNVGGTCSR
"""
verified_microcins_SeqIO = SeqIO.parse(StringIO(verified_microcins), "fasta")
with open("verified_microcins.pep","w") as seq_out:
  SeqIO.write(verified_microcins_SeqIO, seq_out,"fasta")

SAMPLES, = glob_wildcards("{sample}_cinfulOut/")

print(SAMPLES)
rule final:
	input:
		expand("{sample}_cinfulOut/duomolog_microcin/summary_out.txt", sample = SAMPLES)

rule filter_input:
	input:
		"{sample}_cinfulOut/{sample}.faa"
	output:
		"{sample}_cinfulOut/{sample}.30_150.fa"
	shell:
		"seqkit seq -m 30 -M 150 {input} | seqkit rmdup -s > {output}"

rule makeblastdb:
	input:
		"verified_microcins.pep"
	output:
		"verified_microcins.pep.phr"
	shell:
		"makeblastdb -dbtype prot -in {input}"

rule blast:
	input:
		verified_microcins = "verified_microcins.pep",
		blastdb = "verified_microcins.pep.phr",
		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa"
	output:
		"{sample}_cinfulOut/{sample}.verified_microcins.blast.txt"
	shell:
		"blastp -db {input.verified_microcins} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

rule verified_microcinsMSA:
	input:
		"verified_microcins.pep"
	output:
		"verified_microcins.aln"
	shell:
		"mafft --auto {input} > {output}"


rule buildhmm:
	input:
		"verified_microcins.aln"
	output:
		"verified_microcins.hmm"
	shell:
		"hmmbuild {output} {input}"


rule duomolog:
	input:
		verified_microcins = "verified_microcins.pep",
		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa",
		blastout="{sample}_cinfulOut/{sample}.verified_microcins.blast.txt",
		hmm="verified_microcins.hmm"
	output:
		"{sample}_cinfulOut/duomolog_microcin/summary_out.txt"
	shell:
		"""duomolog blast_v_hmmer --inFile {input.verified_microcins} --queryFile {input.input_seqs} \
			--blastFile {input.blastout} \
			--intersectOnly \
			--hmmFile {input.hmm}	\
			--summaryOut {output}
		"""


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





		
