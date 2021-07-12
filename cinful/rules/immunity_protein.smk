from io import StringIO
from Bio import SeqIO

verified_immunity_proteins = """
>classIIb_IP:M_mcmI_tr|Q83TS3|Q83TS3_ECOLX Microcin M imunity protein McmI OS=Escherichia coli OX=562 GN=mcmI PE=4 SV=1
MGEVKKDIKITVIAFVINYLFFYIPVSLYLSYYYGYNFFNLYMFFLSLVVTFLSLWLNVN
FYFFTNLIAKVLK
>classIIb_IP:I47_mchS3_tr|Q712P9|Q712P9_ECOLX MchS3 protein OS=Escherichia coli OX=562 GN=mchS3 PE=4 SV=1
MYLTKKIIISMMFILPSAAFSSDPPPLQQSLEKTTYFSIGMNGFIGYQSEGEKLYTHILT
LDNPEEIFKNIIKNRKSTKESKIYAACGLYYLNVENIESLFNENDKQEYVSVLRGDILTK
IKLNDILNSVIINGCNTKLISEHK
>classIIb_IP:H47_mchI_sp|O86200|MCHI_ECOLX Microcin H47 immunity protein MchI OS=Escherichia coli OX=562 GN=mchI PE=1 SV=1
MSYKKLYQLTAIFSLPLTILLVSLSSLRIVGEGNSYVDVFLSFIIFLGFIELIHGIRKIL
VWSGWKNGS
>classIIb_IP:E492_mceB_sp|Q9ZHG0|IM92_KLEPN Microcin E492 immunity protein OS=Klebsiella pneumoniae OX=573 GN=mceB PE=1 SV=1
MTLLSFGFSPVFFSVMAFCIISRSKFYPQRTRNKVIVLILLTFFICFLYPLTKVYLVGSY
GIFDKFYLFCFISTLIAIAINVVILTINGAKNERN
>classIIb_IP:G492_mceM_tr|B4DCT4|B4DCT4_KLEPN MceM OS=Klebsiella pneumoniae OX=573 GN=SAMEA4364603_01604 PE=4 SV=1
MIFLYLDKIPLFILGIGLLTSFALPGSSALDSPKFLCIYSSTILAGISFIYQVFRHGTNT
EFFLAMLITVSFVVMLPVIKMHFAY
>classIIa_IP:V_cvi
MDRKRTKLELLFAFIINATAIYIALAIYDCVFRGKDFLSMHTFCFSALMSAICYFVGDNYYSISDKIKRRSYENSDSK
>classIIa_IP:L_mclI
MKTWQVFFIILPISIIISLIVKQLNSSNLVQSVVSGIAIALMISIFFNRGK
>classIIa_IP:N_mcnI
MKRNKLTRMSFLNFAFSPVFFSIMACYFIVWRNKRNEFVCNRLLSIIIISFLICFIYPWLNYKIEVKYYIFEQFYLFCFLSSLVAVVINLIVYFILYRRCI
>classIIa_IP:PDI_mcpI
MEGATMFIKLLSFICGLLLGFALLSGSSVIDLYWFSLPSEFSKIVVMLITLFSTARFMDYIIEKIRTISAK
>classIIa_IP:S_mcsI
MDERSSQFRYSKYSAIIFLAVVIISTIVTLSPTFTLRYVGLDIAFFIVFITEILISTLVYLFYLKEFPECRIKIRTDSATVKFSALSFLIIILIQLAVYCYRDYLYHYEPSQINWITVLVMTLVVPYYEEIVYRACAFGFLRSIFKENIIIPCVITSLFFSLMHFQYYNVLDQSVLFVVSMLLLGVRIKSRSLFYPMLIHSGMNTFVILLNIQNIL
"""

verified_immunity_proteins_SeqIO = SeqIO.parse(StringIO(verified_immunity_proteins), "fasta")
with open("immunity_proteins.verified.pep","w") as seq_out:
  SeqIO.write(verified_immunity_proteins_SeqIO, seq_out,"fasta")


SAMPLES, = glob_wildcards("cinfulOut/01_orf_homology/{sample}_prodigal/")

# rule final:
	# input:
		# expand("cinfulOut/01_orf_homology/{sample}_prodigal/immunity_proteins/{sample}.filtered.fa", sample = SAMPLES)

rule filter_immunity_protein:
	input:
		"cinfulOut/01_orf_homology/{sample}_prodigal/{sample}.faa"
	output:
		"cinfulOut/01_orf_homology/{sample}_prodigal/immunity_proteins/{sample}.filtered.fa"
	shell:
		"seqkit seq -m 30 -M 250  {input} | seqkit rmdup -s > {output}"


# rule makeblastdb:
# 	input:
# 		"verified_immunity_proteins.pep"
# 	output:
# 		"verified_immunity_proteins.pep.phr"
# 	shell:
# 		"makeblastdb -dbtype prot -in {input}"

# rule blast:
# 	input:
# 		verified_immunity_proteins = "verified_immunity_proteins.pep",
# 		blastdb = "verified_immunity_proteins.pep.phr",
# 		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa"
# 	output:
# 		"{sample}_cinfulOut/{sample}.verified_immunity_proteins.blast.txt"
# 	shell:
# 		"blastp -db {input.verified_immunity_proteins} -query {input.input_seqs} -outfmt 6 -out {output} -evalue 0.001 -max_target_seqs 1"

# rule verified_immunity_proteinsMSA:
# 	input:
# 		"verified_immunity_proteins.pep"
# 	output:
# 		"verified_immunity_proteins.aln"
# 	shell:
# 		"mafft --auto {input} > {output}"

# rule buildhmm:
# 	input:
# 		"verified_immunity_proteins.aln"
# 	output:
# 		"verified_immunity_proteins.hmm"
# 	shell:
# 		"hmmbuild {output} {input}"

# rule duomolog:
# 	input:
# 		verified_immunity_proteins = "verified_immunity_proteins.pep",
# 		input_seqs = "{sample}_cinfulOut/{sample}.30_150.fa",
# 		blastout="{sample}_cinfulOut/{sample}.verified_immunity_proteins.blast.txt",
# 		hmm="verified_immunity_proteins.hmm"
# 	output:
# 		"{sample}_cinfulOut/duomolog_immunity_protein/summary_out.txt"
# 	shell:
# 		"""duomolog blast_v_hmmer --inFile {input.verified_immunity_proteins} --queryFile {input.input_seqs} \
# 			--blastFile {input.blastout} \
# 			--intersectOnly \
# 			--hmmFile {input.hmm}	\
# 			--summaryOut {output}
# 		"""		





# rule getBestHits:
# 	input:
# 		blast_hits = "verified_immunity_proteins.blast.txt",
# 		hmmer_hits = "verified_immunity_proteins.hmmerOut.txt",
# 		seq = "input_seqs.short.fa"
# 	output:
# 		"verified_immunity_proteins.bestHits.fa"
# 	shell:
# 		"touch {output}"

# rule bestHitsMSA:
# 	input:
# 		"verified_immunity_proteins.bestHits.fa"
# 	output:
# 		"verified_immunity_proteins.bestHits.aln"
# 	shell:
# 		"touch {output}"

# rule evaluateMSA:
# 	input:
# 		"verified_immunity_proteins.bestHits.aln"
# 	output:
# 		"verified_immunity_proteins.evaluateMSA.txt"

# rule subcellular_localization:
# 	input:
# 		"verified_immunity_proteins.bestHits.fa"
# 	output:
# 		"verified_immunity_proteins.bestHits.subcellular_localization.txt"
# 	shell:
# 		"touch {output}"

# rule transmembrane_helix:
# 	input:
# 		"verified_immunity_proteins.bestHits.fa"
# 	output:
# 		"verified_immunity_proteins.bestHits.transmembrane_helix.txt"
# 	shell:
# 		"touch {output}"



# rule putative_immunity_proteins:
# 	input:
# 		seqs = "verified_immunity_proteins.bestHits.fa",
# 		evaluateMSA = "verified_immunity_proteins.evaluateMSA.txt",
# 		subcellular_localization = "verified_immunity_proteins.bestHits.subcellular_localization.txt",
# 		transmembrane_helix = "verified_immunity_proteins.bestHits.transmembrane_helix.txt"

# 	output:
# 		"putative_immunity_proteins.txt"