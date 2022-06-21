from io import StringIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO

def fa2hashDF(fa):
	outDict = {"header":[],"pephash":[],"sequeunce":[]}
	for record in SeqIO.parse(fa,"fasta"):
		pephash = seqhash.seqhash(record.seq,dna_type='PROTEIN')
		outDict["header"].append(record.id)
		outDict["pephash"].append(pephash)
		outDict["sequeunce"].append(str(record.seq))
	return pd.DataFrame.from_dict(outDict)

# For now, these will just have sequences hardcoded in them,
# This guarantees that if you have the codebase, you can do the searches.
# What I will also do is use config options to let users define their own input sequences

rule write_CvaB:
#	input:
#		verified_CvaB = workflow.basedir + "/input/CvaB.verified.pep"
	output:
		CvaB_fasta = config["outdir"] + "/00_dbs/CvaB.verified.pep",
		CvaB_pepHash = config["outdir"] + "/00_dbs/CvaB.verified.pephash.csv"
	run:
		verified_CvaB = open(workflow.basedir + "/input/CvaB.verified.pep")
		verified_CvaB_SeqIO = SeqIO.parse(StringIO(verified_CvaB.read()), "fasta")
		verified_CvaB_pepHashDF = fa2hashDF(StringIO(verified_CvaB.read()))
		verified_CvaB_pepHashDF.to_csv(output.CvaB_pepHash)
		with open(output.CvaB_fasta,"w") as seq_out:
			SeqIO.write(verified_CvaB_SeqIO, seq_out,"fasta")

rule write_MFP:
#	input:
#		verified_MFP = workflow.basedir + "/input/MFP.verified.pep"
	output:
		MFP_fasta = config["outdir"] + "/00_dbs/MFP.verified.pep",
		MFP_pepHash = config["outdir"] + "/00_dbs/MFP.verified.pephash.csv"
	run:
		verified_MFP = open(workflow.basedir + "/input/MFP.verified.pep")
		verified_MFP_SeqIO = SeqIO.parse(StringIO(verified_MFP.read()), "fasta")
		verified_MFP_pepHashDF = fa2hashDF(StringIO(verified_MFP.read()))
		verified_MFP_pepHashDF.to_csv(output.MFP_pepHash)
		with open(output.MFP_fasta,"w") as seq_out:
			SeqIO.write(verified_MFP_SeqIO, seq_out,"fasta")

rule write_microcins:
#	input:
#		verified_microcins = workflow.basedir + "/input/microcins.verified.pep"
	output:
		microcin_fasta = config["outdir"] + "/00_dbs/microcins.verified.pep",
		microcin_pepHash = config["outdir"] + "/00_dbs/microcins.verified.pephash.csv"
	run:
		verified_microcins = open(workflow.basedir + "/input/microcins.verified.pep")
		verified_microcins_SeqIO = SeqIO.parse(StringIO(verified_microcins.read()), "fasta")
		verified_microcins_pepHashDF = fa2hashDF(StringIO(verified_microcins.read()))
		verified_microcins_pepHashDF.to_csv(output.microcin_pepHash)
		with open(output.microcin_fasta,"w") as seq_out:
			SeqIO.write(verified_microcins_SeqIO, seq_out,"fasta")

rule write_immunity_proteins:
#	input:
#		verified_immunity_proteins = workflow.basedir + "/input/immunity_proteins.verified.pep"
	output:
		immunity_protein_fasta = config["outdir"] + "/00_dbs/immunity_proteins.verified.pep",
		immunity_protein_pepHash = config["outdir"] + "/00_dbs/immunity_proteins.verified.pephash.csv"
	run:
		verified_immunity_proteins = open(workflow.basedir + "/input/immunity_proteins.verified.pep")
		verified_immunity_proteins_SeqIO = SeqIO.parse(StringIO(verified_immunity_proteins.read()), "fasta")
		verified_immunity_proteins_pepHashDF = fa2hashDF(StringIO(verified_immunity_proteins.read()))
		verified_immunity_proteins_pepHashDF.to_csv(output.immunity_protein_pepHash)
		with open(output.immunity_protein_fasta,"w") as seq_out:
  			SeqIO.write(verified_immunity_proteins_SeqIO, seq_out,"fasta")

rule microcin_signal:
#	input:
#		verified_SP = workflow.basedir + "/input/SP.verified.pep"
	output:
		config["outdir"] + "/00_dbs/verified_SP.aln"
	run:
		verified_SP = open(workflow.basedir + "/input/SP.verified.pep")
		verified_SP_msa = AlignIO.read(StringIO(verified_SP.read()), "fasta")
		with open(output[0],"w") as alignment_out:
  			AlignIO.write(verified_SP_msa,alignment_out,"fasta")
