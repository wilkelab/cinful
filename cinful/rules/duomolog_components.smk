import pandas as pd
import csv


rule merged_results:
	input:
		#Import nonredundant prodigal (all CDS from all genomes)
		nr_csv = config["outdir"] + "/01_orf_homology/prodigal_out.all.nr_expanded.csv",
		#Import pephash files from verified proteins
		microcin_pepHash = config["outdir"] + "/00_dbs/microcins.verified.pephash.csv",
		immunity_protein_pepHash = config["outdir"] + "/00_dbs/immunity_proteins.verified.pephash.csv",
		CvaB_pepHash = config["outdir"] + "/00_dbs/CvaB.verified.pephash.csv",
		MFP_pepHash = config["outdir"] + "/00_dbs/MFP.verified.pephash.csv",
		#Import all BLAST v HMMER
		microcins = config["outdir"] + "/01_orf_homology/microcins/blast_v_hmmer.csv",
		immunity_proteins = config["outdir"] + "/01_orf_homology/immunity_proteins/blast_v_hmmer.csv",
		unfilteredCvaB = config["outdir"] + "/01_orf_homology/CvaB/blast_v_hmmer.csv",
		unfilteredMFP = config["outdir"] + "/01_orf_homology/MFP/blast_v_hmmer.csv",
		#Import QC files for CvaB and MFP
		QC_Cvab = config["outdir"] + "/01_orf_homology/CvaB/QC.csv",
		QC_MFP = config["outdir"] + "/01_orf_homology/MFP/QC.csv"
	output:
		csv = config["outdir"] + "/02_homology_results/all_merged.csv"
	threads:threads_max
	run:
		with open(output.csv, "w") as csv_file:
			nrDF = pd.read_csv(input.nr_csv)

			microcinDF = pd.read_csv(input.microcins)
			immunity_proteinDF = pd.read_csv(input.immunity_proteins)
			QC_CvabDF = pd.read_csv(input.QC_Cvab)
			QC_MFPDF = pd.read_csv(input.QC_MFP)

			unfilteredCvaBDF = pd.read_csv(input.unfilteredCvaB)
			unfilteredMFPDF = pd.read_csv(input.unfilteredMFP)

			microcin_pepHashDF = pd.read_csv(input.microcin_pepHash)
			immunity_protein_pepHashDF = pd.read_csv(input.immunity_protein_pepHash)
			CvaB_pepHashDF = pd.read_csv(input.CvaB_pepHash)
			MFP_pepHashDF = pd.read_csv(input.MFP_pepHash)

			bestCvaB = unfilteredCvaBDF[unfilteredCvaBDF["qseqid"].isin(QC_CvabDF["id"]) ]
			bestMFP = unfilteredMFPDF[unfilteredMFPDF["qseqid"].isin(QC_MFPDF["id"]) ]

			microcinDF["verified"] = microcinDF["qseqid"].isin(microcin_pepHashDF["pephash"])
			immunity_proteinDF["verified"] = immunity_proteinDF["qseqid"].isin(immunity_protein_pepHashDF["pephash"])
			bestCvaB["verified"] = bestCvaB["qseqid"].isin(CvaB_pepHashDF["pephash"])
			bestMFP["verified"] = bestMFP["qseqid"].isin(MFP_pepHashDF["pephash"])

			componentDFs = [microcinDF, immunity_proteinDF, bestCvaB, bestMFP]
			mergedDf = pd.concat(componentDFs)
#			This is what causes the non-redudancy
			mergedDf_nr = mergedDf.merge(nrDF, left_on ="qseqid", right_on= "pephash")
 			mergedDf_nr.to_csv(csv_file, index = None)
