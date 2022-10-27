import pandas as pd
from Bio import SeqIO
import pyTMHMM

SAMPLES, = glob_wildcards("{sample}.fna")

def prodigalFa2DF(fa):
    prodigalDict = {"id":[],"seq":[],"contig":[],"start":[],"stop":[],"strand":[]}
    for seq_record in SeqIO.parse(fa,"fasta"):
        descriptionParts = seq_record.description.split("#")
    start = descriptionParts[1].strip()
    stop = descriptionParts[2].strip()
    strand = descriptionParts[3].strip()
    contig = '_'.join(seq_record.id.split("_")[:-1])
    prodigalDict["id"].append(seq_record.id)
    prodigalDict["seq"].append(seq_record.seq)
    prodigalDict["contig"].append(contig)
    prodigalDict["stop"].append(stop)
    prodigalDict["start"].append(start)
    prodigalDict["strand"].append(strand)
    return pd.DataFrame.from_dict(prodigalDict)

def homology_withProdigal(homology_results_file, prodigalDF):
	homology_resultsDF = pd.read_csv(homology_results_file)
	return homology_resultsDF.merge(prodigalDF, left_on="qseqid",right_on="pephash")

def componentDFs(homology_withProdigalDF):
    microcinDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="microcins.verified"]
    immunity_proteinDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="immunity_proteins.verified"]
    CvaBDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="CvaB.verified"]
    MFPDF = homology_withProdigalDF[homology_withProdigalDF["component"]=="MFP.verified"]
    return microcinDF, immunity_proteinDF, CvaBDF, MFPDF

def bestHits(microcinDF, immunity_proteinDF, CvaBDF, MFPDF):
    best_microcin_idx = microcinDF.groupby(['contig'])['bitscore'].transform(max) == microcinDF['bitscore']
    best_microcinDF = microcinDF[best_microcin_idx]
    best_immunity_protein_idx = immunity_proteinDF.groupby(['contig'])['bitscore'].transform(max) == immunity_proteinDF['bitscore']
    best_immunity_proteinDF = immunity_proteinDF[best_immunity_protein_idx]
    best_CvaB_idx = CvaBDF.groupby(['contig'])['bitscore'].transform(max) == CvaBDF['bitscore']
    best_CvaBDF = CvaBDF[best_CvaB_idx]
    best_MFP_idx = MFPDF.groupby(['contig'])['bitscore'].transform(max) == MFPDF['bitscore']
    best_MFPDF = MFPDF[best_MFP_idx]
    return pd.concat([best_microcinDF, best_immunity_proteinDF, best_CvaBDF, best_MFPDF])

def contigs_wMicrocins(best_hitsFile):
    best_hitsDF = pd.read_csv(best_hitsFile)
    needed_columns = ["cinful_id","contig","start","stop","strand"]
    renameDict_microcin = {"cinful_id":"microcin_id","start":"microcin_start","stop":"microcin_stop","strand":"microcin_strand"}
    microcin_contigs = best_hitsDF[best_hitsDF["component"] == "microcins.verified"][needed_columns].rename(columns=renameDict_microcin)
    microcin_contig_componentDF = best_hitsDF[best_hitsDF["contig"].isin(microcin_contigs["contig"])]

    renameDict_CvaB = {"cinful_id":"CvaB_id","start":"CvaB_start","stop":"CvaB_stop","strand":"CvaB_strand"}
    microcin_contig_CvaB = microcin_contig_componentDF[microcin_contig_componentDF["component"] == "CvaB.verified"][needed_columns].rename(columns=renameDict_CvaB)

    renameDict_MFP = {"cinful_id":"MFP_id","start":"MFP_start","stop":"MFP_stop","strand":"MFP_strand"}
    microcin_contig_MFP = microcin_contig_componentDF[microcin_contig_componentDF["component"] == "MFP.verified"][needed_columns].rename(columns=renameDict_MFP)

    return microcin_contigs.merge(microcin_contig_CvaB, left_on="contig",right_on="contig").merge(microcin_contig_MFP, left_on="contig",right_on="contig")

def nearestImmunityProtein(immunityDB, bestMicrocinHits):
	candidateList = []
	for row in bestMicrocinHits.to_dict(orient="records"):
		immunitySubset = immunityDB[immunityDB["sample"] + immunityDB["contig"] == row["sample"]+row["contig"] ].sort_values("start") #Finds all prodigal proteins with the same sample and contig as the microcin hits
		immunitySubset["microcinHit"] = row["cinful_id"] #creates a column for Microcin Hits
		immunityUpstream = immunitySubset[immunitySubset["start"] < row["start"]].tail(3) #searches up stream from microcin for 3 CDS
		immunityDownstream = immunitySubset[immunitySubset["start"] > row["stop"]].head(3) #searches downstream from microcin for 3 CDS
		nearestCandidates = pd.concat([immunityUpstream,immunityDownstream]) #concats up and downstream
		candidateList.append(nearestCandidates) #concats nearestCandidates to the empty candidateList created above
	return pd.concat(candidateList)

def tmhmmCol(df,seqCol="seq"):
    tmhmmAnnotations = []
    for seq in df["seq"]:
        tmhmmAnnotation = pyTMHMM.predict(seq.strip("*"), compute_posterior=False)
        tmhmmAnnotations.append(tmhmmAnnotation)
    return tmhmmAnnotations

def cvab_mfp_neighbor(CvaBDF, prodigalDF, mfp_hmmFile):
	seqLen = prodigalDF["seq"].str.len()
	prodigalDF_mfpLen = prodigalDF[(seqLen <= 450) & (seqLen >= 375)]
	hmm = load_hmm(mfp_hmmFile)
	candidateList = []
	for row in CvaBDF.to_dict(orient="records"):
		prodigalDF_mfpLen_cvab_contig = prodigalDF_mfpLen[prodigalDF_mfpLen["sample"] + prodigalDF_mfpLen["contig"]== row["sample"]+row["contig"] ].sort_values("start")
		prodigalDF_mfpLen_cvab_contig["cvab_hit"] = row["cinful_id"]
		mfpUpstream = prodigalDF_mfpLen_cvab_contig[prodigalDF_mfpLen_cvab_contig["start"] < row["start"]].tail(5)
		mfpDownstream = prodigalDF_mfpLen_cvab_contig[prodigalDF_mfpLen_cvab_contig["start"] > row["stop"]].head(5)
		nearestCandidates = pd.concat([mfpUpstream,mfpDownstream])
		candidateList.append(nearestCandidates)
	mfp_candidates =  pd.concat(candidateList)
	sequences = []
	for row in mfp_candidates.to_dict(orient="records"):
		text_seq = pyhmmer.easel.TextSequence(name =bytes(row["pephash"],"utf-8"), sequence = row["seq"])
		sequences.append(text_seq.digitize(hmm.alphabet))
	pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
	hits = pipeline.search_hmm(hmm, sequences)
	return mfp_candidates[mfp_candidates["pephash"].isin([hit.name.decode() for hit in hits])]

rule best_hits:
	input:
		merged_homology_results = config["outdir"] + "/02_homology_results/all_merged.csv",
		signalSeq = config["outdir"] + "/01_orf_homology/microcins/signalSeq.hit.csv"
	output:
		config["outdir"] + "/03_best_hits/best_hits.csv"
	run:
		homologyDF = pd.read_csv(input.merged_homology_results)
		signalSeqDF = pd.read_csv(input.signalSeq)
		microcinDF, immunity_proteinDF, CvaBDF, MFPDF = componentDFs(homologyDF)
		if not microcinDF.empty:
			best_hitsDF = bestHits(microcinDF, immunity_proteinDF, CvaBDF, MFPDF)
			best_hitsDF["signalMatch"] = best_hitsDF["pephash"].isin(signalSeqDF["signalMatch"])
			best_hitsDF.to_csv(output[0], index = False)

rule bestHitsContigs:
	input:
		config["outdir"] + "/03_best_hits/best_hits.csv"
	output:
		config["outdir"] + "/03_best_hits/best_hit_contigs.csv"
	threads:threads_max
	run:
		microcinContigsDF = contigs_wMicrocins(input[0])
		microcinContigsDF.to_csv(output[0],index = False)

rule candidate_immunity:
	input:
		bestHits = config["outdir"] + "/03_best_hits/best_hits.csv",
		prodigalNR = config["outdir"] + "/01_orf_homology/prodigal_out.all.nr_expanded.csv",
		immunityHomologs = config["outdir"] + "/01_orf_homology/immunity_proteins/blast_v_hmmer.csv"
	output:
		config["outdir"] + "/03_best_hits/best_immunity_protein_candidates.csv"
	threads:threads_max
	run:

		immunityDB = pd.read_csv(input.prodigalNR) #inputs all CDS from prodigal non-redundant
		immunityDB = immunityDB[immunityDB["allStandardAA"]] #ensures all AA are standard 20

		bestHits = pd.read_csv(input.bestHits) #reads in best_hits.csv file
		bestMicrocinHits = bestHits[bestHits["component"] == "microcins.verified"] #Filters best_hits.csv for only microcins
		nearestImmunityDF = nearestImmunityProtein(immunityDB, bestMicrocinHits) #function finds all prodigal hits with a matching sample/contig to microcin hits and then filters them to +/- 3 CDS surrounding a microcin hit

		seqLen = nearestImmunityDF["seq"].str.len() #creates an object 'seqLen' that stores all sequence lengths
		nearestImmunityDF = nearestImmunityDF[(seqLen <= 250 ) & (seqLen >= 30 )] #limits immunityDB to those sequences that are between 30 and 250 AA

		nearestImmunityDF["tmhmm"] = tmhmmCol(nearestImmunityDF) #Adds tmhmm prediction column

		immunityHomologs = pd.read_csv(input.immunityHomologs) #reads in blast_v_hmmer.csv for immunity proteins
		nearestImmunityDF["homologyHit"] = nearestImmunityDF["pephash"].isin(immunityHomologs["qseqid"]) #adds homology (True/False) hit if the nearest immunity protein is in the immunityHomologs

		nearestImmunityDF.to_csv(output[0], index = None)

rule candidate_MFP:
	input:
		bestHits = config["outdir"] + "/03_best_hits/best_hits.csv",
		mfp_hmm = config["outdir"] + "/00_dbs/MFP.verified.hmm",
		prodigalDB = config["outdir"] + "/01_orf_homology/prodigal_out.all.nr_expanded.csv"
	output:
		config["outdir"] + "/03_best_hits/best_MFP_candidates.csv"
	threads:threads_max
	run:
		bestHits = pd.read_csv(input.bestHits)
		CvaBDF = bestHits[bestHits["component"] == "CvaB.verified"]
		prodigalDF = pd.read_csv(input.prodigalDB)
		best_MFP_candidates = cvab_mfp_neighbor(CvaBDF, prodigalDF, input.mfp_hmm)
		best_MFP_candidates.to_csv(output[0], index = None)
