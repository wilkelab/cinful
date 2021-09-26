
# mkdir .tests/align_with_verifiedCvab
cp -r .tests/unit/align_with_verifiedCvab/data/ .tests/align_with_verifiedCvab
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.with_verified.aln -F -j1 --keep-target-files --allowed-rules align_with_verifiedCvab --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/align_with_verifiedCvab > .tests/align_with_verifiedCvab.benchmark.txt 
# elapsed time: 9.022s
# peak rss: 155.63 MB

# mkdir .tests/bestHitsContigs
cp -r .tests/unit/bestHitsContigs/data/ .tests/bestHitsContigs
memusg -t  python -m snakemake cinfulOut/03_best_hits/best_hit_contigs.csv -F -j1 --keep-target-files --allowed-rules bestHitsContigs --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/bestHitsContigs > .tests/bestHitsContigs.benchmark.txt 
# elapsed time: 5.596s
# peak rss: 209.04 MB

# mkdir .tests/best_Cvab_fasta
cp -r .tests/unit/best_Cvab_fasta/data/ .tests/best_Cvab_fasta
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.fa -F -j1 --keep-target-files --allowed-rules best_Cvab_fasta --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/best_Cvab_fasta > .tests/best_Cvab_fasta.benchmark.txt 
# elapsed time: 3.437s
# peak rss: 112.46 MB

# mkdir .tests/best_Cvab_headers
cp -r .tests/unit/best_Cvab_headers/data/ .tests/best_Cvab_headers
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.headers -F -j1 --keep-target-files --allowed-rules best_Cvab_headers --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/best_Cvab_headers > .tests/best_Cvab_headers.benchmark.txt 
# elapsed time: 2.957s
# peak rss: 101.21 MB

# mkdir .tests/best_hits
cp -r .tests/unit/best_hits/data/ .tests/best_hits
memusg -t  python -m snakemake cinfulOut/03_best_hits/best_hits.csv -F -j1 --keep-target-files --allowed-rules best_hits --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/best_hits > .tests/best_hits.benchmark.txt 
# elapsed time: 42.838s
# peak rss: 3.66 GB

# mkdir .tests/blast_CvaB
cp -r .tests/unit/blast_CvaB/data/ .tests/blast_CvaB
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/blast.txt -F -j1 --keep-target-files --allowed-rules blast_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/blast_CvaB > .tests/blast_CvaB.benchmark.txt 
# elapsed time: 16.044s
# peak rss: 272.02 MB

## This one failed, said it couldn't find the database file
## Works now, just needed to add the other blast files that are too annoying to add to the rule output
# mkdir .tests/blast_immunity_protein
# cp cinfulOut/00_dbs/immunity_proteins.verified.pep* .tests/blast_immunity_protein/cinfulOut/00_dbs/
cp -r .tests/unit/blast_immunity_protein/data/ .tests/blast_immunity_protein
memusg -t  python -m snakemake cinfulOut/01_orf_homology/immunity_proteins/blast.txt -F -j1 --keep-target-files --allowed-rules blast_immunity_protein --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/blast_immunity_protein > .tests/blast_immunity_protein.benchmark.txt 
# elapsed time: 4m:37s
# peak rss: 142.74 MB

# mkdir .tests/blast_microcin
# cp cinfulOut/00_dbs/microcins.verified.pep* .tests/blast_microcin/cinfulOut/00_dbs/
cp -r .tests/unit/blast_microcin/data/ .tests/blast_microcin
memusg -t  python -m snakemake cinfulOut/01_orf_homology/microcins/blast.txt -F -j1 --keep-target-files --allowed-rules blast_microcin --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/blast_microcin > .tests/blast_microcin.benchmark.txt 
# elapsed time: 2m:23s
# peak rss: 143.12 MB

# mkdir .tests/blast_v_hmmer_CvaB
cp -r .tests/unit/blast_v_hmmer_CvaB/data/ .tests/blast_v_hmmer_CvaB
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/blast_v_hmmer_CvaB > .tests/blast_v_hmmer_CvaB.benchmark.txt 
# elapsed time: 18.601s
# peak rss: 276.37 MB

# mkdir .tests/blast_v_hmmer_immunity_protein
cp -r .tests/unit/blast_v_hmmer_immunity_protein/data/ .tests/blast_v_hmmer_immunity_protein
memusg -t  python -m snakemake cinfulOut/01_orf_homology/immunity_proteins/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_immunity_protein --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/blast_v_hmmer_immunity_protein > .tests/blast_v_hmmer_immunity_protein.benchmark.txt 
# elapsed time: 8.351s
# peak rss: 497.12 MB

# mkdir .tests/blast_v_hmmer_microcin
cp -r .tests/unit/blast_v_hmmer_microcin/data/ .tests/blast_v_hmmer_microcin
memusg -t  python -m snakemake cinfulOut/01_orf_homology/microcins/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_microcin --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/blast_v_hmmer_microcin > .tests/blast_v_hmmer_microcin.benchmark.txt 
# elapsed time: 6.772s
# peak rss: 376.83 MB

# mkdir .tests/buildhmm_CvaB
cp -r .tests/unit/buildhmm_CvaB/data/ .tests/buildhmm_CvaB
memusg -t  python -m snakemake cinfulOut/00_dbs/CvaB.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/buildhmm_CvaB > .tests/buildhmm_CvaB.benchmark.txt 
# elapsed time: 4.976s
# peak rss: 107.07 MB

# mkdir .tests/buildhmm_immunity_protein
cp -r .tests/unit/buildhmm_immunity_protein/data/ .tests/buildhmm_immunity_protein
memusg -t  python -m snakemake cinfulOut/00_dbs/immunity_proteins.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_immunity_protein --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/buildhmm_immunity_protein > .tests/buildhmm_immunity_protein.benchmark.txt 
# elapsed time: 2.980s
# peak rss: 100.6 MB

# mkdir .tests/buildhmm_microcin
cp -r .tests/unit/buildhmm_microcin/data/ .tests/buildhmm_microcin
memusg -t  python -m snakemake cinfulOut/00_dbs/microcins.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_microcin --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/buildhmm_microcin > .tests/buildhmm_microcin.benchmark.txt 
# elapsed time: 3.478s
# peak rss: 105.77 MB

# mkdir .tests/candidate_immunity
cp -r .tests/unit/candidate_immunity/data/ .tests/candidate_immunity
memusg -t  python -m snakemake cinfulOut/03_best_hits/best_immunity_protein_candidates.csv -F -j1 --keep-target-files --allowed-rules candidate_immunity --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/candidate_immunity > .tests/candidate_immunity.benchmark.txt 
# elapsed time: 6m:02s
# peak rss: 3.66 GB

# mkdir .tests/filter_CvaB
cp -r .tests/unit/filter_CvaB/data/ .tests/filter_CvaB
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/filter_CvaB > .tests/filter_CvaB.benchmark.txt 
# elapsed time: 4.319s
# peak rss: 143.75 MB

# mkdir .tests/filter_CvaB_hits
cp -r .tests/unit/filter_CvaB_hits/data/ .tests/filter_CvaB_hits
memusg -t  python -m snakemake cinfulOut/01_orf_homology/CvaB/preQC.csv cinfulOut/01_orf_homology/CvaB/QC.csv -F -j1 --keep-target-files --allowed-rules filter_CvaB_hits --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/filter_CvaB_hits > .tests/filter_CvaB_hits.benchmark.txt 
# elapsed time: 46.550s
# peak rss: 3.66 GB

# mkdir .tests/filter_immunity_protein
cp -r .tests/unit/filter_immunity_protein/data/ .tests/filter_immunity_protein
memusg -t  python -m snakemake cinfulOut/01_orf_homology/immunity_proteins/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_immunity_protein --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/filter_immunity_protein > .tests/filter_immunity_protein.benchmark.txt 
# elapsed time: 7.529s
# peak rss: 172.1 MB

# mkdir .tests/filter_microcin
cp -r .tests/unit/filter_microcin/data/ .tests/filter_microcin
memusg -t  python -m snakemake cinfulOut/01_orf_homology/microcins/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_microcin --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/filter_microcin > .tests/filter_microcin.benchmark.txt 
# elapsed time: 4.999s
# peak rss: 155.56 MB


# mkdir .tests/final
# cp -r .tests/unit/final/data/ .tests/final
# memusg -t  # python -m snakemake final -F -j1 --keep-target-files --allowed-rules final --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/final > .tests/final.benchmark.txt 

# mkdir .tests/makeblastdb_CvaB
cp -r .tests/unit/makeblastdb_CvaB/data/ .tests/makeblastdb_CvaB
memusg -t  python -m snakemake cinfulOut/00_dbs/CvaB.verified.pep.dmnd -F -j1 --keep-target-files --allowed-rules makeblastdb_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/makeblastdb_CvaB > .tests/makeblastdb_CvaB.benchmark.txt 
# elapsed time: 3.002s
# peak rss: 101.07 MB

# mkdir .tests/makeblastdb_immunity_protein
cp -r .tests/unit/makeblastdb_immunity_protein/data/ .tests/makeblastdb_immunity_protein
memusg -t  python -m snakemake cinfulOut/00_dbs/immunity_proteins.verified.pep.phr -F -j1 --keep-target-files --allowed-rules makeblastdb_immunity_protein --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/makeblastdb_immunity_protein > .tests/makeblastdb_immunity_protein.benchmark.txt 
# elapsed time: 3.634s
# peak rss: 106.1 MB

# mkdir .tests/makeblastdb_microcin
cp -r .tests/unit/makeblastdb_microcin/data/ .tests/makeblastdb_microcin
memusg -t  python -m snakemake cinfulOut/00_dbs/microcins.verified.pep.phr -F -j1 --keep-target-files --allowed-rules makeblastdb_microcin --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/makeblastdb_microcin > .tests/makeblastdb_microcin.benchmark.txt 
# elapsed time: 2.832s
# peak rss: 120.99 MB

# mkdir .tests/merged_results
cp -r .tests/unit/merged_results/data/ .tests/merged_results
memusg -t  python -m snakemake cinfulOut/02_homology_results/all_merged.csv -F -j1 --keep-target-files --allowed-rules merged_results --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/merged_results > .tests/merged_results.benchmark.txt 
# elapsed time: 6.232s
# peak rss: 204.87 MB

# mkdir .tests/msa_CvaB
cp -r .tests/unit/msa_CvaB/data/ .tests/msa_CvaB
memusg -t  python -m snakemake cinfulOut/00_dbs/CvaB.verified.aln -F -j1 --keep-target-files --allowed-rules msa_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/msa_CvaB > .tests/msa_CvaB.benchmark.txt 
# elapsed time: 3.567s
# peak rss: 120.98 MB

# mkdir .tests/msa_immunity_protein
cp -r .tests/unit/msa_immunity_protein/data/ .tests/msa_immunity_protein
memusg -t  python -m snakemake cinfulOut/00_dbs/immunity_proteins.verified.aln -F -j1 --keep-target-files --allowed-rules msa_immunity_protein --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/msa_immunity_protein > .tests/msa_immunity_protein.benchmark.txt 
# elapsed time: 3.184s
# peak rss: 121.91 MB

# mkdir .tests/msa_microcin
cp -r .tests/unit/msa_microcin/data/ .tests/msa_microcin
memusg -t  python -m snakemake cinfulOut/00_dbs/microcins.verified.aln -F -j1 --keep-target-files --allowed-rules msa_microcin --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/msa_microcin > .tests/msa_microcin.benchmark.txt 
# elapsed time: 3.232s
# peak rss: 107.75 MB

# So this one gave an error :( I think the most recent version may not give these issues. Though that would mean doing this whole song and dance all over...
# mkdir .tests/nonredundant_prodigal
cp -r .tests/unit/nonredundant_prodigal/data/ .tests/nonredundant_prodigal
memusg -t  python -m snakemake cinfulOut/01_orf_homology/prodigal_out.all.nr.faa cinfulOut/01_orf_homology/prodigal_out.all.nr_expanded.csv -F -j1 --keep-target-files --allowed-rules nonredundant_prodigal --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/nonredundant_prodigal > .tests/nonredundant_prodigal.benchmark.txt 

# mkdir .tests/prodigal
cp -r .tests/unit/prodigal/data/ .tests/prodigal
memusg -t  python -m snakemake cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.gff3 cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.cds cinfulOut/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.faa -F -j1 --keep-target-files --allowed-rules prodigal --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/prodigal > .tests/prodigal.benchmark.txt 

# mkdir .tests/write_CvaB
cp -r .tests/unit/write_CvaB/data/ .tests/write_CvaB
memusg -t  python -m snakemake cinfulOut/00_dbs/CvaB.verified.pep -F -j1 --keep-target-files --allowed-rules write_CvaB --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/write_CvaB > .tests/write_CvaB.benchmark.txt 

# mkdir .tests/write_immunity_proteins
cp -r .tests/unit/write_immunity_proteins/data/ .tests/write_immunity_proteins
memusg -t  python -m snakemake cinfulOut/00_dbs/immunity_proteins.verified.pep -F -j1 --keep-target-files --allowed-rules write_immunity_proteins --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/write_immunity_proteins > .tests/write_immunity_proteins.benchmark.txt 

# mkdir .tests/write_microcins
cp -r .tests/unit/write_microcins/data/ .tests/write_microcins
memusg -t  python -m snakemake cinfulOut/00_dbs/microcins.verified.pep -F -j1 --keep-target-files --allowed-rules write_microcins --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile -d .tests/write_microcins > .tests/write_microcins.benchmark.txt 


# align_with_verifiedCvab	cinful_out/01_orf_homology/CvaB/blast_v_hmmer.with_verified.aln
# bestHitsContigs	cinful_out/03_best_hits/best_hit_contigs.csv
# best_Cvab_fasta	cinful_out/01_orf_homology/CvaB/blast_v_hmmer.fa
# best_Cvab_headers	cinful_out/01_orf_homology/CvaB/blast_v_hmmer.headers
# best_hits	cinful_out/03_best_hits/best_hits.csv
# blast_CvaB	cinful_out/01_orf_homology/CvaB/blast.txt
# blast_immunity_protein	cinful_out/01_orf_homology/immunity_proteins/blast.txt
# blast_microcin	cinful_out/01_orf_homology/microcins/blast.txt
# blast_v_hmmer_CvaB	cinful_out/01_orf_homology/CvaB/blast_v_hmmer.csv
# blast_v_hmmer_immunity_protein	cinful_out/01_orf_homology/immunity_proteins/blast_v_hmmer.csv
# blast_v_hmmer_microcin	cinful_out/01_orf_homology/microcins/blast_v_hmmer.csv
# buildhmm_CvaB	cinful_out/00_dbs/CvaB.verified.hmm
# buildhmm_immunity_protein	cinful_out/00_dbs/immunity_proteins.verified.hmm
# buildhmm_microcin	cinful_out/00_dbs/microcins.verified.hmm
# candidate_immunity	cinful_out/03_best_hits/best_immunity_protein_candidates.csv
# filter_CvaB	cinful_out/01_orf_homology/CvaB/filtered_nr.fa
# filter_CvaB_hits	cinful_out/01_orf_homology/CvaB/preQC.csv cinful_out/01_orf_homology/CvaB/QC.csv
# filter_immunity_protein	cinful_out/01_orf_homology/immunity_proteins/filtered_nr.fa
# filter_microcin	cinful_out/01_orf_homology/microcins/filtered_nr.fa
# final	final
# makeblastdb_CvaB	cinful_out/00_dbs/CvaB.verified.pep.dmnd
# makeblastdb_immunity_protein	cinful_out/00_dbs/immunity_proteins.verified.pep.phr
# makeblastdb_microcin	cinful_out/00_dbs/microcins.verified.pep.phr
# merged_results	cinful_out/02_homology_results/all_merged.csv
# msa_CvaB	cinful_out/00_dbs/CvaB.verified.aln
# msa_immunity_protein	cinful_out/00_dbs/immunity_proteins.verified.aln
# msa_microcin	cinful_out/00_dbs/microcins.verified.aln
# nonredundant_prodigal	cinful_out/01_orf_homology/prodigal_out.all.nr.faa cinful_out/01_orf_homology/prodigal_out.all.nr_expanded.csv
# prodigal	cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.gff3 cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.cds cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.faa
# write_CvaB	cinful_out/00_dbs/CvaB.verified.pep
# write_immunity_proteins	cinful_out/00_dbs/immunity_proteins.verified.pep
# write_microcins	cinful_out/00_dbs/microcins.verified.pep


python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.with_verified.aln -F -j1 --keep-target-files --allowed-rules align_with_verifiedCvab -d .tests/align_with_verifiedCvab  
python -m snakemake cinful_out/03_best_hits/best_hit_contigs.csv -F -j1 --keep-target-files --allowed-rules bestHitsContigs -d .tests/bestHitsContigs  
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.fa -F -j1 --keep-target-files --allowed-rules best_Cvab_fasta -d .tests/best_Cvab_fasta  
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.headers -F -j1 --keep-target-files --allowed-rules best_Cvab_headers -d .tests/best_Cvab_headers  
python -m snakemake cinful_out/03_best_hits/best_hits.csv -F -j1 --keep-target-files --allowed-rules best_hits -d .tests/best_hits  
python -m snakemake cinful_out/01_orf_homology/CvaB/blast.txt -F -j1 --keep-target-files --allowed-rules blast_CvaB -d .tests/blast_CvaB  
python -m snakemake cinful_out/01_orf_homology/immunity_proteins/blast.txt -F -j1 --keep-target-files --allowed-rules blast_immunity_protein -d .tests/blast_immunity_protein  
python -m snakemake cinful_out/01_orf_homology/microcins/blast.txt -F -j1 --keep-target-files --allowed-rules blast_microcin -d .tests/blast_microcin  
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_CvaB -d .tests/blast_v_hmmer_CvaB  
python -m snakemake cinful_out/01_orf_homology/immunity_proteins/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_immunity_protein -d .tests/blast_v_hmmer_immunity_protein  
python -m snakemake cinful_out/01_orf_homology/microcins/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_microcin -d .tests/blast_v_hmmer_microcin  
python -m snakemake cinful_out/00_dbs/CvaB.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_CvaB -d .tests/buildhmm_CvaB  
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_immunity_protein -d .tests/buildhmm_immunity_protein  
python -m snakemake cinful_out/00_dbs/microcins.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_microcin -d .tests/buildhmm_microcin  
python -m snakemake cinful_out/03_best_hits/best_immunity_protein_candidates.csv -F -j1 --keep-target-files --allowed-rules candidate_immunity -d .tests/candidate_immunity  
python -m snakemake cinful_out/01_orf_homology/CvaB/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_CvaB -d .tests/filter_CvaB  
python -m snakemake cinful_out/01_orf_homology/CvaB/preQC.csv cinful_out/01_orf_homology/CvaB/QC.csv -F -j1 --keep-target-files --allowed-rules filter_CvaB_hits -d .tests/filter_CvaB_hits  
python -m snakemake cinful_out/01_orf_homology/immunity_proteins/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_immunity_protein -d .tests/filter_immunity_protein  
python -m snakemake cinful_out/00_dbs/CvaB.verified.pep.dmnd -F -j1 --keep-target-files --allowed-rules makeblastdb_CvaB -d .tests/makeblastdb_CvaB  
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.pep.phr -F -j1 --keep-target-files --allowed-rules makeblastdb_immunity_protein -d .tests/makeblastdb_immunity_protein  
python -m snakemake cinful_out/00_dbs/microcins.verified.pep.phr -F -j1 --keep-target-files --allowed-rules makeblastdb_microcin -d .tests/makeblastdb_microcin  
python -m snakemake cinful_out/02_homology_results/all_merged.csv -F -j1 --keep-target-files --allowed-rules merged_results -d .tests/merged_results  
python -m snakemake cinful_out/00_dbs/CvaB.verified.aln -F -j1 --keep-target-files --allowed-rules msa_CvaB -d .tests/msa_CvaB  
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.aln -F -j1 --keep-target-files --allowed-rules msa_immunity_protein -d .tests/msa_immunity_protein  
python -m snakemake cinful_out/00_dbs/microcins.verified.aln -F -j1 --keep-target-files --allowed-rules msa_microcin -d .tests/msa_microcin  
python -m snakemake cinful_out/01_orf_homology/prodigal_out.all.nr.faa cinful_out/01_orf_homology/prodigal_out.all.nr_expanded.csv -F -j1 --keep-target-files --allowed-rules nonredundant_prodigal -d .tests/nonredundant_prodigal  
python -m snakemake cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.gff3 cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.cds cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.faa -F -j1 --keep-target-files --allowed-rules prodigal -d .tests/prodigal  
python -m snakemake cinful_out/00_dbs/CvaB.verified.pep -F -j1 --keep-target-files --allowed-rules write_CvaB -d .tests/write_CvaB  
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.pep -F -j1 --keep-target-files --allowed-rules write_immunity_proteins -d .tests/write_immunity_proteins  
python -m snakemake cinful_out/00_dbs/microcins.verified.pep -F -j1 --keep-target-files --allowed-rules write_microcins -d --snakefile ~/github.com/tijeco/cinful/cinful/Snakefile .tests/write_microcins  


 memusg -t 