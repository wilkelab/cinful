
mkdir .tests/align_with_verifiedCvab
cp -r .tests/unit/align_with_verifiedCvab/data/ .tests/align_with_verifiedCvab/
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.with_verified.aln -F -j1 --keep-target-files --allowed-rules align_with_verifiedCvab -d .tests/align_with_verifiedCvab

mkdir .tests/bestHitsContigs
cp -r .tests/unit/bestHitsContigs/data/ .tests/bestHitsContigs/
python -m snakemake cinful_out/03_best_hits/best_hit_contigs.csv -F -j1 --keep-target-files --allowed-rules bestHitsContigs -d .tests/bestHitsContigs

mkdir .tests/best_Cvab_fasta
cp -r .tests/unit/best_Cvab_fasta/data/ .tests/best_Cvab_fasta/
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.fa -F -j1 --keep-target-files --allowed-rules best_Cvab_fasta -d .tests/best_Cvab_fasta

mkdir .tests/best_Cvab_headers
cp -r .tests/unit/best_Cvab_headers/data/ .tests/best_Cvab_headers/
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.headers -F -j1 --keep-target-files --allowed-rules best_Cvab_headers -d .tests/best_Cvab_headers

mkdir .tests/best_hits
cp -r .tests/unit/best_hits/data/ .tests/best_hits/
python -m snakemake cinful_out/03_best_hits/best_hits.csv -F -j1 --keep-target-files --allowed-rules best_hits -d .tests/best_hits

mkdir .tests/blast_CvaB
cp -r .tests/unit/blast_CvaB/data/ .tests/blast_CvaB/
python -m snakemake cinful_out/01_orf_homology/CvaB/blast.txt -F -j1 --keep-target-files --allowed-rules blast_CvaB -d .tests/blast_CvaB

## This one failed, said it couldn't find the database file
## Works now, just needed to add the other blast files that are too annoying to add to the rule output
mkdir .tests/blast_immunity_protein
cp -r .tests/unit/blast_immunity_protein/data/ .tests/blast_immunity_protein/
python -m snakemake cinful_out/01_orf_homology/immunity_proteins/blast.txt -F -j1 --keep-target-files --allowed-rules blast_immunity_protein -d .tests/blast_immunity_protein

mkdir .tests/blast_microcin
cp -r .tests/unit/blast_microcin/data/ .tests/blast_microcin/
python -m snakemake cinful_out/01_orf_homology/microcins/blast.txt -F -j1 --keep-target-files --allowed-rules blast_microcin -d .tests/blast_microcin

mkdir .tests/blast_v_hmmer_CvaB
cp -r .tests/unit/blast_v_hmmer_CvaB/data/ .tests/blast_v_hmmer_CvaB/
python -m snakemake cinful_out/01_orf_homology/CvaB/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_CvaB -d .tests/blast_v_hmmer_CvaB

mkdir .tests/blast_v_hmmer_immunity_protein
cp -r .tests/unit/blast_v_hmmer_immunity_protein/data/ .tests/blast_v_hmmer_immunity_protein/
python -m snakemake cinful_out/01_orf_homology/immunity_proteins/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_immunity_protein -d .tests/blast_v_hmmer_immunity_protein

mkdir .tests/blast_v_hmmer_microcin
cp -r .tests/unit/blast_v_hmmer_microcin/data/ .tests/blast_v_hmmer_microcin/
python -m snakemake cinful_out/01_orf_homology/microcins/blast_v_hmmer.csv -F -j1 --keep-target-files --allowed-rules blast_v_hmmer_microcin -d .tests/blast_v_hmmer_microcin

mkdir .tests/buildhmm_CvaB
cp -r .tests/unit/buildhmm_CvaB/data/ .tests/buildhmm_CvaB/
python -m snakemake cinful_out/00_dbs/CvaB.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_CvaB -d .tests/buildhmm_CvaB

mkdir .tests/buildhmm_immunity_protein
cp -r .tests/unit/buildhmm_immunity_protein/data/ .tests/buildhmm_immunity_protein/
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_immunity_protein -d .tests/buildhmm_immunity_protein

mkdir .tests/buildhmm_microcin
cp -r .tests/unit/buildhmm_microcin/data/ .tests/buildhmm_microcin/
python -m snakemake cinful_out/00_dbs/microcins.verified.hmm -F -j1 --keep-target-files --allowed-rules buildhmm_microcin -d .tests/buildhmm_microcin

mkdir .tests/candidate_immunity
cp -r .tests/unit/candidate_immunity/data/ .tests/candidate_immunity/
python -m snakemake cinful_out/03_best_hits/best_immunity_protein_candidates.csv -F -j1 --keep-target-files --allowed-rules candidate_immunity -d .tests/candidate_immunity

mkdir .tests/filter_CvaB
cp -r .tests/unit/filter_CvaB/data/ .tests/filter_CvaB/
python -m snakemake cinful_out/01_orf_homology/CvaB/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_CvaB -d .tests/filter_CvaB

mkdir .tests/filter_CvaB_hits
cp -r .tests/unit/filter_CvaB_hits/data/ .tests/filter_CvaB_hits/
python -m snakemake cinful_out/01_orf_homology/CvaB/preQC.csv cinful_out/01_orf_homology/CvaB/QC.csv -F -j1 --keep-target-files --allowed-rules filter_CvaB_hits -d .tests/filter_CvaB_hits

mkdir .tests/filter_immunity_protein
cp -r .tests/unit/filter_immunity_protein/data/ .tests/filter_immunity_protein/
python -m snakemake cinful_out/01_orf_homology/immunity_proteins/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_immunity_protein -d .tests/filter_immunity_protein

# mkdir .tests/filter_microcin
# cp -r .tests/unit/filter_microcin/data/ .tests/filter_microcin/
# python -m snakemake cinful_out/01_orf_homology/microcins/filtered_nr.fa -F -j1 --keep-target-files --allowed-rules filter_microcin -d .tests/filter_microcin

# mkdir .tests/final
# cp -r .tests/unit/final/data/ .tests/final/
# python -m snakemake final -F -j1 --keep-target-files --allowed-rules final -d .tests/final

mkdir .tests/makeblastdb_CvaB
cp -r .tests/unit/makeblastdb_CvaB/data/ .tests/makeblastdb_CvaB/
python -m snakemake cinful_out/00_dbs/CvaB.verified.pep.dmnd -F -j1 --keep-target-files --allowed-rules makeblastdb_CvaB -d .tests/makeblastdb_CvaB

mkdir .tests/makeblastdb_immunity_protein
cp -r .tests/unit/makeblastdb_immunity_protein/data/ .tests/makeblastdb_immunity_protein/
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.pep.phr -F -j1 --keep-target-files --allowed-rules makeblastdb_immunity_protein -d .tests/makeblastdb_immunity_protein

mkdir .tests/makeblastdb_microcin
cp -r .tests/unit/makeblastdb_microcin/data/ .tests/makeblastdb_microcin/
python -m snakemake cinful_out/00_dbs/microcins.verified.pep.phr -F -j1 --keep-target-files --allowed-rules makeblastdb_microcin -d .tests/makeblastdb_microcin

mkdir .tests/merged_results
cp -r .tests/unit/merged_results/data/ .tests/merged_results/
python -m snakemake cinful_out/02_homology_results/all_merged.csv -F -j1 --keep-target-files --allowed-rules merged_results -d .tests/merged_results

mkdir .tests/msa_CvaB
cp -r .tests/unit/msa_CvaB/data/ .tests/msa_CvaB/
python -m snakemake cinful_out/00_dbs/CvaB.verified.aln -F -j1 --keep-target-files --allowed-rules msa_CvaB -d .tests/msa_CvaB

mkdir .tests/msa_immunity_protein
cp -r .tests/unit/msa_immunity_protein/data/ .tests/msa_immunity_protein/
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.aln -F -j1 --keep-target-files --allowed-rules msa_immunity_protein -d .tests/msa_immunity_protein

mkdir .tests/msa_microcin
cp -r .tests/unit/msa_microcin/data/ .tests/msa_microcin/
python -m snakemake cinful_out/00_dbs/microcins.verified.aln -F -j1 --keep-target-files --allowed-rules msa_microcin -d .tests/msa_microcin

# mkdir .tests/nonredundant_prodigal
# cp -r .tests/unit/nonredundant_prodigal/data/ .tests/nonredundant_prodigal/
# python -m snakemake cinful_out/01_orf_homology/prodigal_out.all.nr.faa cinful_out/01_orf_homology/prodigal_out.all.nr_expanded.csv -F -j1 --keep-target-files --allowed-rules nonredundant_prodigal -d .tests/nonredundant_prodigal

mkdir .tests/prodigal
cp -r .tests/unit/prodigal/data/ .tests/prodigal/
python -m snakemake cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.gff3 cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.cds cinful_out/01_orf_homology/prodigal_out/GCF_000026285.1_ASM2628v2_genomic.faa -F -j1 --keep-target-files --allowed-rules prodigal -d .tests/prodigal

mkdir .tests/write_CvaB
cp -r .tests/unit/write_CvaB/data/ .tests/write_CvaB/
python -m snakemake cinful_out/00_dbs/CvaB.verified.pep -F -j1 --keep-target-files --allowed-rules write_CvaB -d .tests/write_CvaB

mkdir .tests/write_immunity_proteins
cp -r .tests/unit/write_immunity_proteins/data/ .tests/write_immunity_proteins/
python -m snakemake cinful_out/00_dbs/immunity_proteins.verified.pep -F -j1 --keep-target-files --allowed-rules write_immunity_proteins -d .tests/write_immunity_proteins

mkdir .tests/write_microcins
cp -r .tests/unit/write_microcins/data/ .tests/write_microcins/
python -m snakemake cinful_out/00_dbs/microcins.verified.pep -F -j1 --keep-target-files --allowed-rules write_microcins -d .tests/write_microcins


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
python -m snakemake cinful_out/00_dbs/microcins.verified.pep -F -j1 --keep-target-files --allowed-rules write_microcins -d .tests/write_microcins  