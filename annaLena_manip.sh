contigs_dir="/mnt/chaelab/shared/anna_lena/contigs"
for file in $contigs_dir/*.fasta; do
    accessionNum="$(basename "$file" .contigs.fasta)"
    sed -i "s/>/>$accessionNum./g" "$file"
done

## cat all contigs into a single fasta file to build blast database
all_contigs="$contigs_dir/../anna_lena70.contigs.fasta"
cat $contigs_dir/*.fasta > $all_contigs

db_dir="/mnt/chaelab/shared/anna_lena/blastdb"
makeblastdb -in $all_contigs -dbtype nucl -logfile anna_lena70.contigs.blastdb.log -out anna_lena70.contigs

## blast AT5G58120
TE_SV_dir="/mnt/chaelab/rachelle/TE_SV"
blast_results="$TE_SV_dir/results/anna_lena70_blastn/anna_lena70_AT5G58120.80acc.Col-0.blastn.out"
blastn -db $db_dir/anna_lena70.contigs -query $TE_SV_dir/data/AT5G58120_20kbflank_80acc_Col-0_ungapped.fa -out $blast_results

## get contig hits names
contig_names="$TE_SV_dir/results/anna_lena70_blastn/anna_lena70_AT5G58120.hits.contigs.names.txt"
grep '^> [0-9]' $blast_results | sort | uniq | cut -d ' ' -f2 > $contig_names

## get contig hits sequence
contig_hits="$TE_SV_dir/results/anna_lena70_blastn/anna_lena70_AT5G58120.hits.contigs.fasta"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/TE_SV/src'); import get_contigs; get_contigs.get_contigs('${all_contigs}','${contig_names}','${contig_hits}')"

## cat contigs with Col-0 reference
contig_wCol0="$TE_SV_dir/results/anna_lena70_blastn/anna_lena70_AT5G58120.hits.contigs_Col-0.fasta"
cat $TE_SV_dir/data/AT5G58120_Col-0.fasta $contig_hits > $contig_wCol0

## align and view
contig_aligned="$TE_SV_dir/results/anna_lena70_AT5G58120.hits.contigs_Col-0_mafft.fa"
mafft --large --adjustdirectionaccurately $contig_wCol0 > $contig_aligned

## Col-0 doesn't align as well as the rest of the contigs, trying to revcomp appropriate sequences
rev_comp_list="$TE_SV_dir/results/anna_lena70_blastn/contigs_to_revcomp.txt"
fout_revcomp="$TE_SV_dir/results/anna_lena70_blastn/anna_lena70_AT5G58120.hits.contigs_Col-0.revcompsome.fasta"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/TE_SV/src'); import get_contigs; get_contigs.rc_selected_contigs('${contig_wCol0}','${rev_comp_list}','${fout_revcomp}')"
## realign and view
contig_revcomp_aligned="$TE_SV_dir/results/anna_lena70_AT5G58120.hits.contigs_Col-0.revcompsome_mafft.fa"
mafft --large $fout_revcomp > $contig_revcomp_aligned
## EVEN WORSE OTL


## different blast output
json_dir="$TE_SV_dir/results/anna_lena70_blastn/json"
mkdir $json_dir
blast_results_json="$json_dir/anna_lena70_AT5G58120.80acc.Col-0.blastn_json.out"
blastn -outfmt 13 -db $db_dir/anna_lena70.contigs -query $TE_SV_dir/data/AT5G58120_20kbflank_80acc_Col-0_ungapped.fa -out $blast_results_json

col0_json="$json_dir/anna_lena70_AT5G58120.80acc.Col-0.blastn_json.out_1.json"
mkdir $TE_SV_dir/results/anna_lena70_blastn/json_Col0
col0_json_seqs="$TE_SV_dir/results/anna_lena70_blastn/json_Col0_20kbflank/anna_lena70_AT5G58120.Col-0.blastn_seqs.fasta"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import blast_fasta; hit_ranges=blast_fasta.get_hits_minlen('${col0_json}',0,False); import re; tmp = [(re.findall('\d+\.tig\d+',k)[0],v) for k,v in list(hit_ranges.items())]; tmp_dict = dict(tmp); blast_fasta.get_fastaseq_from_range('${contig_wCol0}',tmp_dict,'${col0_json_seqs}')"

## cat Col-0 sequence and realign
contig_col0="$TE_SV_dir/results/anna_lena70_blastn/json_Col0_20kbflank/anna_lena70_AT5G58120.hits.contigs_Col-0_blastnranges.fa"
cat /mnt/chaelab/rachelle/TE_SV/data/AT5G58120_20kbflank_Col-0.fasta $col0_json_seqs > $contig_col0
contig_col0_aligned="$TE_SV_dir/results/anna_lena70_blastn/json_Col0_20kbflank/anna_lena70_AT5G58120.hits.contigs_Col-0_blastnranges_mafft.fa"
mafft --adjustdirectionaccurately --large $contig_col0 > $contig_col0_aligned


## blast contigs against NLR sequences
## get NLR bedfile entries
bed_tair10="/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_gene.bed"
geneID_nlr="/mnt/chaelab/shared/NLR/geneID/all_NLR.txt"
geneID_col=3 # index starts from 0 for python
bed_nlr="/mnt/chaelab/shared/NLR/bed/all_NLR.bed"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.filter_column('${bed_tair10}','${geneID_nlr}','${geneID_col}','${bed_nlr}')"
bed_chr="/mnt/chaelab/rachelle/data/features/gff_bed_Chr_TAIR10.bed"
bed_nlr_flank="/mnt/chaelab/shared/NLR/bed/all_NLR_10kbflank.bed"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import bed_manip; bed_manip.extend('${bed_nlr}','${bed_nlr_flank}','${bed_chr}',10000)"

## get NLR sequences
fasta_tair10="/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta"
fasta_nlr="/mnt/chaelab/shared/NLR/sequence/all_NLR.fasta"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.intersect_bed_fasta('${bed_nlr}','${fasta_tair10}','${fasta_nlr}')"

fasta_nlr_flank="/mnt/chaelab/shared/NLR/sequence/all_NLR_10kbflank.fasta"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.intersect_bed_fasta('${bed_nlr_flank}','${fasta_tair10}','${fasta_nlr_flank}')"

## make blast database with NLR sequences
nlr_db_dir="/mnt/chaelab/rachelle/blastdb/nlr164"
makeblastdb -in $fasta_nlr -dbtype nucl -logfile ${nlr_db_dir}/all_NLR.164_db.log -out ${nlr_db_dir}/all_NLR.164.db

nlr_db_dir_flank="/mnt/chaelab/rachelle/blastdb/nlr164_10kbflank"
makeblastdb -in $fasta_nlr_flank -dbtype nucl -logfile ${nlr_db_dir_flank}/all_NLR.164.10kbflank_db.log -out ${nlr_db_dir_flank}/all_NLR.164.10kbflank.db

## blast DM10 vector against NLR sequences (to find gene)
dm10="/mnt/chaelab/rachelle/data/NLR/pMLBartRgeneTulscha_DM10.fasta"
dm10_blastn_results="$TE_SV_dir/results/DM10/pMLBartRgeneTulscha_DM10.Col-0.blastn.asn"
blastn -outfmt 11 -db $nlr_db_dir/all_NLR.164.db -query $dm10 -out $dm10_blastn_results
dm10_blast_results_sum="$TE_SV_dir/results/DM10/pMLBartRgeneTulscha_DM10.Col-0.blastn_summary.txt"
blast_formatter -archive $dm10_blastn_results -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen evalue bitscore" > $dm10_blast_results_sum
dm10_blast_results_sum_tsv="$TE_SV_dir/results/DM10/pMLBartRgeneTulscha_DM10.Col-0.blastn_summary.tsv"
printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\tevalue\tbit score\n$(grep -v '^#' ${dm10_blast_results_sum})" > $dm10_blast_results_sum_tsv
## blast DM10 vector against anna lena contigs
dm10="/mnt/chaelab/rachelle/data/NLR/pMLBartRgeneTulscha_DM10only.fasta"
dm10_blast_results="$TE_SV_dir/results/DM10/db_annaLena70/pMLBartRgeneTulscha_DM10.annaLena70.blastn.asn"
blastn -outfmt 11 -db $db_dir/anna_lena70.contigs -query $dm10 -out $dm10_blast_results
dm10_blast_results_sum="$TE_SV_dir/results/DM10/db_annaLena70/pMLBartRgeneTulscha_DM10.annaLena70.blastn_summary.txt"
blast_formatter -archive $dm10_blast_results -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen evalue bitscore" > $dm10_blast_results_sum
dm10_blast_results_sum_tsv="$TE_SV_dir/results/DM10/db_annaLena70/pMLBartRgeneTulscha_DM10.annaLena70.blastn_summary.tsv"
printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\tevalue\tbit score\n$(grep -v '^#' ${dm10_blast_results_sum})" > $dm10_blast_results_sum_tsv

## blast contigs against NLR sequences (unflanked)
al_nlr164_blast="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0"
al_nlr164_blast_results="$al_nlr164_blast/anna_lena70.Col-0.blastn.asn"
blastn -outfmt 11 -db $nlr_db_dir/all_NLR.164.db -query $all_contigs -out $al_nlr164_blast_results
al_nlr164_blast_results_sum="$al_nlr164_blast/anna_lena70.Col-0.blastn_summary.txt"
blast_formatter -archive $al_nlr164_blast_results -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen evalue bitscore" > $al_nlr164_blast_results_sum
al_nlr164_blast_results_sum_tsv="$al_nlr164_blast/anna_lena70.Col-0.blastn_summary.tsv"
printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\tevalue\tbit score\n$(grep -v '^#' ${al_nlr164_blast_results_sum})" > $al_nlr164_blast_results_sum_tsv

## blast contigs against NLR sequences (flanked)
al_nlr164_blast_flank="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0_10kbflank"
al_nlr164_blast_results_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn.asn"
blastn -outfmt 11 -db $nlr_db_dir_flank/all_NLR.164.10kbflank.db -query $all_contigs -out $al_nlr164_blast_results_flank
al_nlr164_blast_results_sum_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_summary.txt"
blast_formatter -archive $al_nlr164_blast_results_flank -outfmt "7 qacc sacc pident length mismatch gapopen qstart qend sstart send qlen evalue bitscore" > $al_nlr164_blast_results_sum_flank
al_nlr164_blast_results_sum_tsv_flank="$al_nlr164_blast/anna_lena70.Col-0_10kbflank.blastn_summary.tsv"
printf "query acc.\tsubject acc.\t%% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tq. len\tevalue\tbit score\n$(grep -v '^#' ${al_nlr164_blast_results_sum_flank})" > $al_nlr164_blast_results_sum_tsv_flank

## retain only best scoring (bitscore) hit per query-subject pair (to reduce dataset size)
al_nlr164_blast_topHits="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topHits.tsv"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.keep_best_hit_per_query_subject_pair('${al_nlr164_blast_results_sum_tsv}','${al_nlr164_blast_topHits}',-1,min_length=500)"

al_nlr164_blast_topHits_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topHits.tsv"
python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.keep_best_hit_per_query_subject_pair('${al_nlr164_blast_results_sum_tsv_flank}','${al_nlr164_blast_topHits_flank}',-1,min_length=500)"

# ## keep only query-subject pairs where query has hit in AT5G58120
# al_nlr164_blast_at5g58120="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topHits_AT5G58120.tsv"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.keep_query_aligned_to_selected_subjects('${al_nlr164_blast_topHits}','${al_nlr164_blast_at5g58120}',['AT5G58120_Chr5:23517491-23521067'])"

# al_nlr164_blast_at5g58120_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10flank.blastn_topHits_AT5G58120.tsv"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import file_manip; file_manip.keep_query_aligned_to_selected_subjects('${al_nlr164_blast_topHits_flank}','${al_nlr164_blast_at5g58120_flank}',['AT5G58120_Chr5:23507491-23531067'])"

# ## manually identify high scoring contigs (using LibreOffice Calc lol) for alignment and further analysis on AT5G58120

# ## get sequences of high scoring contigs with hit in AT5G58120
# al_nlr164_blast_at5g58120_contigs_seq="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topContigs_AT5G58120.fasta"
# top_contig_titles="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topContigs_AT5G58120.txt"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.filter_title('${

# all_contigs}','${al_nlr164_blast_at5g58120_contigs_seq}','${top_contig_titles}')"
# al_nlr164_blast_at5g58120_contigs_seq_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_AT5G58120.fasta"
# top_contig_titles_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_AT5G58120.txt"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.filter_title('${all_contigs}','${al_nlr164_blast_at5g58120_contigs_seq_flank}','${top_contig_titles_flank}')"

# ## align for viewing
# al_contigs_at5g58120_col0="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topContigs_Col-0_AT5G58120.fasta"
# al_contigs_at5g58120_col0_aligned="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topContigs_Col-0_AT5G58120_mafft.fa"
# cat $TE_SV_dir/data/AT5G58120_Col-0.fasta $al_nlr164_blast_at5g58120_contigs_seq > $al_contigs_at5g58120_col0
# mafft --adjustdirectionaccurately --large $al_contigs_at5g58120_col0 > $al_contigs_at5g58120_col0_aligned

# al_contigs_at5g58120_col0_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120.fasta"
# al_contigs_at5g58120_col0_aligned_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft.fa"
# cat $TE_SV_dir/data/AT5G58120_10kbflank_Col-0.fasta $al_nlr164_blast_at5g58120_contigs_seq_flank > $al_contigs_at5g58120_col0_flank
# mafft --adjustdirectionaccurately --large $al_contigs_at5g58120_col0_flank > $al_contigs_at5g58120_col0_aligned_flank

# ## trim alignment to reference
# al_contigs_at5g58120_col0_aligned_trim="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0/anna_lena70.Col-0.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed.fa"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.trim_alignment('${al_contigs_at5g58120_col0_aligned}','NC_003076','${al_contigs_at5g58120_col0_aligned_trim}')"

# al_contigs_at5g58120_col0_aligned_trim_flank="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed.fa"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.trim_alignment('${al_contigs_at5g58120_col0_aligned_flank}','NC_003076','${al_contigs_at5g58120_col0_aligned_trim_flank}')"

# ## keep only CDS
# gff_bed="/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed"
# geneID="AT5G58120"
# start_pos="23507255"
# ref_title="NC_003076"
# al_contigs_at5g58120_col0_aligned_cds="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_CDSonly.fa"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.extract_CDS('${al_contigs_at5g58120_col0_aligned_trim_flank}','${gff_bed}','${ref_title}','${geneID}',int('${start_pos}'),'${al_contigs_at5g58120_col0_aligned_cds}')"
# ## translate CDS and realign
# al_contigs_at5g58120_col0_cds_translated="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_CDStranslated.fasta"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.translate_aligned_CDS('${al_contigs_at5g58120_col0_aligned_cds}','${ref_title}','${al_contigs_at5g58120_col0_cds_translated}', stop_symbol='.')"
# al_contigs_at5g58120_col0_cds_translated_aligned="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_CDStranslated_mafft.fa"
# mafft --large $al_contigs_at5g58120_col0_cds_translated > $al_contigs_at5g58120_col0_cds_translated_aligned

# ## keep entire gene including introns
# gff_bed="/mnt/chaelab/shared/genomes/TAIR10/features/TAIR10_GFF3_genes.bed"
# geneID="AT5G58120"
# start_pos="23507255"
# ref_title="NC_003076"
# al_contigs_at5g58120_col0_aligned_gene="$al_nlr164_blast_flank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_gene.fa"
# python3 -c "import os; os.chdir('/mnt/chaelab/rachelle/src'); import fasta_manip; fasta_manip.extract_feature('${al_contigs_at5g58120_col0_aligned_trim_flank}','${gff_bed}','${ref_title}','${geneID}',int('${start_pos}'),'${al_contigs_at5g58120_col0_aligned_gene}', 'gene')"

# ## find repeats
# fa_10kb="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0_10kbflank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed.fa"
# fa_ungapped="/mnt/chaelab/rachelle/TE_SV/results/anna_lena70_blastn/json_Col0_10kbflank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed_ungap.fa"
# sed 's/-//g' $fa_10kb > $fa_ungap
# ## manually put the '-' back into Col-0's coordinates
# repeat_out_dir="/mnt/chaelab/rachelle/TE_SV/results/RepeatMasker_20190531_AT5G58120_10kbflank_annalena70"
# RepeatMasker -species arabidopsis -dir $repeat_out_dir $fa_ungapped
# repeat_out="/mnt/chaelab/rachelle/TE_SV/results/RepeatMasker_20190531_AT5G58120_10kbflank_annalena70/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed_ungap.fa.out"
# repeat_out_tsv="/mnt/chaelab/rachelle/TE_SV/results/RepeatMasker_20190531_AT5G58120_10kbflank_annalena70/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed_ungap.fa.out.tsv"
# awk 'FNR > 3' $repeat_out | tr -s ' ' | tr ' ' '\t' > $repeat_out_tsv
# ## files were copied over to desktop and visualised with R using the following commands
# # RMresults = "/home/rachelle/cey/TE_SV/results/RepeatMasker_20190531_AT5G58120_10kbflank_annalena70/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed_ungap.fa.out.tsv" ## $repeat_out_tsv
# # fasta = "/home/rachelle/cey/TE_SV/results/anna_lena70_blastn/json_Col0_10kbflank/anna_lena70.Col-0_10kbflank.blastn_topContigs_Col-0_AT5G58120_mafft_trimmed.fa"
# # plot_RepeatMasker(RMresults, fasta, legend_size = 0.2) ## from /home/rachelle/cey/TE_SV/src/TE_SV_plot.R


# ## TODO: visualise
