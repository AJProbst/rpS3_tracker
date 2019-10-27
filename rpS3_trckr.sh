#!/bin/bash

if [ "$#" -lt 1 ]
then
  echo "usage: rpS3_trckr.sh <list>"
  echo "where <list> is a /n delimited list of paths to prodigal-predicted gene files, e.g., sample_min1000.fa.genes.faa but the folder containing this file ending with "_min1000.fa.genes.faa" should also contain the original assembly as sample_min1000.fa. The headers of the assembly file need to be specifically formatted:"
  echo ">sequence_name read_length_XXX read_count_YYY"
  echo "where XXX is the read length in bps and YYY is the number of reads mapped."
  echo "output folder is current folder"
  exit 1
fi

# path to tlbx
tblx=/home/ajp/apps/fasta
# path to folder with HMMs
hmm=/home/ajp/apps/.databases/rpS3_HMM

# Input variable
die_liste=$1

## ++++++++++++++ FUNCTIONS ++++++++++++++ ##
# 1. Extract S3
extr_rp(){
  
  # retrive sample name from file
  sample_name=$(ls $sample | awk -F \/ '{print$NF}')
  echo -e "... working on sample ${sample_name} ..."
  for domain in Bacteria Archaea Eukaryotes; do
    hmmsearch --tblout /dev/stdout -o /dev/null --cpu 1 --notextw ${hmm}/${domain}_90_trimmed.hmm ${sample} | grep -v "^#" > ${sample_name}.${domain}.hmm_results.txt
    awk '{print$1}' ${sample_name}.$domain.hmm_results.txt >> ${sample_name}.$domain.all_hits.txt
    pullseq -i ${sample} -n ${sample_name}.$domain.all_hits.txt > ${sample_name}.$domain.all_hits.fasta
    ruby ${tblx}/remove_linebreaks_from_fasta.rb -f ${sample_name}.$domain.all_hits.fasta > ${sample_name}.$domain.all_hits_woLB.fasta
    for gene in $(cat ${sample_name}.$domain.all_hits.txt); do grep -w -A 1 $gene ${sample_name}.$domain.all_hits_woLB.fasta | sed 1d | wc; done | awk '{print$3}' > ${sample_name}.$domain.length.tmp
    awk '{print$6}' ${sample_name}.$domain.hmm_results.txt > ${sample_name}.$domain.scores.tmp
    paste ${sample_name}.$domain.all_hits.txt ${sample_name}.$domain.scores.tmp ${sample_name}.$domain.length.tmp > ${sample_name}.${domain}_hits.txt
  done

  # filter tables by score and length
    # BACTERIA
  awk '{ if ($2 >= 111) print $1 }' ${sample_name}.Bacteria_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}.Bacteria_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}.Bacteria_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 > ${sample_name}_final_rpS3.ids
    # ARCHAEA
  awk '{ if ($2 >= 172) print $1 }' ${sample_name}.Archaea_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}.Archaea_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}.Archaea_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 >> ${sample_name}_final_rpS3.ids
    # EUKARYOTES
  awk '{ if ($2 >= 175) print $1 }' ${sample_name}.Eukaryotes_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}.Eukaryotes_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}.Eukaryotes_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 >> ${sample_name}_final_rpS3.ids
  # retrieve final set of rpS3 genes for sample
  cat ${sample_name}_final_rpS3.ids | pullseq -i ${sample} -N > ${sample_name}_final_rpS3.fa
  # clean up
  rm ${sample_name}.*.all_hits.fasta ${sample_name}.*.all_hits_woLB.fasta ${sample_name}.*.scores.tmp ${sample_name}.*.length.tmp ${sample_name}*_hits.txt ${sample_name}_selected1 ${sample_name}_selected2 ${sample_name}_selected3 ${sample_name}_crit1 ${sample_name}_final_rpS3.ids #${sample_name}.*.all_hits.txt
  echo -e "... done with ${sample_name} ..."
}

# 2. Extract scaffolds
extr_scaff(){
  scaff_file=$(ls $sample | sed "s/\.genes\.faa//g")
  pullseq -i ${scaff_file} -n all_rpS3_scaffold.ids >> all_rpS3_scaffold.fna
}

# 3. Extract clusters
extr_cluster(){
  grep "^>" ${cluster} | awk '{print$1}' | awk -F_ 'NF{NF-=1};1' | sed "s/\ /_/g" | sed "s/>//g" | pullseq -i ../all_rpS3_scaffold.fna -N > ${cluster}_scaffolds.fna
  usearch -sortbylength ${cluster}_scaffolds.fna -fastaout ${cluster}_scaffolds.sorted.faa
  grep -m 1 ">" ${cluster}_scaffolds.sorted.faa | sed "s/>//g" | pullseq -i ${cluster}_scaffolds.sorted.faa -N >> ../final_rpS3_scaffolds.fna
}

## ++++++++++++++ END OF FUNCTIONS ++++++++++++++ ##

# run rpS3 search
for sample in $(cat $die_liste); do
  extr_rp
done

# determine the best clustering parameter and run the culstering
cat *_final_rpS3.fa > all_rpS3.fa
mkdir -p rpS3_clusters
usearch -sortbylength all_rpS3.fa -fastaout all_rpS3_sorted.faa
max=$(ruby ${tblx}/fasta_length_individual.rb -f all_rpS3.fa | sort -k1,1nr | head -n 1)
min=$(ruby ${tblx}/fasta_length_individual.rb -f all_rpS3.fa | sort -k1,1nr | tail -n 1)
var=$(echo "scale=3; $min/$max" | bc)
usearch -cluster_fast all_rpS3_sorted.faa -id 0.99 -centroids all_rps3_centroids.faa -uc all_rps3_99_clusters.uc -clusters rpS3_clusters/rpS3_ -query_cov $var -target_cov $var

# extract scaffolds
echo -e "... extracting scaffolds with rpS3 genes ..."
grep ">" all_rpS3.fa | sed "s/>//g" | awk '{print$1}' | awk -F_ 'NF{NF-=1};1' | sed "s/\ /_/g" > all_rpS3_scaffold.ids
for sample in $(cat $die_liste); do
  extr_scaff
done
rm all_rpS3_scaffold.ids

# extract scaffolds of clusters
cd rpS3_clusters
ls -1 rpS3_* > all_genes.txt
split all_genes.txt -l 6 batch. -d -a 5
for batch in $(ls -1 batch.*); do
  for cluster in $(cat $batch); do
    extr_cluster
  done
  wait
done
rm all_genes.txt batch.*
ruby ${tblx}/ScaffCov2OTU.rb
cd ..

# create a scaff2OTU file
awk -F\, '{print$1}' OTU_table_from_assembly_coverage.csv | sed 1d > OTUs
for i in $(cat OTUs); do head -n 1 rpS3_clusters/${i}_scaffolds.sorted.faa | sed "s/>//g" | awk '{print$1}' >> scaffs; done
paste scaffs OTUs > scaff2cluster.txt
rm OTUs scaffs

echo -e "... all done! output file is final_rpS3_scaffolds.fna"
echo -e "... there also an OTU clustering table from scaffold abundances as given in the assembly: OTU_table_from_assembly_coverage.csv"
echo -e "... there is also a scaff2OTU file called scaff2bin.txt"

# store some intermediate files
mkdir -p intermediate_files
mv *_final_rpS3.fa intermediate_files/.
mv *hmm_results.txt intermediate_files/.
mv all_rpS3_sorted.faa intermediate_files/.
mv all_rps3_centroids.faa intermediate_files/.
mv all_rps3_99_clusters.uc intermediate_files/.
mv all_rpS3_scaffold.fna intermediate_files/.
