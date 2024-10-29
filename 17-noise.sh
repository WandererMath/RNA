#export ref="../ref_files_for_genome_viewer/genomic.gtf"
#export data_folder="../bowtie_output/"

export ref=$1
export data_folder=$2
export out_folder=$3

echo $1 $2 $3

for file in "$data_folder/"*; do 
    featureCounts -s 2 -M -Q 20 -T 5 -t gene -g gene_id -a $ref -o "$out_folder/$(basename $file)" $file &

done

wait

#featureCounts -M -Q 20 -T 5 -t gene -g gene_id -a $ref -o $filename.txt $file