rowdata=/path/clean_data
output=/path/tri_clean_data
for ((i=1;i<=6;i++));
do
trimmomatic PE -phred33 $rowdata/M${i}*R1.fq.gz $rowdata/M${i}*R2.fq.gz -baseout $output/M${i}_clean.fq.gz  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 HEADCROP:11 
done
