##
RSEM_mm10_index=/path/RSEM_mm10_GRCm38.101_index
star_index=/path/RSEM_mm10_GRCm38.101_index
star=/path/star-2.7.6a-0/bin
tri_data=/path/
align=/path/
for ((i=1;i<=6;i++));
do
rsem-calculate-expression \
				        --paired-end \
                                         -p 16  \
                                         --estimate-rspd \
                                         --append-names \
                                         --star --star-path $star \
                                         $tri_data/M${i}_clean_1P.fq  \
                                         $tri_data/M${i}_clean_2P.fq   \
                                         $star_index/mouse_star_ref \
                                         $align/STAR_M${i}
done              
