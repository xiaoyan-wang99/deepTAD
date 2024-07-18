#!/bin/bash 
chromList=($(seq 20 22))
resolutions=(10000 25000 50000 100000)

for chrom in ${chromList[@]}; do
for resolution in ${resolutions[@]}; do
  display_reso=`expr $(($resolution/1000))`
  window_size=$(echo "scale=0; 200 / $display_reso" | bc)
  Rscript assemble_filterTAD_run.r -i /home/wangxiaoyan/deepTAD/work/deepTAD/KR_normalization/DATA/HIC002_deepTAD_${display_reso}k.chr${chrom} -p /home/wangxiaoyan/deepTAD/work_public2/cnn_transformer_predict_results_boundary/hic002_${display_reso}k_chr${chrom}-predbound.txt -o /home/wangxiaoyan/deepTAD/work/deepTAD/KR_normalization/DOMAINS/HIC002_deepTAD_${display_reso}k.chr${chrom}  -w  ${window_size}

done
done


