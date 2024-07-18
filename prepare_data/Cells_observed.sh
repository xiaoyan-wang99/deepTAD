#!/bin/bash

juicer_tool="/home/wangxiaoyan/deepTAD/prepare_data/juicer_tools.jar"
chromList=($(seq 1 22))
chromList[${#chromList[*]}]=X
resolutions="10000 25000 50000 100000"

DPATH="/public_data/wulili/data/HiC/Rwa_Cell/"
CELL="NHEK/GSE63525_NHEK_combined_30.hic IMR90/GSE63525_IMR90_combined_30.hic HMEC/GSE63525_HMEC_combined_30.hic CH12-LX/GSE63525_CH12-LX_combined_30.hic HUVEC/GSE63525_HUVEC_combined_30.hic K562/GSE63525_K562_combined_30.hic"

outputDir="/home/wangxiaoyan/deepTAD/prepare_data"

for cell in $CELL; do
    cell_name=$(echo $cell | cut -d'/' -f1)
    echo $cell_name
    
    #hic_file=$(echo $cell | cut -d'/' -f2)
    cell_outputDir="${outputDir}/${cell_name}"
    mkdir -p $cell_outputDir
    for resolution in $resolutions; do
        display_reso=$((resolution / 1000))
        mkdir -p "${cell_outputDir}/${display_reso}kb"
        for chrom in ${chromList[@]}; do
            java -jar $juicer_tool dump observed KR $DPATH$cell chr$chrom chr$chrom BP $resolution ${cell_outputDir}/${display_reso}kb/${display_reso}k_KR.chr${chrom}_tmp -d
            python remove_nan.py ${cell_outputDir}/${display_reso}kb/${display_reso}k_KR.chr${chrom}_tmp ${cell_outputDir}/${display_reso}kb/${display_reso}k_KR.chr${chrom}
        done
    done
done