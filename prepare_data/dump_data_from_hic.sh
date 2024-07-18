#!/bin/bash

checkMakeDirectory(){
        echo -e "checking directory: $1"
        if [ ! -e "$1" ]; then
                echo -e "\tmakedir $1"
                mkdir -p "$1"
        fi
}

dumpdata() {

    list=($(seq $1 $2))
    chromList=($(seq 1 22))
    chromList[${#chromList[*]}]=X

    for li in ${list[@]}; do
        j=$((1551550+li-1))
        previous_name="GSM"$j
        latter_name=$(printf "HIC%03d" $li)
        dataset=${previous_name}_${latter_name}.hic
        #echo $dataset
        resolution=(10000 25000 50000 100000)
        mkdir -p $latter_name

        for chrom in ${chromList[@]}; do
            for reso in ${resolution[@]}; do
                display_reso=$((reso / 1000)
                java -jar juicer_tools.jar dump observed KR ./${dataset} chr${chrom} chr${chrom} BP ${reso} ${latter_name}/${latter_name}_${display_reso}k_KR.chr${chrom}_tmp -d
                python remove_nan.py ${latter_name}/${latter_name}_${display_reso}k_KR.chr${chrom}_tmp ${latter_name}/${latter_name}_${display_reso}k_KR.chr${chrom}
            done
        done
    done
}
dumpdata 2 2
dumpdata  50 56
dumpdata  69 74

