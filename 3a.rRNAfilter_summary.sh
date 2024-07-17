#!/bin/bash
# generate summary of rRNA filtering

cd /vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/scripts/

SOURCE_DIR='../logs/rRNA_globinfilter/'
LOGFILES=("$SOURCE_DIR"*.log)


echo -e 'file\t# read pairs\t# concordant once\t% concordant once\t# concordant more than once\t% concordant more than once\t# of aligned discordant\t% of aligned disdorant\t# of mates condordant once\t% of mates concordant once\t# of mates concordant more than once\t% of mates concordant more than once\toverall % alignment' > rRNA_summary_tg.txt

echo "Number of files: ${#LOGFILES[*]}"

for item in ${LOGFILES[@]}
do
        printf "   %s\n" $item

	#item="$item"_rRNA_globinfilter
        OUTPUT1=$(cat $item | grep -E 'reads' | sed 's/ reads.*//')
        OUTPUT2=$(cat $item | grep -E 'concordantly exactly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT3=$(cat $item | grep -E 'concordantly exactly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT4=$(cat $item | grep -E 'aligned concordantly >1 times' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT5=$(cat $item | grep -E 'aligned concordantly >1 times' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT6=$(cat $item | grep -E 'aligned discordantly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT7=$(cat $item | grep -E 'aligned discordantly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT8=$(cat $item | grep -E 'aligned exactly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT9=$(cat $item | grep -E 'aligned exactly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT10=$(cat $item | grep -E 'aligned >1 times' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT11=$(cat $item | grep -E 'aligned >1 times' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT12=$(cat $item | grep -E 'overall alignment rate' | sed 's/%.*//')

        echo -e ''$item'\t'${OUTPUT1}'\t'${OUTPUT2}'\t'${OUTPUT3}'\t'${OUTPUT4}'\t'${OUTPUT5}'\t'${OUTPUT6}'\t'${OUTPUT7}'\t'${OUTPUT8}'\t'${OUTPUT9}'\t'${OUTPUT10}'\t'${OUTPUT11}'\t'${OUTPUT12}'' >> rRNA_summary_tg.txt
         
done
