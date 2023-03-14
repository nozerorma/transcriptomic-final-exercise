#### PRE-PROCESSING SCRIPT ####

f_path=$2
r_path=$3
outdir=$4
logdir=$5

f_name=$(basename "$f_path")
r_name=$(basename "$r_path")
base_sid=$(basename "$f_path" | cut -d"_" -f1)
f_sid=$(basename "$f_name" .fastq)
r_sid=$(basename "$r_name" .fastq)

if [ "$1" == "cutadapt" ]; then
    
    echo -e "\nRunning cutadapt...\n" 
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Trim already performed for $f_sid, skipping...\n"
    
    else
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
            -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
            -a CTGTCTCTTATACACATCT...AGATGTGTATAAGAGACAG -a TGGAATTCTCGGGTGCCAAGG \
            -a "G{100}" -g "G{100}" -a "A{100}" -g "A{100}" -q 20 \
            -o "$outdir/${f_sid}_trimmed.fastq" -p "$outdir/${r_sid}_trimmed.fastq" \
            "$f_path" "$r_path" --cores 6 > "$logdir"/log.txt
    fi
fi

read -rp "Would you like to re-run QC repor for your trimmed samples? (Y/n) " runqc

for trimmed_sid in $(find $outdir -type f -name "*_trimmed.fastq")
    case $runqc in
        [Yy]* )
            scripts/qc.sh $trimmed_sid "out/qc/fastqc" "out/qc/fastq_screen"
        ;;
        [Nn]* )
            echo -e "\nSkipping QC analysis...\n" 
        ;;
	esac
done


# Uncomment and modify conveniently for using Trimmomatic
# elif [ "$1" == "trimmomatic" ]; then
    
#     trimmomatic PE -phred33 "$f_path" "$r_path" \
#             "$outdir"/"$f_name" "$outdir"/"$f_name"_unpaired \
#             "$outdir"/"$r_name" "$outdir"/"$r_name"_unpaired \
#             TRAILING:30 SLIDINGWINDOW:4:30 >> "$logdir"/log.txt


# 