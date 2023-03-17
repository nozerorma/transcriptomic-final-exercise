#### PRE-PROCESSING SCRIPT ####

source "scripts/spinner.sh"
RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

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
        # adapters infered from most common Illumina, listed in here
        # https://tinyurl.com/illumina-adapters
        # https://tinyurl.com/illumina-adapters-2
        (cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
            -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
            -a CTGTCTCTTATACACATCT...AGATGTGTATAAGAGACAG -a TGGAATTCTCGGGTGCCAAGG \
            -a "G{100}" -g "G{100}" -a "A{100}" -g "A{100}" -q 20 \
            -o "$outdir/${f_sid}_trimmed.fastq" -p "$outdir/${r_sid}_trimmed.fastq" \
            "$f_path" "$r_path" --cores 14 -m 1 --discard-untrimmed > "$logdir"/log.txt)
    fi
fi

# Re-run quality control
read -rp "Would you like to re-run QC report for your trimmed samples? (Y/n) " runqc

for trimmed_sid in $(find $outdir -type f -name "*_trimmed.fastq");do
    case $runqc in
        [Yy]* )
            bash scripts/qc.sh $trimmed_sid "out/qc/fastqc" "out/qc/fastq_screen"
        ;;
	esac
done

case $runqc in
[Yy]* )
	printf "${GREEN}\nNow, take your time to give a look to the QC analysis.\n${NC}"
	echo "When you are ready, press any key to continue..."
	read -n 1 -s -r -p ""
;;
esac

# Uncomment and modify conveniently for using Trimmomatic
# elif [ "$1" == "trimmomatic" ]; then
    
#     trimmomatic PE -phred33 "$f_path" "$r_path" \
#             "$outdir"/"$f_name" "$outdir"/"$f_name"_unpaired \
#             "$outdir"/"$r_name" "$outdir"/"$r_name"_unpaired \
#             TRAILING:30 SLIDINGWINDOW:4:30 >> "$logdir"/log.txt


# 

