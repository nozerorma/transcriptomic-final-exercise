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

# Comprobation for independent use of script
if [ "$#" -ne 5 ]
then
    printf "${RED}Usage: $1 <tool> $2 <read1> $3 <read2> $4 <outdir> $5 <logdir>${NC}\n"
    echo -e 'tool: "cutadapt", "fastp"\n'
	exit 1
fi

# Cutadapt
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
            -a "G{10}" -A "G{10}" -g "G{10}" -G "G{10}" -g "A{10}" -G "A{10}" -q 15 \
            -o "$outdir/${f_sid}_trimmed.fastq" -p "$outdir/${r_sid}_trimmed.fastq" \
            "$f_path" "$r_path" --cores 14 -m 15 --discard-untrimmed > "$logdir"/log.txt)
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
fi

# FASTP (no need to re-run QC as it already performs it. Much powerful tool although it may miss some adapters)
if [ "$1" == "fastp" ]; then
    
    echo -e "\nRunning fastp...\n" 
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Pre-processing already performed for $f_sid, skipping...\n"
    
    else
        (fastp -i $f_path -o "$outdir/${f_sid}_trimmed.fastq" -I $r_path -O "$outdir/${r_sid}_trimmed.fastq" \
       		--detect_adapter_for_pe -5 -3 -D --dup_calc_accuracy 3 -c -p -P 20 \
            -h "$outdir/${base_sid}.html" -w 14 > "$logdir"/log.txt) & spinner $!
    fi
fi

