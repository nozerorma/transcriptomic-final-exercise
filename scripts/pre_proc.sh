#### PRE-PROCESSING SCRIPT ####

f_path=$3
r_path=$4
outdir=$5
logdir=$6

f_name=$(basename "$f_path")
r_name=$(basename "$r_path")
base_sid=$(basename "$f_path" | cut -d"_" -f1)
f_sid=$(basename "$f_name" .fastq)
r_sid=$(basename "$r_name" .fastq)

if [ "$1" == "fastq-screen" ]; then
    
    mkdir -p "log/qc/fastq_screen"
    mkdir -p "out/qc/fastq_screen"
    echo -e "\nRunning FastQScreen...\n"
    
    # a esto le tengo que dar la vuelta, redundante ("$(ls -A $outdir)")
    read -rp "Would you like to download genome indexes from database? \
It takes long AND it takes space... (Y/n): " genDownload
    case $genDownload in
        [Yy]* )
            if [ ! -d "res/fastq_screen_samples" ]; then
                echo -e "\nDownloading genomes...\n" 
                fastq_screen --get_genomes --outdir "res/fastq_screen_samples" > log/qc/fastq_screen/fastq_screen.log
            else
                echo -e "\nGenomes already present in folder, skipping...\n"
            fi
            ;;
        [Nn]* )
            echo -e "\nSkipping genome download...\n"
            ;;
    esac
    
    read -rp "Would you like to run FastQScreen analysis? (Y/n): " isFastqs
    case $isFastqs in
        [Yy]* )
            if [ ! -d "out/qc/fastq_screen/$base_sid" ]; then # probablemente debería hacerlo tb con -f
            fastq_screen --conf "res/fastq_screen_samples/FastQ_Screen_Genomes/fastq_screen.conf" \
            --tag --aligner bowtie2 --subset 100000 --threads 6 \
            --outdir "out/qc/fastq_screen/$base_sid" "$f_path" "$r_path" > log/qc/fastq_screen/fastq_screen.log
            fi
            ;;

        [Nn]* )
            echo -e "Skipping FastQScreen analysis...\n"
            ;;
    esac
fi

if [ "$2" == "cutadapt" ]; then
    
    echo -e "\nRunning cutadapt...\n" 
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Trim already performed for $f_sid, skipping...\n"
    else
        #took out -q 10 and --discard-untrimmed until further notice
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
            -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
            -a CTGTCTCTTATACACATCT...AGATGTGTATAAGAGACAG -a TGGAATTCTCGGGTGCCAAGG \
            -a "G{100}" -g "G{100}" -a "A{100}" -g "A{100}" -q 20 \
            -o "$outdir/${f_sid}_trimmed.fastq" -p "$outdir/${r_sid}_trimmed.fastq" \
            "$f_path" "$r_path" --cores 6 > "$logdir"/log.txt
    fi

# Uncomment and modify conveniently for using Trimmomatic
# elif [ "$2" == "trimmomatic" ]; then
    
#     trimmomatic PE -phred33 "$f_path" "$r_path" \
#             "$outdir"/"$f_name" "$outdir"/"$f_name"_unpaired \
#             "$outdir"/"$r_name" "$outdir"/"$r_name"_unpaired \
#             TRAILING:30 SLIDINGWINDOW:4:30 >> "$logdir"/log.txt

fi

# read -a "Would you like to re-run fastqc for your trimmed samples? (Y/n) " runfastqc