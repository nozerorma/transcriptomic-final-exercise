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
    # a esto le tengo que dar la vuelta
    read -rp "Would you like to download genome indexes from database? \
It takes long AND it takes space... (Y/n): " genDownload
    case $genDownload in
        [Yy]* )
            if [ ! -d "res/fastq_screen_samples" ]; then
                echo -e "\nDownloading genomes...\n"
                #echo -e "\nFastQC log as of $(date +'%x                %H:%M:%S')" >> log/qc/fastq_screen/fastq_screen.log
                #echo -e "___________________________________________________________\n">> log/qc/fastq_screen/fastq_screen.log 
                fastq_screen --get_genomes --outdir "res/fastq_screen_samples" #>> log/qc/fastq_screen/fastq_screen.log
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
            if [ ! -d "out/qc/fastq_screen/$base_sid" ]; then
            fastq_screen --conf res/fastq_screen_samples/FastQ_Screen_Genomes/fastq_screen.conf \
            --aligner bowtie2 --subset 100000 --threads 6 \
            --outdir "out/qc/fastq_screen/$base_sid" "$f_path" "$r_path" #>> log/qc/fastq_screen/fastq_screen.log
            fi
            ;;

        [Nn]* )
            echo -e "Skipping FastQScreen analysis...\n"
            ;;
    esac
fi

echo -e"
if [ "$2" == "cutadapt" ]; then
    # mkdir -p $outdir/$f_name
    # mkdir -p $outdir/$r_name
    echo -e "\nRunning cutadapt...\n" 
    if [ -f "$outdir/$f_sid\_trimmed.fastq" ]; then
        echo -e "Trim already performed for $f_sid, skipping..."
    else
        cutadapt -o "$outdir"/"$f_sid"_trimmed.fastq -p \
        "$outdir"/"$r_sid"_trimmed.fastq \
        "$f_path" "$r_path" --cores 6 #>> "$logdir"/log.txt
    fi

elif [ "$2" == "trimmomatic" ]; then
    trimmomatic PE -phred33 "$f_path" "$r_path" \
    "$outdir"/"$f_name" "$outdir"/"$f_name"_unpaired \
    "$outdir"/"$r_name" "$outdir"/"$r_name"_unpaired \
    TRAILING:30 SLIDINGWINDOW:4:30 #>> "$logdir"/log.txt
fi"