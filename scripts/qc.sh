# QC script

# maybe insert logs?

mkdir -p "out/qc/fastqc"
mkdir -p "out/qc/fastqscreen"

sid=$2
base_sid=$(basename $sid)
fastqc_dir=$3
fastqscreen_dir=$4

echo -e "\nPerforming QC analysis...\n"


# FASTQC
if [ "$(ls -A $fastqc_dir)" ]; then
	echo -e "FastQC analysis already performed for $base_sid, skipping analysis.\n" 

else
	fastqc -o $fastqc_dir/$sid $sid

fi

# FASTQSCREEN
if  [ "$(ls -A $fastqscreen_dir)" ]; then
	echo -e "FastQScreen analysis already performed for $base_sid, skipping analysis.\n" 
	
else
	screen_gen="res/fastq_screen_samples/FastQ_Screen_Genomes/"
	if [ ! "$(ls -A $screen_gen)" ]; then
		echo "No reference genomes could be found."
		read -rp "Would you like to download genome indexes from database? \
		It takes long AND it takes space... (Y/n): " genDownload
		
		case $genDownload in
			[Yy]* )
				echo -e "\nDownloading genomes...\n" 
				fastq_screen --get_genomes --outdir "res/fastq_screen_samples"
			;;
			[Nn]* )
				echo -e "\nSkipping genome download...\n"
			;;
		esac
	
	else
		fastq_screen --conf "$screen_gen/fastq_screen.conf" \
            --tag --aligner bowtie2 --subset 100000 --threads 6 \
            --outdir "$fastqscreen_dir/$base_sid" $sid
    fi
fi

# Results visualization
if [ "$1" == "vis" ]; then
	read -rp "Visualize FastQC results in default browser (Y/n): " visFastqc
	case $visFastqc in
		[Yy]* )
			# Show fastqc reports in firefox
			echo -e "\nOpening firefox for QC report visualization...\n"
			
			for fastqc_rep in $(find out/qc -type f -name *.html); do
					xdg-open $fastqc_rep >/dev/null 2>&1
			done
			;;
		
		[Nn]* )
			# Skip fastqc visualization
			echo -e "\nSkipping FastQC report visualization...\n"
			;;	
	esac
fi
