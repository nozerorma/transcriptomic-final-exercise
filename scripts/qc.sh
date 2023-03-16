#### QC script ####

source $(dirname "$0")/spinner.sh
RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

mkdir -p "out/qc/fastqc"
mkdir -p "out/qc/fastq_screen"

sid=$2
base_sid=$(basename "$sid")
nofastq_sid=$(basename "$sid" .fastq)
cut_sid=$(basename "$sid" | cut -d"_" -f1)
fastqc_dir=$3
fastqscreen_dir=$4

# FASTQC
if [ "$1" == "run" ]; then
	
	echo -e "${YELLOW}\nPerforming QC analysis for $nofastq_sid...
___________________________________________________________ ${NC}\n"

	if [ -f $fastqc_dir/$cut_sid/${nofastq_sid}_fastqc.html ]; then
		echo -e "FastQC analysis already performed for $nofastq_sid, skipping analysis.\n" 

	else
		echo -e "Running FastQC analysis..."
		mkdir -p $fastqc_dir/$cut_sid
		mkdir -p log/qc
		(fastqc -o $fastqc_dir/$cut_sid $sid 2>&1 >/dev/null | tail -n +2 2>&1 log/qc/fastqc.log) & spinner $!
		echo
	fi

	# FASTQSCREEN
	if [ -f $fastqscreen_dir/$cut_sid/${nofastq_sid}_screen.html ]; then
		echo -e "FastQScreen analysis already performed for $nofastq_sid, skipping analysis.\n" 
		
	else
		screen_gen="res/fastq_screen_samples/FastQ_Screen_Genomes"
		if [ ! "$(ls -A $screen_gen)" ]; then
			echo "No reference genomes could be found."
			echo -e "\nWould you like to download genome indexes from database?"
			read -rp "It takes long AND it takes space... (Y/n): " genDownload
			
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
			echo -e "\nRunning FastQScreen analysis..."
			mkdir -p $fastqscreen_dir/$cut_sid
			mkdir -p log/qc
			(fastq_screen --conf "$screen_gen/fastq_screen.conf" \
				--tag --aligner bowtie2 --subset 100000 --threads 14 \
				--outdir "$fastqscreen_dir/$cut_sid" $sid 2>&1 log/qc/fastqscreen.log) & spinner $!
			echo
		fi
	fi
fi