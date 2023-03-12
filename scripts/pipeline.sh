echo -e "#### RNA-SEQ PIPELINE for FINAL MODULE PROJECT ####\n"
echo -e "Pipeline started at $(date +'%H:%M:%S')\n\n\n"

# Cleanup script
bash scripts/cleanup.sh

# Stop execution when having a non-zero status and trap errors giving line number
set -e
trap 'echo Error at about $LINENO' ERR

# Download reference genome GRCh38, GTF annotation and reference transcriptome
# not performed as they are included with the repo

# Download samples
mkdir -p res/samples	
bash scripts/download.sh

# Read quality control
echo -e "\nPerforming read QC...\n"
mkdir -p out/qc/fastqc
fastqc_dir="out/qc/fastqc"
for sid in $(find res/samples/dumped_fastq -type f -name *.fastq | sort -u); do
        base_sid=$(basename $sid .fastq | cut -d"_" -f1)
        mkdir -p $fastqc_dir/$base_sid
        if [ -d /$fastqc_dir/$base_sid ]; then
		echo -e "FastQC already performed for $base_sid, skipping.\n"
	else
		echo "Started analysis of $sid"
        	fastqc -o $fastqc_dir/$base_sid $sid
        fi
done
# Selectively show FastQC results
read -rp "Would you like to visualize your FastQC results? This is recommended for choosing the most appropriate pipeline for your workflow ('Y'/'N'): " visFastqc   
case $visFastqc in
	[Yy]* )
		# Show fastqc reports in firefox
       		echo -e "\nOpening firefox for FastQC report visualization...\n"
		for fastqc_rep in $(find $fastqc_dir/$base_sid -type f -name *.html | sort -u); do
			firefox fastqc_rep
		done
	break;;
	
	[Nn]* )
		# Skip fastqc visualization
		echo -e "\nSkipping FastQC report visualization...\n"

	break;;	
		
esac

# Different pipelines
while true; do
	echo -e "\nWorkflow pipeline menu.\n"
	echo -e "\nPlease choose an option according to your workflow, or enter 'E' to exit : \n
		\tPRE-PROCESSING TOOL\tS
		
