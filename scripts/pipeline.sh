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
	if [ -d "$fastqc_dir/$base_sid" ]; then
		echo -e "FastQC already performed for $base_sid, skipping.\n"
		continue
	else
		mkdir -p $fastqc_dir/$base_sid
		echo "Started analysis of $sid"
		fastqc -o $fastqc_dir/$base_sid $sid
		continue
	fi
done
# Selectively show FastQC results
read -rp "Would you like to visualize your FastQC results? This is recommended for choosing the most appropriate pipeline for your workflow ('Y'/'N'): " visFastqc   
case $visFastqc in
	[Yy]* )
		# Show fastqc reports in firefox
       	echo -e "\nOpening firefox for FastQC report visualization...\n"
		for sid in $fastqc_dir; do
			for fastqc_rep in $(find $sid -type f -name *.html); do
				echo $fastqc_rep
				xdg-open $fastqc_rep 2>/dev/null
			done
		done;;
	
	[Nn]* )
		# Skip fastqc visualization
		echo -e "\nSkipping FastQC report visualization...\n";;	
esac

# Different pipelines for different workflows
echo -e "\nWORKFLOW PIPELINE MENU\n"
echo -e "
BASIC WORKFLOWS FOR HIGH QUALITY READS, NO PRE-PROCESSING\n
\tWORKFLOW 1 (Using Salmon pseudo-aligner for basic error-correction capacities)\n
\tWORKFLOW 2 (Using STAR aligner)\n
\tWORKFLOW 3 (Using HISAT2 aligner)\n
ADVANCED WORKFLOWS INCLUDING PREPROCESSING\n
\tWORKFLOW 4 (Toolset: FastQScreen, Cutadapt, STAR)\n
\tWORKFLOW 5 (Toolset: FastQScreen, Trimmomatic, STAR)\n
\tWORKFLOW 6 (Toolset: FastQScreen, Cutadapt, HISAT2)\n
\tWORKFLOW 7 (Toolset: FastQScreen, Trimmomatic, HISAT2)\n
\tWORKFLOW 8 (Toolset: FastQScreen, Cutadapt, Salmon)\n
\tWORKFLOW 9 (Toolset: FastQScreen, Trimmomatic, Salmon)\n
\tWORKFLOW 10 (Toolset: FastQScreen, Cutadapt, Kallisto)\n
\tWORKFLOW 11 (Toolset: FastQScreen, Trimmomatic, Kallisto)\n"

read -rp "Option: " menuOp

# for samples found in dumped_fastq dir
for sid in $(find res/samples/dumped_fastq -type f -name *.fastq | sort -u); do
	base_sid=$(basename $sid .fastq)
	NOW=$(date "+%Y-%m-%d")

	case $menuOp in
		1 )
			outdir="aligned/salmon/salmon-$NOW"
			bash scripts/align.sh salmon $sid out/$outdir log/$outdir 
		;;
		2 )
			outdir="aligned/star/star-$NOW" 
			bash scripts/align.sh star $sid out/$outdir log/$outdir 
		;;
		3 )
			outdir="aligned/hisat2/hisat2-$NOW"
			bash scripts/align.sh hisat $sid out/$outdir log/$outdir
		;;
		4 )
			outdir="trimmed/cutadapt/star/star-$NOW"
			bash scripts/align.sh star cutadapt $sid out/$outdir log/$outdir
		;;
		5 )
			outdir="trimmed/trimmomatic/star/star-$NOW"
			bash scripts/align.sh star trimmomatic $sid out/$outdir log/$outdir
		;;
		6 )
			outdir="trimmed/cutadapt/hisat2/hisat2-$NOW"
			bash scripts/align.sh hisat cutadapt $sid out/$outdir log/$outdir
		;;
		7 )
			outdir="trimmed/trimmomatic/hisat2/hisat2-$NOW"
			bash scripts/align.sh hisat trimmomatic $sid out/$outdir log/$outdir 
		;;
		8 )
			outdir="trimmed/cutadapt/salmon/salmon-$NOW"
			bash scripts/align.sh salmon cutadapt $sid out/$outdir log/$outdir
		;;
		9 )
			outdir="trimmed/trimmomatic/salmon/salmon-$NOW"
			bash scripts/align.sh salmon trimmomatic $sid out/$outdir log/$outdir
		;;
		10 )
			outdir="trimmed/trimmomatic/kallisto/kallisto-$NOW"
			bash scripts/align.sh kallisto cutadapt $sid out/$outdir log/$outdir
		;;
		11 )
			outdir="trimmed/trimmomatic/kallisto/kallisto-$NOW"
			bash scripts/align.sh kallisto trimmomatic $sid out/$outdir log/$outdir
		;;
	esac
done
echo -e "\n\n############ Pipeline finished at $(date +'%H:%M:%S') ##############\n"