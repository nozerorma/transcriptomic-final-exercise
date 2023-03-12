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

	case $menuOp in
		# for basic workflows
		1 )
			mkdir -p out/aligned/salmon/salmon log/aligned/salmon/salmon
			outdir="aligned/salmon/salmon"
			bash scripts/align.sh salmon $sid out/$outdir log/$outdir 
		;;
		2 )
			mkdir -p out/aligned/star/star log/aligned/star/star
			outdir="aligned/star/star" 
			bash scripts/align.sh star $sid out/$outdir log/$outdir 
		;;
		3 )
			mkdir -p out/aligned/hisat2/hisat2 log/aligned/hisat2/hisat2
			outdir="aligned/hisat2/hisat2"
			bash scripts/align.sh hisat $sid out/$outdir log/$outdir
		;;

		# for cutadapt workflows
		4 | 6 | 8 | 10 )
			mkdir -p out/trimmed/cutadapt/$base_sid log/trimmed/cutadapt/$base_sid
			trimdir="trimmed/cutadapt/$base_sid"
			bash scripts/pre_proc.sh fastq-screen cutadapt $sid out/$outdir log/$outdir
			case $menuOp in
				4 )
					for trimmed_sid in $(find $trimdir -type f -name \*); do
						mkdir -p out/star/star_cutadapt log/star/star_cutadapt
						outdir="star/star_cutadapt"
						bash scripts/align.sh star $sid out/$outdir log/$outdir
					done
				;;
				6 )
					for trimmed_sid in $(find $trimdir -type f -name \*); do
						mkdir -p out/hisat2/hisat2_cutadapt log/hisat2/hisat2_cutadapt
						outdir="hisat2/hisat2_cutadapt"
						bash scripts/align.sh hisat $sid out/$outdir log/$outdir
					done
				;;
				8 )
					for trimmed_sid in $(find $trimdir -type f -name \*); do
						mkdir -p out/salmon/salmon_cutadapt log/salmon/salmon_cutadapt
						outdir="salmon/salmon_cutadapt"
						bash scripts/align.sh salmon $sid out/$outdir log/$outdir
					done
				;;
				10 )
					for trimmed_sid in $(find $trimdir -type f -name \*); do
						mkdir -p out/kallisto/kallisto_cutadapt log/kallisto/kallisto_cutadapt
						outdir="kallisto/kallisto_cutadapt"
						bash scripts/align.sh kallisto $sid out/$outdir log/$outdir
					done
				;;
		;;

		# for trimmomatic workflows
		5 | 7 | 9 | 11 )
			mkdir -p out/trimmed/trimmomatic/$base_sid log/trimmed/trimmomatic/$base_sid
			trimdir="trimmed/trimmomatic/$base_sid"
			bash scripts/pre_proc.sh fastq-screen trimmomatic $sid out/$outdir log/$outdir
			case $menuOp in	
			5 )
				for trimmed_sid in $(find $trimdir -type f -name \*); do
					mkdir -p out/star/star_trimmomatic log/star/star_trimmomatic
					outdir="star/star_trimmomatic"
					bash scripts/align.sh star $sid out/$outdir log/$outdir
				done
			;;
			7 )
				for trimmed_sid in $(find $trimdir -type f -name \*); do
					mkdir -p out/hisat2/hisat2_trimmomatic log/hisat2/hisat2_trimmomatic
					outdir="hisat2/hisat2_trimmomatic"
					bash scripts/align.sh  hisat $sid out/$outdir log/$outdir
				done
			;;
			9 )
				for trimmed_sid in $(find $trimdir -type f -name \*); do
					mkdir -p out/salmon/salmon_trimmomatic log/salmon/salmon_trimmomatic
					outdir="salmon/salmon_trimmomatic"
					bash scripts/align.sh salmon $sid out/$outdir log/$outdir
				done
			;;
			11 )
				for trimmed_sid in $(find $trimdir -type f -name \*); do				
					mkdir -p out/kallisto/kallisto_trimmomatic log/kallisto/kallisto_trimmomatic
					outdir="kallisto/kallisto_trimmomatic"
					bash scripts/align.sh kallisto $sid out/$outdir log/$outdir
				done
			;;
		;;
	esac

	# Postprocessing with SAMtools, htseq and deeptools
	echo -e "\nPerforming post-alignment steps...\n"
	bash post_proc.sh $sid 
done

 
echo -e "\n\n############ Pipeline finished at $(date +'%H:%M:%S') ##############\n"