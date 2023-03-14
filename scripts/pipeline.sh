PIPEVERSION="1.0 - sra_automatization_pipeline"
STARTTIME=`date +'%y-%m-%d %H:%M:%S'`
RUN_ID=`date +"%Y%m%d%H%M%S"`

RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

# This script's design guidelines are based on
# the script available at https://github.com/giuliospinozzi/arpir.

echo "
  +--------------------------------------------------------+
  |                                                        |
  |        Illumina Pipeline for SRA automatization        |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Miguel Ramon                                |
  |  Date:     March 2022                                  |
  |  Contact:  miralnso@proton.me	                   |
  |  Version:  1.0 - SRA - RNAseq alignment                |
  +--------------------------------------------------------+
"
printf "${GREEN}\n#### RNA-SEQ PIPELINE for FINAL MODULE PROJECT ####${NC}\n\n"
printf "${YELLOW}Pipeline started at $(date +'%H:%M:%S')
___________________________________________________________ ${NC}\n"

# Stop execution when having a non-zero status and trap errors giving line number
set -e
trap 'printf ${RED}Error at line $LINENO${NC}' ERR

# Download samples and references if not exist
bash scripts/download.sh

# Different pipelines for different workflows
printf "${GREEEN}\nWORKFLOW PIPELINE MENU\n${NC}"

printf "${YELLOW}BASIC WORKFLOWS FOR HIGH QUALITY READS, NO PRE-PROCESSING\n${NC}"
echo -e "
\tWORKFLOW 1 (Using Salmon pseudo-aligner for basic error-correction capacities)\n
\tWORKFLOW 2 (Using STAR aligner)\n
\tWORKFLOW 3 (Using HISAT2 aligner)\n"
printf "${YELLOW}ADVANCED WORKFLOWS INCLUDING PREPROCESSING\n${NC}"
echo -e "
\tWORKFLOW 4 (Toolset: FastQScreen, Cutadapt, STAR)\n
\tWORKFLOW 5 (Toolset: FastQScreen, Cutadapt, HISAT2)\n
\tWORKFLOW 6 (Toolset: FastQScreen, Cutadapt, Salmon)\n
\tWORKFLOW 7 (Toolset: FastQScreen, Cutadapt, Kallisto)\n
\t(Trimmomatic workflows not working as of now)\n
\tWORKFLOW 8 (Toolset: FastQScreen, Trimmomatic, STAR)\n
\tWORKFLOW 9 (Toolset: FastQScreen, Trimmomatic, HISAT2)\n
\tWORKFLOW 10 (Toolset: FastQScreen, Trimmomatic, Salmon)\n
\tWORKFLOW 11 (Toolset: FastQScreen, Trimmomatic, Kallisto)\n"

read -rp "Option: " menuOp
# de repente me ha cambiado la dirección de analisis por la cara!
# for samples found in dumped_fastq dir
for f_path in $(find res/samples/dumped_fastq -mindepth 2 -type f -name "*_1.fastq"); do \
    r_path=${f_path/_1/_2} \
    f_name=$(basename "$f_path") \
    r_name=$(basename "$r_path") \
    sid=${f_name%_*}  # extract the SRA entry ID from the filename
    echo -e "\nForward file for $sid: $f_path \n
Reverse file for $sid: $r_path"

	case $menuOp in
		# for basic workflows
		1 )
 			mkdir -p out/aligned/star/star log/aligned/star/star
			outdir="aligned/star/star"
			bash scripts/index.sh "STAR"
			bash scripts/align.sh "STAR" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"
		;;
		2 )
			mkdir -p out/aligned/hisat2/hisat2 log/aligned/hisat2/hisat2
			outdir="aligned/hisat2/hisat2"
			bash scripts/index.sh "HISAT2"			
			bash scripts/align.sh "HISAT2" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"
		;;
		3 )
			mkdir -p out/aligned/salmon/salmon log/aligned/salmon/salmon
			outdir="aligned/salmon/salmon"
			bash scripts/index.sh "SALMON"
			bash scripts/align.sh "SALMON" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"			
		;;

		# for cutadapt workflows
		4 | 5 | 6 | 7 )
			
			mkdir -p out/trimmed/cutadapt/$sid log/trimmed/cutadapt/$sid
			trimdir="trimmed/cutadapt/$sid"
			bash scripts/pre_proc.sh "cutadapt" "$f_path" "$r_path" "out/$trimdir" "log/$trimdir"
			
			if [ "$menuOp" == "4" ]; then
				
				bash scripts/index.sh "STAR" 
				
				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/star/star_cutadapt log/star/star_cutadapt
					outdir="star/star_cutadapt"
					bash scripts/align.sh "STAR" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "cutadapt"
				done

			elif [ "$menuOp" == "5" ]; then
				
				bash scripts/index.sh "HISAT2"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/hisat2/hisat2_cutadapt log/hisat2/hisat2_cutadapt
					outdir="hisat2/hisat2_cutadapt"
					bash scripts/align.sh "HISAT2" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "cutadapt"
				done

			elif [ "$menuOp" == "6" ]; then

				bash scripts/index.sh "SALMON"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/salmon/salmon_cutadapt log/salmon/salmon_cutadapt
					outdir="salmon/salmon_cutadapt"
					bash scripts/align.sh "SALMON" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "cutadapt"
				done

			elif [ "$menuOp" == "7" ]; then

				bash scripts/index.sh "KALLISTO"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/kallisto/kallisto_cutadapt log/kallisto/kallisto_cutadapt
					outdir="kallisto/kallisto_cutadapt"
					bash scripts/align.sh "KALLISTO" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "cutadapt"
				done
			fi
		;;
		# for trimmomatic workflows (might remove them or comment them, let's see)
		8 | 9 | 10 | 11 )
			
			mkdir -p out/trimmed/trimmomatic/$sid log/trimmed/trimmomatic/$sid
			trimdir="trimmed/trimmomatic/$sid"
			bash scripts/pre_proc.sh "trimmomatic" "$f_path" "$r_path" "out/$trimdir" "log/$trimdir"
			
			if [ "$menuOp" == "8" ]; then

				bash scripts/index.sh "STAR"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/star/star_trimmomatic log/star/star_trimmomatic
					outdir="star/star_trimmomatic"
					bash scripts/align.sh "STAR" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "trimmomatic"
				done

			elif [ "$menuOp" == "9" ]; then

				bash scripts/index.sh "HISAT2"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/hisat2/hisat2_trimmomatic log/hisat2/hisat2_trimmomatic
					outdir="hisat2/hisat2_trimmomatic"
					bash scripts/align.sh  "HISAT2" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "trimmomatic"
				done

			elif [ "$menuOp" == "10" ]; then

				bash scripts/index.sh "SALMON"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p out/salmon/salmon_trimmomatic log/salmon/salmon_trimmomatic
					outdir="salmon/salmon_trimmomatic"
					bash scripts/align.sh "SALMON" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "trimmomatic"
				done

			elif [ "$menuOp" == "11" ]; then

				bash scripts/index.sh "KALLISTO"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do				
					mkdir -p out/kallisto/kallisto_trimmomatic log/kallisto/kallisto_trimmomatic
					outdir="kallisto/kallisto_trimmomatic"
					bash scripts/align.sh "KALLISTO" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "trimmomatic"
				done
			fi
		;;
	esac
	
	# Postprocessing with SAMtools, htseq and deeptools
	echo -e "\nPerforming post-alignment steps...\n"
	bash post_proc.sh $sid 

done

 
echo -e "\n\n############ Pipeline finished at $(date +'%H:%M:%S') ##############\n"
