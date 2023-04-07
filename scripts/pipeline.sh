PIPEVERSION="1.0 - sra_automatization_pipeline"
STARTTIME=`date +'%y-%m-%d %H:%M:%S'`
RUN_ID=`date +"%Y%m%d%H%M%S"`

source "scripts/spinner.sh"
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
echo -e "${GREEN}\n#### RNA-SEQ PIPELINE for FINAL MODULE PROJECT ####${NC}\n
${YELLOW}Pipeline started at $(date +'%H:%M:%S')
___________________________________________________________ ${NC}\n"

# Stop execution when having a non-zero status and trap errors giving line number
#set -e
trap 'printf "${RED}Error at line $LINENO${NC}"' ERR

# Check if required conda environment is present, else create

if { conda env list | grep 'SRA_pipeline'; } >/dev/null 2>&1; then
	eval "$(conda shell.bash hook)" >/dev/null 2>&1
	conda activate SRA_pipeline >/dev/null 2>&1
	echo "Running conda environment: $CONDA_DEFAULT_ENV"

else
	if [ -d ~/mambaforge ] || [ -d ~/miniforge ]; then
		(echo "Creating conda environment: SRA_pipeline"
		mamba env create -f envs/SRA_pipeline.yaml >/dev/null 2>&1) & spinner $!
		eval "$(conda shell.bash hook)"
		conda activate SRA_pipeline >/dev/null 2>&1
		echo "Running conda environment: $CONDA_DEFAULT_ENV"
	elif [ -d ~/condaforge/ ] || [ -d ~/miniconda ] || [ -d /anaconda3 ]; then
		(echo "Creating conda environment: SRA_pipeline"
		conda env create -f envs/SRA_pipeline.yaml >/dev/null 2>&1) & spinner $!
		eval "$(conda shell.bash hook)"
		conda activate SRA_pipeline >/dev/null 2>&1
		echo "Running conda environment: $CONDA_DEFAULT_ENV"
	fi
fi


# Download samples and references if not exist
bash scripts/download.sh
#retn_code=$?
if [ $? == "1" ]; then
	echo -e "\n\n############ Pipeline finished at $(date +'%H:%M:%S') ##############\n"
	exit 0
fi

# Different pipelines for different workflows
echo -e "${GREEN}\nWORKFLOW PIPELINE MENU\n${NC}
${RED}Default cores set to 14, change accordingly\n${NC}"

echo -e "
${YELLOW}BASIC WORKFLOWS FOR HIGH QUALITY READS, NO PRE-PROCESSING\n${NC}
\tWORKFLOW 1 (Using STAR aligner)\n
\tWORKFLOW 2 (Using HISAT2 aligner)\n
\tWORKFLOW 3 (Using SALMON pseudo-aligner)\n
\tWORKFLOW 4 (Using KALLISTO pseudo-aligner)\n
${YELLOW}ADVANCED WORKFLOWS INCLUDING PREPROCESSING\n${NC}
\tWORKFLOW 5 (Toolset: FastQScreen, Cutadapt, STAR)\n
\tWORKFLOW 6 (Toolset: FastQScreen, Cutadapt, HISAT2)\n
\tWORKFLOW 7 (Toolset: FastQScreen, Cutadapt, SALMON)\n
\tWORKFLOW 8 (Toolset: FastQScreen, Cutadapt, KALLISTO)\n
\tWORKFLOW 9 (Toolset: FastQScreen, Fastp, STAR)\n
\tWORKFLOW 10 (Toolset: FastQScreen, Fastp, HISAT2)\n
\tWORKFLOW 11 (Toolset: FastQScreen, Fastp, SALMON)\n
\tWORKFLOW 12 (Toolset: FastQScreen, Fastp, KALLISTO)\n"

read -rp "Option: " menuOp
# de repente me ha cambiado la dirección de analisis por la cara!
# for samples found in dumped_fastq dir
for f_path in $(find res/samples/dumped_fastq -mindepth 2 -type f -name "*_1.fastq"); do \
    r_path=${f_path/_1/_2} \
    sid=$(basename "$f_path" _1.fastq)   # extract the SRA entry ID from the filename
    echo -e "\nForward file for $sid: $f_path \n
Reverse file for $sid: $r_path"

	case $menuOp in
		# for basic workflows
		1 )
 			mkdir -p "out/aligned/star/$sid/$RUN_ID" "log/aligned/star/$sid/$RUN_ID"
			outdir="aligned/star/$sid/$RUN_ID"
			bash scripts/index.sh "STAR"
			bash scripts/align.sh "STAR" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"
		;;
		2 )
			mkdir -p "out/aligned/hisat2/$sid/$RUN_ID" "log/aligned/hisat2/$sid/$RUN_ID"
			outdir="aligned/hisat2/$sid/$RUN_ID"
			bash scripts/index.sh "HISAT2"			
			bash scripts/align.sh "HISAT2" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"
		;;
		3 )
			mkdir -p "out/aligned/salmon/$sid/$RUN_ID" "log/aligned/salmon/$sid/$RUN_ID"
			outdir="aligned/salmon/$sid/$RUN_ID"
			bash scripts/index.sh "SALMON"
			bash scripts/align.sh "SALMON" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"			
		;;
		4 )
			mkdir -p "out/aligned/kallisto/$sid/$RUN_ID" "log/aligned/kallisto/$sid/$RUN_ID"
			outdir="aligned/kallisto/$sid/$RUN_ID"
			bash scripts/index.sh "KALLISTO"
			bash scripts/align.sh "KALLISTO" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "untrimmed"		
		;;

		# for cutadapt workflows
		5 | 6 | 7 | 8 )
			
			mkdir -p out/trimmed/cutadapt/$sid log/trimmed/cutadapt/$sid
			trimdir="trimmed/cutadapt/$sid"
			bash scripts/pre_proc.sh "cutadapt" "$f_path" "$r_path" "out/$trimdir" "log/$trimdir"

			for f_trimmed in $(find "out/$trimdir" -type f -name "*_1_trimmed.fastq"); do
				r_trimmed=${f_trimmed/_1/_2}
				sid=$(basename "$f_trimmed" _1_trimmed.fastq) 
				if [ "$menuOp" == "5" ]; then
				
					bash scripts/index.sh "STAR" 
					mkdir -p "out/aligned/star_cutadapt/$sid/$RUN_ID" 
					mkdir -p "log/aligned/star_cutadapt/$sid/$RUN_ID"
					outdir="aligned/star_cutadapt/$sid/$RUN_ID"
					bash scripts/align.sh "STAR" "$f_trimmed" "$r_trimmed" "out/$outdir" "log/$outdir" "cutadapt"
				
				elif [ "$menuOp" == "6" ]; then
					
					bash scripts/index.sh "HISAT2" 

					mkdir -p "out/aligned/hisat2_cutadapt/$sid/$RUN_ID" 
					mkdir -p "log/aligned/hisat2_cutadapt/$sid/$RUN_ID"
					outdir="aligned/hisat2_cutadapt/$sid/$RUN_ID"
					bash scripts/align.sh "HISAT2" "$f_trimmed" "$r_trimmed" "out/$outdir" "log/$outdir" "cutadapt"

				elif [ "$menuOp" == "7" ]; then

					bash scripts/index.sh "SALMON"

					mkdir -p "out/aligned/salmon_cutadapt/$sid/$RUN_ID" 
					mkdir -p "log/aligned/salmon_cutadapt/$sid/$RUN_ID"
					outdir="aligned/salmon_cutadapt/$sid/$RUN_ID"
					bash scripts/align.sh "SALMON" "$f_trimmed" "$r_trimmed" "out/$outdir" "log/$outdir" "cutadapt"

				elif [ "$menuOp" == "8" ]; then

					bash scripts/index.sh "KALLISTO"

					mkdir -p "out/aligned/kallisto_cutadapt/$sid/$RUN_ID" 
					mkdir -p "log/aligned/kallisto_cutadapt/$sid/$RUN_ID"
					outdir="aligned/kallisto_cutadapt/$sid/$RUN_ID"
					bash scripts/align.sh "KALLISTO" "$f_trimmed" "$r_trimmed" "out/$outdir" "log/$outdir" "cutadapt"
				
				fi
			done
		;;

		9 | 10 | 11 | 12 )
			
			mkdir -p out/trimmed/fastp/$sid log/trimmed/fastp/$sid
			trimdir="trimmed/fastp/$sid"
			bash scripts/pre_proc.sh "fastp" "$f_path" "$r_path" "out/$trimdir" "log/$trimdir"
			
			if [ "$menuOp" == "9" ]; then

				bash scripts/index.sh "STAR"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir "out/aligned/star_fastp/$sid/$RUN_ID" "log/aligned/star_fastp/$sid/$RUN_ID"
					outdir="aligned/star_fastp/$sid/$RUN_ID"
					bash scripts/align.sh "STAR" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "fastp"
				done

			elif [ "$menuOp" == "10" ]; then

				bash scripts/index.sh "HISAT2"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p "out/aligned/hisat2_fastp/$sid/$RUN_ID" "log/aligned/hisat2_fastp/$sid/$RUN_ID"
					outdir="aligned/hisat2_fastp/$sid/$RUN_ID"
					bash scripts/align.sh  "HISAT2" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "fastp"
				done

			elif [ "$menuOp" == "11" ]; then

				bash scripts/index.sh "SALMON"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do
					mkdir -p "out/aligned/salmon_fastp/$sid/$RUN_ID" "log/aligned/salmon_fastp/$sid/$RUN_ID"
					outdir="aligned/salmon_fastp/$sid/$RUN_ID"
					bash scripts/align.sh "SALMON" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "fastp"
				done

			elif [ "$menuOp" == "12" ]; then

				bash scripts/index.sh "KALLISTO"

				for trimmed_sid in $(find "out/$trimdir" -type f -name \*); do				
					mkdir -p "out/aligned/kallisto_fastp/$sid/$RUN_ID" "log/aligned/kallisto_fastp/$sid/$RUN_ID"
					outdir="aligned/kallisto_fastp/$sid/$RUN_ID"
					bash scripts/align.sh "KALLISTO" "$f_path" "$r_path" "out/$outdir" "log/$outdir" "fastp"
				done
			fi
		;;
	esac

done

 
echo -e "\n\n############ Pipeline finished at $(date +'%H:%M:%S') ##############\n"
