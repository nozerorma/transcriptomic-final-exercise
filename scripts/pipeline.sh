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

echo -e "\nPerforming sample"
