#### DOWNLOAD SCRIPT ####

source "scripts/spinner.sh"
RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

mkdir -p "res/samples"	
mkdir -p "data/assembly/reference_grch38"
sample_dir="res/samples"
genome_dir="data/assembly/reference_grch38"

# Activate if required or if willing to use whole genome
# echo -e "\Downloading GRCh38 genome, transcriptome and annotations files\n"
# wget -nc -P $genome_dir --content-disposition "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz"
# wget -nc -P $genome_dir --content-disposition "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
# wget -nc -P $genome_dir --content-disposition "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"

echo -e "\nEnter one or more SRA accessions (e.g. SRR123456) separated by spaces, or enter 'S' to stop: \n"

#SRAentry=0


while true; do
    read -a SRAentries  # Read a line of input and split it into an array of accessions
    if [[ "${SRAentries[0]}" =~ ^[SsNn0]$ ]]; then
        # If the first element is 'S', 's', 'N', 'n' or '0', exit the loop
        exit 1
		
    else
        for SRAentry in "${SRAentries[@]}"; do
            if [[ "$SRAentry" =~ ^[Ss]RR[0-9]{6}$ ]]; then
                echo "Matched SRA accession $SRAentry"
                (wget -nc -P $sample_dir/SRA --content-disposition https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRAentry/$SRAentry) & spinner
                echo "Finished downloading $SRAentry"
				
				# Dump fastq from SRA entry
				mkdir -p $sample_dir/dumped_fastq
				echo -e "\nDumping fastq files from $SRAentry...\n"
					if [ "$(ls -A $sample_dir/dumped_fastq/$SRAentry)" ]; then
						echo -e "Fastq already dumped for $SRAentry, skipping.\n"
					
					else
						(fasterq-dump -fp -O $sample_dir/dumped_fastq/$SRAentry $SRAentry) & spinner $!
					fi				
            
			else
                echo "Did not match SRA accession $SRAentry. Skipping download."
            fi
        done
		break
    fi
done


# QC assay

read -rp "Would you like to perform a QC analysis (FastQC and FastQScreen)? (Y/n): " performQC

for SRAentry in "${SRAentries[@]}"; do
	for sid in $(find $sample_dir/dumped_fastq/$SRAentry -type f -name '*.fastq'); do
		base_sid=$(basename $sid | cut -d"_" -f1)

		if [ "$(ls -A $sample_dir/dumped_fastq/$base_sid)" ]; then
			case $performQC in
				[Yy]* )
					bash scripts/qc.sh $sid "out/qc/fastqc" "out/qc/fastq_screen"
				;;
				[Nn]* )
					echo -e "\nSkipping QC analysis...\n"
					break 
				;;
			esac
		fi
	done
done


# Results visualization
case $performQC in
[Yy]* )
	printf "${GREEN}\nNow, take your time to give a look to the QC analysis.\n${NC}"
	echo "When you are ready, press any key to continue..."
	read -n 1 -s -r -p ""
;;
esac