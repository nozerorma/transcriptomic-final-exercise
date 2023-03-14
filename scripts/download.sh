mkdir -p res/samples	
mkdir -p data/assembly/reference_grch38
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
        break
    else
        for SRAentry in "${SRAentries[@]}"; do
            if [[ "$SRAentry" =~ ^[Ss]RR[0-9]{6}$ ]]; then
                echo "Matched SRA accession $SRAentry"
                wget -nc -P $sample_dir/SRA --content-disposition https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRAentry/$SRAentry
                echo "Finished downloading $SRAentry"
            else
                echo "Did not match SRA accession $SRAentry. Skipping download."
            fi
        done
    fi
done

mkdir -p $sample_dir/dumped_fastq

for sra in $sample_dir/SRA/*; do
    SRAentry=$(basename "$sra")
    echo -e "\nDumping fastq files from $SRAentry...\n"
    if [ "$(ls -A $sample_dir/dumped_fastq/$SRAentry)" ]; then # pass this to other similar features
        echo -e "Fastq already dumped for $SRAentry, skipping.\n"
    else
        fasterq-dump -fp -O $sample_dir/dumped_fastq/$SRAentry $sample_dir/SRA/$SRAentry
    fi
done