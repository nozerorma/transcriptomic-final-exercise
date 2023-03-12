echo -e "\nEnter one or more SRA accessions (e.g. SRR123456) separated by spaces, or enter 'S' to stop: \n"

sample_dir="res/samples"
SRAentry=0

while true; do
    	read -a SRAentries  # Read a line of input and split it into an array of accessions
    	if [[ "${SRAentries[0]}" =~ ^[Ss]$ ]]; then
        	# If the first element is 'S', exit the loop
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
	break
	fi
done

echo -e "\nDumping fastq files from SRA...\n"
mkdir -p $sample_dir/dumped_fastq

for sra in $sample_dir/SRA/*; do
    	SRAentry=$(basename "$sra")
    	if [ -d $sample_dir/dumped_fastq/$SRAentry ]; then
        	echo -e "Fastq already dumped for $SRAentry, skipping.\n"
	else
		fasterq-dump -fp -O $sample_dir/dumped_fastq/$SRAentry $sample_dir/SRA/$SRAentry
	fi
done
