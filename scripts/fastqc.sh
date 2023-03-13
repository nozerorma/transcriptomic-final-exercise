# FASTQC script

# this can be much further optimized by inserting vars in calling scripts

# Read quality control
echo -e "\nPerforming read QC...\n"
mkdir -p out/qc/fastqc
fastqc_dir="out/qc/fastqc"
for sid in $(find res/samples/dumped_fastq -type f -name '*.fastq' | sort -u); do
	is_sid=$(basename "$sid" .fastq)
    base_sid=$(basename "$sid" .fastq | cut -d"_" -f1)
	
    if [ -f "$fastqc_dir/$base_sid/${is_sid}_fastqc.html" ]; then
		echo -e "FastQC already performed for "$is_sid", skipping.\n"
    
    else
		mkdir -p "$fastqc_dir"/"$base_sid"
		mkdir -p "log/qc/fastqc"
		echo "Started analysis of $sid"
		#echo -e "\nFastQC log as of $(date +'%x                %H:%M:%S')" >> log/qc/fastqc/$base_sid.log
		#echo -e "___________________________________________________________\n">> log/qc/fastqc/$base_sid.log 
		fastqc -o "$fastqc_dir"/"$base_sid" "$sid" >> log/qc/fastqc/$base_sid.log
	fi

done

# Selectively show FastQC results
read -rp "Would you like to visualize your FastQC results? This is \
recommended for choosing the most appropriate pipeline for your workflow (Y/n): " visFastqc   

case $visFastqc in
	[Yy]* )
		# Show fastqc reports in firefox
       	echo -e "\nOpening firefox for FastQC report visualization...\n"
		for sid in ${fastqc_dir}; do
			for fastqc_rep in $(find $sid -type f -name *.html); do
				echo $fastqc_rep
				xdg-open $fastqc_rep 2>/dev/null
			done
		done
		;;
	
	[Nn]* )
		# Skip fastqc visualization
		echo -e "\nSkipping FastQC report visualization...\n"
		;;	
esac