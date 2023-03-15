#### ALIGNMENT SCRIPT ####

tool=$1
f_path=$2
r_path=$3
outdir=$4
logdir=$5
workflow=$6

f_name=$(basename "$f_path")
r_name=$(basename "$r_path")
base_sid=$(basename "$f_path" | cut -d"_" -f1)
f_sid=$(basename "$f_name" .fastq)
r_sid=$(basename "$r_name" .fastq)

# Maybe change all the ls -A for -f, makes more sense
# This section should be changed accordingly
###### It would be nice of me to glob, but I don't seem to be able
ref_gen="data/assembly/reference_grch38/Homo_sapiens.GRCh38.dna.chromosome.21.fa"
ref_cdna="data/assembly/reference_grch38/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ref_gtf="data/assembly/reference_grch38/Homo_sapiens.GRCh38.109.chr21.gtf"

trimmed_dir="out/trimmed/$workflow"
untrimmed_dir="res/samples/dumped_fastq" # falta completar
index_dir="res/index/$tool"

# First lets see what workflow we are working with

sample_dir=0

if [ "$workflow" == "untrimmed" ]; then
	sample_dir="$untrimmed_dir"

elif [ "$workflow" == "cutadapt" ] || [ "$workflow" == "trimmomatic" ] ; then
	sample_dir="$trimmed_dir"
fi

echo -e "\nAligning $base_sid to reference with $tool...\n"


# overview from https://www.reneshbedre.com/blog/star-aligner.html#mapping-reads-to-genome
# STAR
if [ "$1" == "STAR" ]; then
	if [ "$(ls -A $outdir)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment...\n"
	
	else

		STAR --runThreadN 14 --readFilesIn "$f_path" $r_path  \
				--genomeDir "$index_dir" --outReadsUnmapped Fastx  \
				--outFileNamePrefix "$outdir" \
				--outSAMtype BAM SortedByCoordinate 
			
			# add more params for statistics
			bam_file=$(find "$outdir" -type f -name "*.bam")
			samtools stats "$bam_file" > "$bam_file".txt
			samtools index "$bam_file"
				# index bam here with samtools

			echo -e "\nSample $base_sid aligned using $tool.\n"
	fi
# overview from https://bioinfo-dirty-jobs.github.io/rana2//lectures/07.rnaseq_hisat2/
# HISAT2
elif [ "$1" == "HISAT2" ]; then

	if [ "$(ls -A $outdir)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment...\n"
	
	else
		
		hisat2 --new-summary --summary-file "$outdir.hisat2.summary" \
		-p 14 -x "$index_dir/HISAT2" -1 "$f_path" -2 "$r_path" -k 1 -S "$outdir.sam" 
		
		# problema con el pipe
		# add more params for statistics
		samtools view -bS "$outdir.sam" > "$outdir.bam"
		samtools stats "$outdir.bam" > "$outdir.txt"
		samtools sort --write-index "$outdir.bam" -o "$outdir.sorted.bam"
		samtools stats "$outdir.sorted.bam" > "$outdir.sorted.txt" 
		# was going to run picard but it has a bunch of incompatibilities
		# with the tools I'm already using, aborting
	fi

## ALIGNMENT
### PSEUDO-TOOLS

# overview from https://salmon.readthedocs.io/en/latest/salmon.html
# SALMON (mapping-based mode, using GTF annotations)
# -l set as A for automatic guessing of strandness, change accordingly
elif [ "$1" == "SALMON" ]; then
	if [ "$(ls -A $outdir)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment...\n"
	
	else
		
		salmon quant -i "$index_dir" -l A -1 "$f_path" -2 "$r_path" --validateMappings \
			-o "$outdir" -g "$ref_gtf" -p 14
	fi

# KALLISTO
# tema de los gtf en ambos, no se si sin meterlo en idx tiene sentido
elif [ "$1" == "KALLISTO" ]; then
	if [ "$(ls -A $outdir)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment...\n"
	
	else
		
		base_ref_cdna=$(basename "$ref_cdna" .gz)
		kallisto quant -i "$index_dir/$base_ref_cdna.idx" --bias --fusion \
		"$f_path" "$r_path" -o "$outdir" --pseudobam --genomebam \
		--gtf "$ref_gtf" -t 14 
	fi	
fi