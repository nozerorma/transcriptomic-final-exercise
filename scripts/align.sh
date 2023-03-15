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

# This section should be changed accordingly
# It would be nice of me to glob, but I don't seem to be able
ref_gen="data/assembly/reference_grch38/Homo_sapiens.GRCh38.dna.chromosome.21.fa"
ref_cdna="data/assembly/reference_grch38/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ref_gtf="data/assembly/reference_grch38/Homo_sapiens.GRCh38.109.chr21.gtf"

trimmed_dir="out/trimmed/$workflow"
untrimmed_dir="res/samples/dumped_fastq" # falta completar
index_dir="res/index/$tool"

# First lets see what workflow we are working with

sample_dir=0

if [ $workflow == "untrimmed" ]; then
	sample_dir=$untrimmed_dir

elif [ $workflow == "cutadapt" ] || [ $workflow == "trimmomatic" ] ; then
	sample_dir=$trimmed_dir
fi

echo -e "\nAligning $base_sid to reference with $tool...\n"

# overview from https://www.reneshbedre.com/blog/star-aligner.html#mapping-reads-to-genome
# STAR
if [ "$1" == "STAR" ]; then
	if [ "$(ls -A $outdir/$base_sid)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment.\n"   

	else
		STAR --runThreadN 6 --readFilesIn $f_path $r_path  \
			--genomeDir $index_dir --outReadsUnmapped Fastx  \
			--outFileNamePrefix $outdir/$base_sid \
			--outSAMtype BAM SortedByCoordinate 
		
		# add more params for statistics
		bam_file=$(find $outdir -type f -name "*.bam")
		samtools index $bam_file #> $logdir/$tool.log # - (necessary?)
			# index bam here with samtools

		echo -e "\nSample $base_sid aligned using $tool.\n"

	fi

# overview from https://bioinfo-dirty-jobs.github.io/rana2//lectures/07.rnaseq_hisat2/
# HISAT2
elif [ "$1" == "HISAT2" ]; then
	if [ "$(ls -A $outdir/$base_sid)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment.\n"   

	else
		hisat2 -p 6 --dta -x $index_dir -1 $f_path -2 $r_path \
		-S $outdir/$base_sid/$base_sid.sam 
		
		
		# add more params for statistics
		samtools view -bS $outdir/$base_sid/$base_sid.sam > $base_sid.bam | \
		samtools sort --write-index #> $logdir/$tool.log 
	
	fi

## ALIGNMENT
### PSEUDO-TOOLS

# overview from https://salmon.readthedocs.io/en/latest/salmon.html
# SALMON (mapping-based mode, using GTF annotations)
# -l set as A for automatic guessing of strandness, change accordingly
elif [ "$1" == "SALMON" ]; then
	if [ "$(ls -A $outdir/$base_sid)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment.\n" 
	else
		salmon quant -i $index_dir -l A -1 $f_path -2 $r_path --validateMappings \
			-o $outdir/$base_sid -g $ref_gtf -p 6 #> $logdir/$tool.log
	fi
# KALLISTO
# tema de los gtf en ambos, no se si sin meterlo en idx tiene sentido
elif [ "$1" == "KALLISTO" ]; then
	if [ "$(ls -A $outdir/$base_sid)" ]; then
		echo -e "Alignment already performed for $base_sid, skipping alignment.\n" 
	else
		base_ref_cdna=$(basename $ref_cdna .gz)
		kallisto quant -i $index_dir/$base_ref_cdna.idx --bias --fusion \
		$f_path $r_path -o $outdir/$base_sid --pseudobam --genomebam \
		--gtf $ref_gtf -t 6 
	fi
fi