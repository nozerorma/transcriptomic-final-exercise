#### ALIGNMENT SCRIPT ####

source "scripts/spinner.sh"
RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'


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

# Comprobation for independent use of script
if [ "$#" -ne 6 ]
then
    printf "${RED}Usage: $1 <tool> $2 <f_path> $3 <r_path> $4 <outdir> $5 <logdir> $6 <workflow>${NC}\n"
    echo -e 'tool: "STAR", "HISAT2", "SALMON", "KALLISTO"
workflow: "untrimmed", "cutadapt", "fastp"\n'
	exit 1
fi

# First lets see what workflow we are working with

sample_dir=0

if [ "$workflow" == "untrimmed" ]; then
	sample_dir="$untrimmed_dir"

elif [ "$workflow" == "cutadapt" ] || [ "$workflow" == "fastp" ] ; then
	sample_dir="$trimmed_dir"
fi

echo -e "${YELLOW}\nAligning $base_sid to reference with $tool...
___________________________________________________________ ${NC}\n"


# overview from https://www.reneshbedre.com/blog/star-aligner.html#mapping-reads-to-genome
# STAR
if [ "$1" == "STAR" ]; then

	(STAR --runThreadN 14 --readFilesIn "$f_path" $r_path  \
			--genomeDir "$index_dir" --outReadsUnmapped Fastx  \
			--outFileNamePrefix "$outdir/" \
			--outSAMtype BAM SortedByCoordinate ) & spinner $!
		
		# add more params for statistics (there's a few problems here to solve)
		bam_file=$(find "$outdir" -type f -name "*.bam")
		(samtools stats "$bam_file" > "$bam_file".txt) & spinner $!
		(samtools index "$bam_file") & spinner $!
		(bamCoverage -b "$bam_file" -o "$outdir/$base_sid.bw" --normalizeUsing BPM) & spinner $!
			# index bam here with samtools

		echo -e "\nSample $base_sid aligned using $tool.\n"

	# POST_PROCESSING: READ COUNT
	read -rp "Which tool would you like to use for feature count? (featurecounts/htseq) " counts
	mkdir "$outdir/counts"
	countdir="$outdir/counts"
	bash scripts/post_proc.sh "$counts" "$bam_file" "$countdir"
	

# overview from https://bioinfo-dirty-jobs.github.io/rana2//lectures/07.rnaseq_hisat2/
# HISAT2
elif [ "$1" == "HISAT2" ]; then

	(hisat2 --new-summary --summary-file "$outdir/$base_sid.hisat2.summary" \
	-p 14 -x "$index_dir/HISAT2" -1 "$f_path" -2 "$r_path" -k 1 -S "$outdir/$base_sid.sam") & spinner $!
	
	# problema con el pipe
	# add more params for statistics
	(samtools view -bS "$outdir/$base_sid.sam" > "$outdir/$base_sid.bam") & spinner $!
	(samtools stats "$outdir/$base_sid.bam" > "$outdir/$base_sid.txt") & spinner $!
	(samtools sort --write-index "$outdir/$base_sid.bam" -o "$outdir/$base_sid.sorted.bam") & spinner $!
	(samtools stats "$outdir/$base_sid.sorted.bam" > "$outdir/$base_sid.sorted.txt") & spinner $!
	(bamCoverage -b "$outdir/$base_sid.sorted.bam" -o "$outdir/$base_sid.bw" --normalizeUsing BPM) & spinner $!

	# was going to run picard but it has a bunch of incompatibilities
	# with the tools I'm already using, aborting

	# POST_PROCESSING: READ COUNT
	read -rp "Which tool would you like to use for feature count? (featurecounts/htseq) " counts
	mkdir "$outdir/counts/"
	countdir="$outdir/counts"
	bash scripts/post_proc.sh "$counts" "$outdir/$base_sid.sorted.bam" "$countdir"

## ALIGNMENT
### PSEUDO-TOOLS

# overview from https://salmon.readthedocs.io/en/latest/salmon.html
# SALMON (mapping-based mode, using GTF annotations)
# -l set as A for automatic guessing of strandness, change accordingly
elif [ "$1" == "SALMON" ]; then
		
	salmon quant -i "$index_dir" -l A -1 "$f_path" -2 "$r_path" --validateMappings \
		-o "$outdir" -g "$ref_gtf" -p 14 --writeMappings 	

# KALLISTO
# tema de los gtf en ambos, no se si sin meterlo en idx tiene sentido
elif [ "$1" == "KALLISTO" ]; then
	
	base_ref_cdna=$(basename "$ref_cdna" .gz)

	(kallisto quant -i "$index_dir/$base_ref_cdna.idx" --bias --fusion \
	"$f_path" "$r_path" -o "$outdir" --pseudobam --genomebam \
	--gtf "$ref_gtf" -t 14) & spinner $!

	pseudobam="$outdir/pseudoalignments.bam"
	(bamCoverage -b "$pseudobam" -o "$outdir/$base_sid.pseudoalignments.bw" \
	--normalizeUsing BPM 2>&1) & spinner $!	
	
	# POST_PROCESSING: READ COUNT
	read -rp "Which tool would you like to use for feature count? (featurecounts/htseq) " counts
	mkdir "$outdir/counts/"
	countdir="$outdir/counts"
	bash scripts/post_proc.sh "$counts" "$pseudobam" "$countdir"

fi