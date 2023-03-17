#### POST-PROC SCRIPT ####

source "scripts/spinner.sh"
RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

ref_gen="data/assembly/reference_grch38/Homo_sapiens.GRCh38.dna.chromosome.21.fa"
ref_cdna="data/assembly/reference_grch38/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ref_gtf="data/assembly/reference_grch38/Homo_sapiens.GRCh38.109.chr21.gtf"

tool=$1
bam=$2
countdir=$3

# Comprobation for independent use of script
if [ "$#" -ne 3 ]
then
    printf "${RED}Usage: $1 <tool> $2 <bam_file> $3 <counts_out_dir>${NC}\n"
    echo -e 'tool: "featurecounts", "htseq"\n'
	exit 1
fi

# Counts using FEATURECOUNTS

if [ "$1" == "featurecounts" ]; then

    (featureCounts -T 14 -s 2 -p --countReadPairs -a "$ref_gtf" -t exon -g gene_id \
    -o "$countdir/counts.txt" "$bam") & spinner $!

# Counts using htseq

elif [ "$1" == "htseq" ]; then

    base_bam=$(basename $bam .bam)

    (htseq-count --format=bam --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id    \
    --additional-attr=gene_name "$bam" "$ref_gtf" > "$countdir/$base_bam.htseq") & spinner $!

fi
