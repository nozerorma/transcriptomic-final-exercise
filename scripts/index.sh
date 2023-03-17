# INDEX script

tool=$1
#threads=$3

mkdir -p "log/index/$tool"
logdir="log/index/$tool"
mkdir -p "res/index/$tool"
outdir="res/index/$tool"

# This section should be changed accordingly
# It would be nice of me to glob, but I don't seem to be able
ref_gen="data/assembly/reference_grch38/Homo_sapiens.GRCh38.dna.chromosome.21.fa"
ref_cdna="data/assembly/reference_grch38/Homo_sapiens.GRCh38.cdna.all.fa.gz"

echo -e "${YELLOW}\nBuilding "$tool" index...
___________________________________________________________ ${NC}\n"

# Comprobation for independent use of script
if [ "$#" -ne 1 ]
then
    printf "${RED}Usage: $1 <tool> ${NC}\n"
    echo -e 'tool: "STAR", "HISAT2", "SALMON", "KALLISTO"\n'
	exit 1
fi

# STAR INDEX

if [ "$1" == "STAR" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        # change genomeSAindexNbases to 11 according to own program's advise, due to ref length
        STAR 	--runThreadN 14 --runMode genomeGenerate --genomeDir $outdir \
                --genomeFastaFiles $ref_gen --runRNGseed 1998 --genomeSAindexNbases 11 > $logdir/$tool.log

        echo -e "$tool index built.\n"
    fi

# HISAT2 INDEX
# From RNAseq hands-on session by Jaime Martínez de Villarreal, Epithelial Carcinogenesis Group, CNIO
elif [ "$1" == "HISAT2" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        hisat2-build -p 14 --seed 1998 $ref_gen $outdir/HISAT2 > $logdir/$tool.log
    fi

# SALMON INDEX
elif [ "$1" == "SALMON" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        salmon index -t $ref_cdna -i $outdir --gencode -p 14 > $logdir/$tool.log
    fi

# KALLISTO INDEX
elif [ "$1" == "KALLISTO" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        base_ref_cdna=$(basename $ref_cdna .gz)
        kallisto index -i $outdir/${base_ref_cdna}.idx $ref_cdna > $logdir/$tool.log

    fi
fi