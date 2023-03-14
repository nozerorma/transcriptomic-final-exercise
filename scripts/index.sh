# INDEX script

tool=$1
#threads=$3

mkdir -p "log/index/$tool"
logdir="log/index/$tool"
mkdir -p "res/index/$tool"
outdir="res/index/$tool"
ref_gen="data/assembly/reference_grch38/Homo_sapiens.GRCh38.dna.chromosome.21.fa"
ref_cdna=data/assembly/reference_grch38/Homo_sapiens.GRCh38.cdna.all.fa.gz*

echo -e "\nBuilding "$tool" index...\n"

# STAR INDEX

if [ "$1" == "STAR" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        STAR 	--runThreadN 4 --runMode genomeGenerate --genomeDir $outdir \
                --genomeFastaFiles $ref_gen --runRNGseed 1998 --genomeSAindexNbases 11 > $logdir/log.txt

        echo -e "$tool index built.\n"
    fi

# HISAT2 INDEX
# From RNAseq hands-on session by Jaime Martínez de Villarreal, Epithelial Carcinogenesis Group, CNIO
elif [ "$1" == "HISAT2" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        hisat2-build -p 6 --seed 1998 $ref_gen $outdir/ > $logdir/log.txt
    fi

# SALMON INDEX
elif [ "$1" == "SALMON" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        salmon index -t $ref_cdna -i $outdir --gencode -p 6 > $logdir/log.txt
    fi

# KALLISTO INDEX
elif [ "$1" == "KALISTO" ]; then
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Index already built, skipping...\n"
    else
        base_ref_cdna=$(basename $ref_cdna .fa.gz)
        kallisto --make-unique -i $outdir/${base_ref_cdna}.fa.idx > $logdir/log.txt

    fi
fi