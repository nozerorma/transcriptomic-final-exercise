# INDEX script

tool=$1
#threads=$3

mkdir -p "res/index/$tool"
outdir="res/index/$tool"
ref_gen="data/assembly/reference_grch38/Homo_sapiens.GRCh38.dna*"
ref_cdna="data/assembly/reference_grch38/Homo_sapiens.GRCh38.cdna*"

echo -e "\nBuilding $tool index...\n"

# STAR INDEX

if [ "$1" == "STAR" ]; then
    
    if [ -f "$outdir/*" ]; then
        echo -e "\Index already built, skipping...\n"
    else
        STAR 	--runThreadN 4 --runMode genomeGenerate --genomeDir $outdir \
                --genomeFastaFiles $ref_gen --runRNGseed 1998

        echo -e "$tool index built.\n"
    fi

# HISAT2 INDEX
# From RNAseq hands-on session by Jaime Martínez de Villarreal, Epithelial Carcinogenesis Group, CNIO
elif [ "$1" == "HISAT2" ]; then
    
    if [ -f "$outdir/*" ]; then
        echo -e "\Index already built, skipping...\n"
    else
        hisat2-build -p 6 --seed 1998 $ref_gen $outdir
    fi

# SALMON INDEX
elif [ "$1" == "SALMON" ]; then
    
    if [ -f "$outdir/*" ]; then
        echo -e "\Index already built, skipping...\n"
    else
        echo "echo"
    fi

# KALLISTO INDEX
elif [ "$1" == "KALISTO" ]; then
    
    if [ -f "$outdir/*" ]; then
        echo -e "\Index already built, skipping...\n"
    else
        echo "echo"
    fi
fi