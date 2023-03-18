## Ejercicio final de transcriptómica (curso 2022-2023)

## RNA-seq analysis using an interactive pipeline

### RESOLUCIÓN APARTADO 1

Para responder a las preguntas del apartado 1, se decidió realizar un pipeline interactivo en el cual el usuario, tras introducir una o varias entradas de cualquier muestra en formato SRA, puede ir desarrollando diferentes tareas en función del workflow deseado.

Específicamente, el usuario tiene la posibilidad de:

0. Menú interactivo para navegar a través de los diferentes workflows.
```bash
echo -e "${GREEN}\nWORKFLOW PIPELINE MENU\n${NC}
${RED}Default cores set to 14, change accordingly\n${NC}"

echo -e "
${YELLOW}BASIC WORKFLOWS FOR HIGH QUALITY READS, NO PRE-PROCESSING\n${NC}
\tWORKFLOW 1 (Using STAR aligner)\n
\tWORKFLOW 2 (Using HISAT2 aligner)\n
\tWORKFLOW 3 (Using SALMON pseudo-aligner)\n
\tWORKFLOW 4 (Using KALLISTO pseudo-aligner)\n
${YELLOW}ADVANCED WORKFLOWS INCLUDING PREPROCESSING\n${NC}
\tWORKFLOW 5 (Toolset: FastQScreen, Cutadapt, STAR)\n
\tWORKFLOW 6 (Toolset: FastQScreen, Cutadapt, HISAT2)\n
\tWORKFLOW 7 (Toolset: FastQScreen, Cutadapt, SALMON)\n
\tWORKFLOW 8 (Toolset: FastQScreen, Cutadapt, KALLISTO)\n
\t${RED}(Trimmomatic workflows not working as of now)${NC}\n
\tWORKFLOW 9 (Toolset: FastQScreen, Trimmomatic, STAR)\n
\tWORKFLOW 10 (Toolset: FastQScreen, Trimmomatic, HISAT2)\n
\tWORKFLOW 11 (Toolset: FastQScreen, Trimmomatic, SALMON)\n
\tWORKFLOW 12 (Toolset: FastQScreen, Trimmomatic, KALLISTO)\n"

read -rp "Option: " menuOp
```

1. Input interactivo de entrada SRA, descarga de fichero, y dumping de ficheros forward y reverse utilizando la herramienta 'Fasterq-dump' `download.sh`.
```bash
while true; do
  read -a SRAentries  # Read a line of input and split it into an array of accessions
  if [[ "${SRAentries[0]}" =~ ^[SsNn0]$ ]]; then
    # If the first element is 'S', 's', 'N', 'n' or '0', exit the loop
      exit 1
  
  else
    for SRAentry in "${SRAentries[@]}"; do
      if [[ "$SRAentry" =~ ^[Ss]RR[0-9]{6}$ ]]; then
        echo "Matched SRA accession $SRAentry"
        (wget -nc -P $sample_dir/SRA --content-disposition https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRAentry/$SRAentry) & spinner
        echo "Finished downloading $SRAentry"
				
      # Dump fastq from SRA entry
      mkdir -p $sample_dir/dumped_fastq
      echo -e "\nDumping fastq files from $SRAentry...\n"
        if [ "$(ls -A $sample_dir/dumped_fastq/$SRAentry)" ]; then
          echo -e "Fastq already dumped for $SRAentry, skipping.\n"

        else
          (fasterq-dump -fp -O $sample_dir/dumped_fastq/$SRAentry $SRAentry) & spinner $!
        fi				
            
      else
        echo "Did not match SRA accession $SRAentry. Skipping download."
      fi
    done
    break
  fi
done
```

2. Análisis de calidad de muestras descargadas, utilizando las herramientas 'FastQC' y 'Fastqscreen'. Incluye la posibilidad de descargar genomas de referencia de FastQScreen `qc.sh`.
```bash
# FASTQC
	
echo -e "${YELLOW}\nPerforming QC analysis for $nofastq_sid...
___________________________________________________________ ${NC}\n"

if [ -f $fastqc_dir/$cut_sid/${nofastq_sid}_fastqc.html ]; then
	echo -e "FastQC analysis already performed for $nofastq_sid, skipping analysis.\n" 

else
	echo -e "Running FastQC analysis..."
	mkdir -p $fastqc_dir/$cut_sid
	mkdir -p log/qc
	(fastqc -o $fastqc_dir/$cut_sid $sid 2>&1 >/dev/null | tail -n +2 2>&1 log/qc/fastqc.log) & spinner $!
	echo
fi

# FASTQSCREEN
if [ -f $fastqscreen_dir/$cut_sid/${nofastq_sid}_screen.html ]; then
	echo -e "FastQScreen analysis already performed for $nofastq_sid, skipping analysis.\n" 
	
else
	screen_gen="res/fastq_screen_samples/FastQ_Screen_Genomes"
	if [ ! "$(ls -A $screen_gen)" ]; then
		echo "No reference genomes could be found."
		echo -e "\nWould you like to download genome indexes from database?"
		echo -e "It takes long AND it takes space. Note that if you don't download the references,
you have to manually curate fastqscreen.conf file and point the program to it" 
		read -rp "Continue? (Y/n): " genDownload
		
		case $genDownload in
			[Yy]* )
				echo -e "\nDownloading genomes...\n" 
				fastq_screen --get_genomes --outdir "res/fastq_screen_samples"
			;;
			[Nn]* )
				echo -e "\nSkipping genome download...\n"
			;;
		esac
	
	else
		echo -e "\nRunning FastQScreen analysis..."
		mkdir -p $fastqscreen_dir/$cut_sid
		mkdir -p log/qc
		(fastq_screen --conf "$screen_gen/fastq_screen.conf" \
			--tag --aligner bowtie2 --subset 100000 --threads 14 \
			--outdir "$fastqscreen_dir/$cut_sid" $sid 2>&1 log/qc/fastqscreen.log) & spinner $!
		echo
	fi
fi
```

3. Pre-procesado de las muestras tras observación de análisis de calidad, utilizando las herramientas 'Cutadapt' (incluye adaptadores más comunes en secuenciación Illumina RNA) y 'Trimmomatic' (**TO-DO**) `pre_proc.sh`. Posibilidad de realizar un nuevo QC de las muestras curadas.
```bash
if [ "$1" == "cutadapt" ]; then
    
    echo -e "\nRunning cutadapt...\n" 
    
    if [ "$(ls -A $outdir)" ]; then
        echo -e "Trim already performed for $f_sid, skipping...\n"
    
    else
        # adapters infered from most common Illumina, listed in here
        # https://tinyurl.com/illumina-adapters
        # https://tinyurl.com/illumina-adapters-2
        (cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
            -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
            -a CTGTCTCTTATACACATCT...AGATGTGTATAAGAGACAG -a TGGAATTCTCGGGTGCCAAGG \
            -a "G{100}" -g "G{100}" -a "A{100}" -g "A{100}" -q 20 \
            -o "$outdir/${f_sid}_trimmed.fastq" -p "$outdir/${r_sid}_trimmed.fastq" \
            "$f_path" "$r_path" --cores 14 -m 1 --discard-untrimmed > "$logdir"/log.txt)
    fi
fi

# Re-run quality control
read -rp "Would you like to re-run QC report for your trimmed samples? (Y/n) " runqc

for trimmed_sid in $(find $outdir -type f -name "*_trimmed.fastq");do
    case $runqc in
        [Yy]* )
            bash scripts/qc.sh $trimmed_sid "out/qc/fastqc" "out/qc/fastq_screen"
        ;;
	esac
done

case $runqc in
[Yy]* )
	printf "${GREEN}\nNow, take your time to give a look to the QC analysis.\n${NC}"
	echo "When you are ready, press any key to continue..."
	read -n 1 -s -r -p ""
;;
esac
```

4. Construcción de Index para las herramientas de alineamiento y pseudoalineamiento utilizadas (STAR, HISAT2, SALMON y KALLISTO `index.sh`.
```bash
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
```

5. Alineamiento (o pseudoalineamiento de los reads), con la posibilidad de seguir con workflows anteriores en los que se ha realizado trimming de las mismas. Mismos alineadores anteriormente referidos `align.sh`. Incluye también conversiónd de ficheros SAM a sorted.BAM, así como estadísitcas de los diferentes estadíos de alineamiento (SAMTOOLS), y extracción de coverage en formato .bw (bamCoverage). 
```bash
# First lets see what workflow we are working with

sample_dir=0

if [ "$workflow" == "untrimmed" ]; then
	sample_dir="$untrimmed_dir"

elif [ "$workflow" == "cutadapt" ] || [ "$workflow" == "trimmomatic" ] ; then
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
```

6. Conteo de los reads y elaboración  de matriz de cuentas, utilizando HTSeq o countFeatures. Se puede realizar sobre los ficheros bam ordenados de STAR, HISAT2 y KALLISTO (*post_proc.sh*).
```bash
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
```

**NÓTESE QUE, SI BIEN TODOS LOS DIFERENTES SCRIPTS ESTÁN PENSADOS PARA SER CORRIDOS EN CONJUNTO A TRAVÉS DEL MENÚ INTERACTIVO PRESENTE EN  _pipeline.sh_, PUEDEN SER CORRIDOS INDEPENDIENTEMENTE. PARA VER LOS ARGUMENTOS NECESARIOS BASTA CON TRATAR DE CORRER CUALQUIERA DE LOS SCRIPTS.**


### DESCRIPCIÓN DEL EJERCICIO

El objetivo de este trabajo es que el alumno demuestre los conocimientos obtenidos acerca del análisis de RNA-seq. Para ello deberá **redactar un informe** en el que se expliquen los datos de partida y se extraiga una **conclusión de los resultados**. Así mismo, el alumno debe **detallar el proceso de análisis** indicando el *software* (incluída la versión) empleado, así como los parámetros utilizados en cada uno de los pasos. En caso de que hubiera que eliminar muestras por motivos técnicos o biológicos, el alumno debe indicar y justificar el por qué en cada caso. Se pueden introducir *code chunks* e imágenes para apoyar el informe. De manera alternativa, puede aportarse todo el código en forma de repositorio público. Se han planteado **5 preguntas (10 puntos en total)** para guiar la redacción del informe.

El trabajo consta de dos apartados en los que se utilizarán datos de un experimento en el que disponemos de 24 cultivos primarios de tumores paratiroideos negativos para receptores de estrógenos alfa (ERα). Las muestras, procedentes de 4 pacientes diferentes, se han tratado con dos fármacos diferentes: diarilpropionitrilo (DPN) o 4-hidroxitamoxifeno (OHT) a 24h o 48h. El DPN es un agonista del ERα mientras que el OHT es un inhibidor competitivo de los receptores de estrógenos.

El **primer apartado (3 preguntas)** abarca los pasos de control de calidad y de fuentes de contaminación, *trimming*, alineamiento y cuantificación para obtener cuentas crudas y normalizadas a partir de un **subset de ficheros fastq**. El **segundo apartado (2 preguntas)** parte de la **matriz completa de cuentas crudas** y está enfocado a realizar un control de calidad biológico, detectar los genes diferencialmente expresados entre condiciones y los pathways enriquecidos en cada una de ellas.

#### Apartado 1

El dataset original consta de 27 muestras paired-end depositadas en SRA. Con el fin de poder abordar las cuestiones planteadas a en el primer apartado sólo es necesario descargar 2 muestras (SRR479052 y SRR479054), es decir cuatro ficheros fastq. Además, en este repositorio se proporciona:

- Un fichero fasta con la secuencia de la referencia genómica, en este caso correspondiente al cromosoma 21 humano (ensamblaje GRCh38).
- Un fichero GTF con la anotación génica para los genes del cromosoma 21 (GRCh38.ensembl.109).

Se pide realizar un análisis de dichas muestras similar al realizado en clase, considerando los siguientes puntos:

**Pregunta 1 (1.5 puntos):** Realizar control de calidad de dichas muestras con el programa FastQC, incluir plots más reseñables y comentar cada uno de los apartados. De manera complementaria se podrá realizar un analisis de contaminación con el programa FastQScreen.

**Pregunta 2 (1.5 puntos):** Para poder llevar a cabo el alineamiento de las muestras en vuestros ordenadores será necesario trabajar con archivos reducidos correspondientes al cromosoma 21. Se requiere el indexado de la secuencia de este cromosoma, así como el alineamiento de las muestras a dicha referencia. Para ello se podrá utilizar el alineador HISAT2 utilizado en clase u otros (alineadores o pseudoalineadores). Comentar cada uno de los comandos y parámetros empleados, justificando su uso.

**Pregunta 3 (1.5 puntos):** Una vez generados los archivos alineados se reportarán las estadísticas de alineamiento y se procederá a la cuantificación de la expresión utilizando el archivo GTF correspondiente. Para ello se podrá utilizar HTSeq u otras herramientas de cuantificación. En cualquier caso, detallar y justificar los comandos y parámetros empleados para ello.

#### Apartado 2

En este repositorio se proporcionan todos los inputs del Apartado 2:

- La matriz de cuentas crudas para los 24 cultivos analizados. 
- Data frame con los metadatos asociados al experimento.
- GMT para realizar un GSEA.

**Pregunta 4 (3 puntos):** ¿Qué genes se encuentran diferencialmente expresados entre las muestras pertenecientes al grupo tratado con OHT con respecto al control tras 24h? ¿Y en el caso de las muestras tratadas con DPN, también tras 24h? Como parte de la respuesta a esta pregunta, podéis entregar una o más tablas adjuntas donde se incluyan los genes diferencialmente expresados, indicando el criterio o los criterios que habéis seguido para filtrar los resultados, así como cualquier otro gráfico o gráficos que hayáis elaborado durante el análisis.

**Pregunta 5 (2.5 puntos):** Nuestro colaborador ha comparado las muestras tratadas con DPN tras 48h con las muestras control. Nos ha llamado para contarnos que los cambios de expresión tras 48h son mucho más evidentes y nos preguntamos si el DPN produce algún efecto en las primeras 24h. Para contestar esta pregunta, le pedimos que genere un GMT (input) con los genes más expresados en las muestras tratadas tras 48h (DPN_perturbed) y los genes más expresados en la muestra control (DPN_unperturbed). Realiza un análisis con GSEA para determinar el efecto del tratamiento a las 24h. ¿A qué conclusión llegas? Incluid una tabla con los resultados del análisis, destacando las columnas con los valores utilizados para extraer vuestras conclusiones. También incluid los gráficos característicos de este tipo de análisis.
