## Ejercicio final de transcriptómica (curso 2022-2023)

## RNA-seq analysis using an interactive pipeline

### MÉTODOS

Para dar respuesta a las diferentes preguntas que propone el ejercicio, se decidió realizar un pipeline interactivo, que permitiese al usuario introducir cualquier muestra SRA, la cual sería descargada, desempaquetada en los correspondientes archivos FASTQ, analizada, preprocesada, alineada y post-procesada utilizando diferentes metodologías.

Esto mismo puede encontrarse en el fichero `INFORME.md`presente en este mismo repositorio, donde también se encuentran las discusiones sobre los resultados observados.

#### Descarga del genoma `download.sh`

Se utilizó la herramienta de sistema *wget*, integrando un *input control* que confirme la correcta introducción del formato de SRA por parte del usuario (**if [[ "SRAentry" =~ ^[Ss]if [[ "SRAentry" =~ ^[Ss]RR[0-9]{6}RR[0-9]{6} ]]**), y los parámetros específicos de *wget* **-nc (no clobber, evita descargas en duplicado)**, **-P (designa el directorio de descarga)**, y **--content-disposition (mantiene la integridad del formato de descarga, evitando el que por error descargue la el enlace como HTML)**.

```shell
 if [[ "$SRAentry" =~ ^[Ss]RR[0-9]{6}$ ]]; then
                echo "Matched SRA accession $SRAentry"
                (wget -nc -P $sample_dir/SRA --content-disposition https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRAentry/$SRAentry) & spinner
                echo "Finished downloading $SRAentry"
```

Al igual que para el resto de procesos que se realizan a lo largo del pipeline, se añadió un pequeño spinner o reloj para embellecer la espera y para dar al usuario seguridad de que el programa no se ha congelado (ver `& spinner` y `spinner.sh`. )

#### Desempaquetado de ficheros FASTQ `download.sh`

Se utilizó la herramienta fasterq-dump para desempaquetar el fichero SRA en los consecuentes fFASTQ y rFASTQ, introduciendo una comprobación de existencia del fichero. Los parámetros utilizados fueron **-fp** (forzar extracción y mostrar barra de procesamiento).

```shell
if [ "$(ls -A $sample_dir/dumped_fastq/$SRAentry)" ]; then
                        echo -e "Fastq already dumped for $SRAentry, skipping.\n"

                    else
                        (fasterq-dump -fp -O $sample_dir/dumped_fastq/$SRAentry $SRAentry) & spinner $!
```

Una vez descargados, se le pregunta al usuario si desea realizar un control de calidad de los FASTQ brutos. En caso afirmativo se procede al paso de QC; en caso contrario, se procede directamente al menú general.

```shell
read -rp "Would you like to perform a QC analysis (FastQC and FastQScreen)? (Y/n): " performQC

for SRAentry in "${SRAentries[@]}"; do
    for sid in $(find $sample_dir/dumped_fastq/$SRAentry -type f -name '*.fastq'); do
        base_sid=$(basename $sid | cut -d"_" -f1)

        if [ "$(ls -A $sample_dir/dumped_fastq/$base_sid)" ]; then
            case $performQC in
                [Yy]* )
                    bash scripts/qc.sh $sid "out/qc/fastqc" "out/qc/fastq_screen"
                ;;
                [Nn]* )
                    echo -e "\nSkipping QC analysis...\n"
                    break 
                ;;
            esac
        fi
    done
done


# Results visualization
case $performQC in
[Yy]* )
    printf "${GREEN}\nNow, take your time to give a look to the QC analysis.\n${NC}"
    echo "When you are ready, press any key to continue..."
    read -n 1 -s -r -p ""
;;
esac
```

#### OPT 1: Realizar control de calidad de FASTQ brutos `qc.sh`

Para el control de calidad de las muestras, se utilizó la combinación de programas *fastQC* y *fastq-screen*. El script tiene dos posibilidades de uso; por una parte, está pensado para ser utilizado dentro de la ejecución del pipeline `pipeline.sh`, como se ha descrito hasta ahora en la metodología; por otra, y al igual que el resto de herramientas que se describan en adelante, permite un uso independiente con cualquier fichero correspondiente que disponga el usuario, mediante la siguiente porción de código:

```shell
# Comprobation for independent use of script
if [ "$#" -ne 3 ]
then
    printf "${RED}Usage: $1 <fastq_file> $2 <fastqc_out_dir> $3 <fastqscreen_out_dir> ${NC}\n"
    exit 1
fi 
```

En este caso, el programa comprueba la inserción de 3 parámetros junto con la ejecución del código (sucede por defecto con el pipeline general), y en caso contrario, avisa al usuario del uso correcto de la herramienta. El ${RED} hace referencia a una serie de variables declaradas en la cabecera de cada script, indicando diferentes colores que van a ser usados para dar mayor comprehensibilidad al código.

##### FASTQC

Realiza las comprobaciones pertinentes y después ejecuta el programa *fastQC* con los parámetros **-o (directorio de salida del informe generado)** y **2>&1 >/dev/null | tail -n +2 2>&1 log/qc/fastqc.log (redirige la salida de consola a un fichero externo y elimina una línea innecesaria que puede dar lugar a confusión).**

```shell
fastqc -o $fastqc_dir/$cut_sid $sid 2>&1 >/dev/null | tail -n +2 2>&1 log/qc/fastqc.log
```

##### FASTQ-SCREEN

Realiza las comprobaciones pertinentes

```shell
if [ -f $fastqscreen_dir/$cut_sid/${nofastq_sid}_screen.html ]; then
    echo -e "FastQScreen analysis already performed for $nofastq_sid, skipping analysis.\n" 
```

y otorga al usuario la posibilidad de descargar los índices de referencia utilizados por la herramienta para el uso de la misma, añadiendo un aviso para el usuario sobre el peso de los mismos y el tiempo de descarga que requieren

```shell
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
```

para después correr la herramienta utilizando los parámetros **--conf (apunta al fichero de configuración)**, **--tag (etiqueta las lecturas en función del indice con el que casen)**, -**-aligner (especifica el alineador a utilizar, el cual debe corresponderse con el usado para la generación de los índices; en el caso de los usados de serie por *fastqscreen*, bowtie2 )**, **--subset (tamaño de sampleo, para no usar todo el fastq)**,**--threads (numero de cores a utilizar)**, y **--outdir (indica el directorio de salida).**

```shell
fastq_screen --conf "$screen_gen/fastq_screen.conf" \
    --tag --aligner bowtie2 --subset 100000 --threads 14 \
    --outdir "$fastqscreen_dir/$cut_sid" $sid 2>&1 log/qc/fastqscreen.log) & spinner $!
```

#### OPT 2: Menú de Opciones y elección de workflow `pipeline.sh`

Tanto como si se ha decidido no realizar un control de calidad de las muestras como si se ha hecho, el siguiente paso consiste en el procesamiento de la muestra a través de numeras opciones de workflow en función de los intereses del investigador, a través de un menú interactivo.

Éste se encuentra dividido en dos categorías; por una parte, workflows sin preprocesamiento para alineamiento rápido en muestras de muy alta calidad; y por otra, diferentes workflows que incluyen pre-procesado previo de las mismas.

```shell
# Different pipelines for different workflows
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
\tWORKFLOW 9 (Toolset: FastQScreen, Fastp, STAR)\n
\tWORKFLOW 10 (Toolset: FastQScreen, Fastp, HISAT2)\n
\tWORKFLOW 11 (Toolset: FastQScreen, Fastp, SALMON)\n
\tWORKFLOW 12 (Toolset: FastQScreen, Fastp, KALLISTO)\n"

read -rp "Option: " menuOp
```

##### Herramientas de pre-procesamiento

###### Cutadapt

Cutadapt es una herramienta que permite el cribado de las lecturas por calidad, así como la remoción de adaptadores usados en la secuenciación. Debido al tipo de experimento realizado, se utilizó con los siguientes parámetros (los cuales deben ser modificados acordes al experimento a analizar):

- **-a, -A, -g y -G** representan respectivamente adaptadores hallados en los extremos 5' y 3' de las muestras forward (minúscula) y reverse (mayúscula). Se eliminaron los diferentes adaptadores encontrados en el análisis FasQC + típicos en secuenciación Illumina PE (enlaces), así como residuos PolyA y PolyG.

- **-q** es el umbral de calidad que debe superar una lectura para ser considerada válida. Se fijó en 15 por concordancia con el umbral por defecto de FastP.

- **-o** y **-p** son las rutas de salida del fichero procesado, seguidas en este caso de los ficheros a procesar (f_path y r_path).

- **--cores** fija el número de nucleos

- **-m** fija la longitud mínima que debe tener una lectura para ser considerada válida. Se fijó en 15 por concordancia con el umbral por defecto de FastP.

- **--discard-untrimmed** toma como válidas solo aquellas lecturas que hayan sido preprocesadas.

```shell
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
            -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG \
            -a CTGTCTCTTATACACATCT...AGATGTGTATAAGAGACAG -a TGGAATTCTCGGGTGCCAAGG \
            -a "G{10}" -A "G{10}" -g "G{10}" -G "G{10}" -g "A{10}" -G "A{10}" -q 15 \
            -o "$outdir/${f_sid}_trimmed.fastq" -p "$outdir/${r_sid}_trimmed.fastq" \
            "$f_path" "$r_path" --cores 14 -m 15 --discard-untrimmed > "$logdir"/log.txt
```

Tras esto, se le da la posibilidad al usuario de realizar un nuevo control de calidad de las muestras procesadas.

```shell
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

###### FastP

Herramienta de control de calidad y preprocesado muy completa, que permite automatizar la eliminación de adaptadores y cribado de las lecturas sin necesidad de hardcodear las secuencias. Se usó con los siguientes parámetros:

- **-i, -o, -I** y **-O** representan los ficheros de entrada (i) y salida (o), para las lecturas forward (minúscula) y reverse (mayúscula).

- **--detect_adapter_for_pe** detecta automáticamente adaptadores usados en PE.

- **5** y **3** realizan los cortes en ambos extremos de la lectura, en función de si las bases en estas superan un umbral de calidad determinado.

- **-D** realiza deduplicación de las muestras.

- **--dup_calc_accuracy** fijado en 3 como calidad por defecto en deduplicación.

- **-c** corrige errores de solapamiento en muestras PE.

- **-p** realiza un análisis de secuencias sobrerrepresentadas.

- **-h** exporta los resultados en formato HTML para mejor visualización que en JSON por defecto.

- **-w** fija el número de cores a utilizar.

```shell
fastp -i $f_path -o "$outdir/${f_sid}_trimmed.fastq" -I $r_path -O "$outdir/${r_sid}_trimmed.fastq" \
 --detect_adapter_for_pe -5 -3 -D --dup_calc_accuracy 3 -c -p -P 20 \
 -h "$outdir/${base_sid}.html" -w 14 > "$logdir"/log.txt
```

##### Herramientas de indexado y alineamiento

Para ver una descripción detallada de las diferencias entre las distintas herramientas de alineamiento, echar un vistazo al manual incluído con las mismas.

Uso independeinte del script de indexado

```shell
# Comprobation for independent use of script
if [ "$#" -ne 1 ]
then
    printf "${RED}Usage: $1 <tool> ${NC}\n"
    echo -e 'tool: "STAR", "HISAT2", "SALMON", "KALLISTO"\n'
    exit 1
fi
```

y de alineamiento:

```shell
# Comprobation for independent use of script
if [ "$#" -ne 6 ]
then
    printf "${RED}Usage: $1 <tool> $2 <f_path> $3 <r_path> $4 <outdir> $5 <logdir> $6 <workflow>${NC}\n"
    echo -e 'tool: "STAR", "HISAT2", "SALMON", "KALLISTO"
workflow: "untrimmed", "cutadapt", "fastp"\n'
    exit 1
fi
```

###### STAR

Parámetros utilizados para el indexado de las muestras:

- El único parámetro a remarcar es **-genomeSAindex**, el cual fue fijado a partir de las recomendaciones dadas por la propia herramienta en su primera ejecución. Modificar acorde a las condiciones del experimento.

```shell
STAR     --runThreadN 14 --runMode genomeGenerate --genomeDir $outdir \
                --genomeFastaFiles $ref_gen --runRNGseed 1998 --genomeSAindexNbases 11 > $logdir/$tool.log
```

Parámetros utilizados para el alineamiento:

A remarcar el parámetro **--outSAMtype** fijado en sorted BAM que nos permite ahorrarnos pasos en la conversión del SAM al BAM. Ejecución posterior de *samtools* para generar index del fichero de alineamiento y diferentes estadísticas del mismo (generadas con *samtools* y *bamCoverage*).

```shell
STAR --runThreadN 14 --readFilesIn "$f_path" $r_path  \
            --genomeDir "$index_dir" --outReadsUnmapped Fastx  \
            --outFileNamePrefix "$outdir/" \
            --outSAMtype BAM SortedByCoordinate

# add more params for statistics (there's a few problems here to solve)
bam_file=$(find "$outdir" -type f -name "*.bam")
(samtools stats "$bam_file" > "$bam_file".txt) & spinner $!
(samtools index "$bam_file") & spinner $!
(bamCoverage -b "$bam_file" -o "$outdir/$base_sid.bw" --normalizeUsing BPM) & spinner $!
```

###### HISAT2

Ningún parámetro a recalcar en la generación del índice:

```shell
hisat2-build -p 14 --seed 1998 $ref_gen $outdir/HISAT2
```

El alineamiento de la muestra requiere de un mayor número de pasos al utilizar *hisat2* en su conversión a BAM.

```shell
(hisat2 --new-summary --summary-file "$outdir/$base_sid.hisat2.summary" \
    -p 14 -x "$index_dir/HISAT2" -1 "$f_path" -2 "$r_path" -k 1 -S "$outdir/$base_sid.sam") & spinner $!

# problema con el pipe
# add more params for statistics
(samtools view -bS "$outdir/$base_sid.sam" > "$outdir/$base_sid.bam") & spinner $!
(samtools stats "$outdir/$base_sid.bam" > "$outdir/$base_sid.txt") & spinner $!
(samtools sort --write-index "$outdir/$base_sid.bam" -o "$outdir/$base_sid.sorted.bam") & spinner $!
(samtools stats "$outdir/$base_sid.sorted.bam" > "$outdir/$base_sid.sorted.txt") & spinner $!
(bamCoverage -b "$outdir/$base_sid.sorted.bam" -o "$outdir/$base_sid.bw" --normalizeUsing BPM) & spinner $!
```

###### SALMON

Tanto *SALMON* como *KALLISTO* son pseudoalineadores que trabajan de una manera diferente a los dos previamente descritos.

Generación del índice (utilizando exoma de referencia en vez de genoma de referencia):

```shell
salmon index -t $ref_cdna -i $outdir --gencode -p 14 > $logdir/$tool.log
```

Alineamiento de las muestras:

- **-l** es el tipo de stranding de la muestra. Fijado en A para determinación automática.
- **-1 y -2** son los ficheros de entrada.
- **-g** da la opción de añadir fichero de anotación.
- **--writeMappings** escribe los alineamientos en formato SAM.

```shell
salmon quant -i "$index_dir" -l A -1 "$f_path" -2 "$r_path" --validateMappings \
    -o "$outdir" -g "$ref_gtf" -p 14 --writeMappings 
```

###### KALLISTO

Generación de índice:

```shell
kallisto index -i $outdir/${base_ref_cdna}.idx $ref_cdna > $logdir/$tool.log
```

Pseudoalineamiento:

- **--bias** trata de corrigir errores de sesgo.

- **--fusion** busca fusiones en Pizzly.

- **--pseudobam** escribe los pseudoalinamientos con el transcriptoma a un fichero con formato BAM.

- **--genomebam** projecta los pseudoalineamientos frente al exoma de referencia.

- **--gtf** permite la inserción de anotaciones en formato GTF.

```shell
kallisto quant -i "$index_dir/$base_ref_cdna.idx" --bias --fusion \
    "$f_path" "$r_path" -o "$outdir" --pseudobam --genomebam \
    --gtf "$ref_gtf" -t 14

    pseudobam="$outdir/pseudoalignments.bam"
    bamCoverage -b "$pseudobam" -o "$outdir/$base_sid.pseudoalignments.bw" \
    --normalizeUsing BPM 2>&1
```

Post-procesado y conteo de muestras (read summarization programs)

Uso independiente del script:

```shell
if [ "$#" -ne 3 ]
then
    printf "${RED}Usage: $1 <tool> $2 <bam_file> $3 <counts_out_dir>${NC}\n"
    echo -e 'tool: "featurecounts", "htseq"\n'
    exit 1
fi
```

###### FEATURECOUNTS

```shell
featureCounts -T 14 -s 2 -p --countReadPairs -a "$ref_gtf" -t exon -g gene_id \
    -o "$countdir/counts.txt" "$bam"
```

###### HTSEQ

```shell
htseq-count --format=bam --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id    \
    --additional-attr=gene_name "$bam" "$ref_gtf" > "$countdir/$base_bam.htseq"
```



## APARTADO 2

Para la realización del apartado 2, se recomienda encarecidamente utilizar el script `deseq2.sh`, el cual realiza una comprobación automática de los paquetes necesarios y abre una nueva sesión de RStudio en el entorno de R adecuado.

Los scripts usados para la generación de los ficheros necesarios para el GSEA se encuentran también presentes en este entorno, mientras que la herramienta como tal debe ser corrida utilizando el script `gsea.sh`, ya que también realiza una comprobación de los paquetes necesarios para su uso.
