#!/bin/bash
#!/bin/bash

if { conda env list | grep 'SRA_pipeline'; } >/dev/null 2>&1; then
	eval "$(conda shell.bash hook)" >/dev/null 2>&1
	conda activate SRA_pipeline >/dev/null 2>&1
	echo "Running conda environment: $CONDA_DEFAULT_ENV"

else
	if [ -d ~/mambaforge ] || [ -d ~/miniforge ]; then
		(echo "Creating conda environment: SRA_pipeline"
		mamba env create -f envs/SRA_pipeline.yaml >/dev/null 2>&1) & spinner $!
		eval "$(conda shell.bash hook)"
		conda activate SRA_pipeline >/dev/null 2>&1
		echo "Running conda environment: $CONDA_DEFAULT_ENV"
	elif [ -d ~/condaforge/ ] || [ -d ~/miniconda ] || [ -d /anaconda3 ]; then
		(echo "Creating conda environment: SRA_pipeline"
		conda env create -f envs/SRA_pipeline.yaml >/dev/null 2>&1) & spinner $!
		eval "$(conda shell.bash hook)"
		conda activate SRA_pipeline >/dev/null 2>&1
		echo "Running conda environment: $CONDA_DEFAULT_ENV"
	fi
fi

# Open-up R env

rstudio DESeq2_differential_analysis.Rproj
