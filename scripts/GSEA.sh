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

gsea-cli.sh GSEAPreranked -gmx data/input-2/DPN_response.gmt -collapse No_Collapse -mode Median_of_probes -norm meandiv -nperm 1000 -rnd_seed 1998 -rnk data/input-2/ranked_DPN.rnk -scoring_scheme weighted -rpt_label dpn_24_analysis_gsea -create_svgs true -include_only_symbols false -make_sets true -plot_top_x 20 -set_max 500 -set_min 15 -zip_report false -out dge_reports/GSEA
