#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the environment
conda activate fibrosis_preprocessing_env

# convert Jupyter notebook(s) to script
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run Python scripts for preprocessing
python nbconverted/0.convert_cytotable.py
python nbconverted/1.sc_quality_control.py
python nbconverted/2.single_cell_processing.py
