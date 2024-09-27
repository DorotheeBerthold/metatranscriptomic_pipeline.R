#### Metatranscriptomics pipeline - Parkinson lab

# Working directory
cd /Users/dobertho/Desktop/Dorothee_homeoffice/Python_analysis/metatranscriptomics_tutorial

# Get scripts & untar files
wget https://github.com/ParkinsonLab/2017-Microbiome-Workshop/releases/download/Extra/precomputed_files.tar.gz
tar -xvf precomputed_files.tar.gz *.py
tar -xvf precomputed_files.tar.gz mouse1.fastq

# Install fastqc
conda install fastqc


# create fastqc file from raw data
fastqc mouse1.fastq
