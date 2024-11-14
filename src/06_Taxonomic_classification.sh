### Infer origins of mRNA reads

# Kaiju generation of taxonomic classifications using proGenomes DB
###########################################################################################

# PACKAGE: Kaiju (https://github.com/bioinformatics-centre/kaiju)
# -t: hierarchical representation of taxonomy IDs
# -f: precomputed index for kaiju
# -i: input reads
# -z: number of threads supported on your system
# -o: output file for kaiju taxonomic classifications

# Generation of taxonomic classification

kaiju -t nodes.dmp \
-f kaiju_db.fmi \
-i mouse1_mRNA.fastq \
-z 4 \
-o mouse1_classification.tsv

# Run small python script to restrict specificity of classifications to genus-level taxa
###########################################################################################

# Command input: 3_Constrain_Classification.py <Minimum_Taxonomic_Rank> 
# <kaiju_Classification> <nodes_file> <names_file> <Output_Classifications>

# Activate permission
chmod +x 3_Constrain_Classification_DB.py

# Run script
./3_Constrain_Classification_DB.py genus mouse1_classification.tsv \
kaijudb/nodes.dmp kaijudb/names.dmp mouse1_classification_genus.tsv

# Generate readable summary using Kaiju
###########################################################################################

# PACKAGE: Kaiju 
# -t: hierarchical representation of taxonomy IDs
# -n: taxonomic names corresponding to each taxonomy ID
# -i: Kaiju taxonomic classifications
# -o: summary report output file
# -r: taxonomic rank for which summary report will be produced

kaiju2table \
-t nodes.dmp \
-n names.dmp \
-r genus \
-o mouse1_classification_Summary.tsv \
mouse1_classification_genus.tsv


# Create overview pie chart with Krona
###########################################################################################

# Create output file readable by Krona with Kaiju
kaiju2krona \
-t nodes.dmp \
-n names.dmp \
-i mouse1_classification_genus.tsv \
-o mouse1_classification_Krona.tsv


# PACKAGE: Krona (https://github.com/marbl/Krona/wiki)
# Install Krona
tar -xzf precomputed_files.tar.gz KronaTools
sudo KronaTools/install.pl

# Run Krona
# -o: output file (html)
# input file (output from kaiju2krona)
KronaTools/scripts/ImportText.pl \
-o mouse1_classification.html \
mouse1_classification_Krona.tsv
