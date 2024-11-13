### Remove abundant rRNA sequences

# Covariance model based on Inference of RNA alignment Rfam database
###########################################################################################

# Convert FASTQ to FASTA

vsearch --fastq_filter mouse1_mouse_blat.fastq --fastaout mouse1_mouse_blate.fasta

# PACKAGE: Infernal
# -o: infernal output log file
# --tblout: simple tabular output file
# --noali: omit alignment section --> greatly reduced output volume
# --anytrunc: relaxes tresholds for truncated alignments
# --rfam: uses strict filtering strategy for large DB
# -E: report target sequence with E-value of 0.001
# This can take about 4h for 100'000 reads on single core

cmsearch -o mouse1_rRNA.log --tblout mouse1_rRNA.infernalout \
--anytrunc --rfam -E 0.001 Rfam.cm mouse1_mouse_blat.fasta

# Run small python script to filter out rRNA reads
###########################################################################################

# Command input: 2_Infernal_Filter.py <Input_Reads.fq> <Infernal_Output_File> 
# <mRNA_Reads_Output> <rRNA_Reads_Output>

# Activate permission
chmod +x 2_Infernal_Filter_DB.py

# Run script
./2_Infernal_Filter_DB.py mouse1_mouse_blat.fastq \
mouse1_rRNA.infernalout mouse1_unique_mRNA.fastq mouse1_unique_rRNA.fastq
