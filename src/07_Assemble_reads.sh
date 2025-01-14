### Assemble reads to larger contigs to increase annotation 

# SPAdes genome assembly for generation of putative reads
###########################################################################################

# PACKAGE: SPAdes (https://github.com/ablab/spades)
# --rna: uses mRNA transcript assembly algorithm
# -o: output directory
# -s: single-end input reads

SPAdes-4.0.0-Darwin/bin/spades.py --rna \
-s mouse1_unique_mRNA.fastq \
-o mouse1_spades \
mv mouse1_spades/transcripts.fasta mouse1_contigs.fasta # rename file

# Build an index to allow BWA to search against our set of contigs (named mouse1_contigs.fasta)

bwa index -a bwtsw mouse1_contigs.fasta

# First, make sure the mRNA file has only unique FASTA sequences:

seqkit rmdup -s mouse1_unique_mRNA.fastq -o mouse1_unique_mRNA_deduplicated.fastq

# Map entire set of putative mRNA reads to this newly generated DB

bwa mem -t 4 mouse1_contigs.fasta mouse1_unique_mRNA_deduplicated.fastq > mouse1_contigs.sam

# Run small python script to extract unmapped reads & generate mapping table for each contig
###########################################################################################

# Command input: 5_Contig_Map.py <Reads_Used_In_Alignment> <Output_SAM_From_BWA> \
# <Output_File_For_Unassembed_Reads> <Output_File_For_Contig_Map>

./5_Contig_Map.py mouse1_unique_mRNA_deduplicated.fastq mouse1_contigs.sam \
mouse1_unassembled.fastq mouse1_contigs_map.tsv
