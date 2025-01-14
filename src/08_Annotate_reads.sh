### Annotate reads to known genes/proteins

# BWA annotation of putative reads to known genes
# SUBMIT TO EULER
###########################################################################################

# BWA has high stringency, but sequence diversity at nucleotide level results in fewer matches
# microbial genome database curated from NCBI: microbial_all_cds.fasta (9GB)
# followed by BLAT (as in 04 & 05) to reduce using DIAMOND which uses much larger DB (>60GB)

# Get all bacterial sequences from NCBI & concatenate them into a single file
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.ffn.tar.gz
tar -xzvf all.ffn.tar.gz
cat *.ffn > microbial_all_cds.fasta

# Index the microbial genome database
bwa index -a bwtsw microbial_all_cds.fasta
samtools faidx microbial_all_cds.fasta

# Search all assembled contigs against microbial genome database & generate annotation.sam
# mouse1_contigs.fasta is the index built with BWA from 08
bwa mem -t 4 microbial_all_cds.fasta mouse1_contigs.fasta > mouse1_contigs_annotation_bwa.sam

# Convert FASTQ to FASTA
vsearch --fastq_filter mouse1_unassembled.fastq --fastaout mouse1_unassembled.fasta

# Search all unassembled contigs against microbial genome database & generate unassembled_annotation.sam
bwa mem -t 4 microbial_all_cds.fasta mouse1_unassembled.fasta > mouse1_unassembled_annotation_bwa.sam

# Perform additional annotation with BLAT
# SUBMIT TO EULER
###########################################################################################

# Convert to FASTA
samtools fasta mouse1_contigs_annotation_bwa.sam > mouse1_contigs_bwa.fasta
samtools fasta mouse1_unassembled_annotation_bwa.sam > mouse1_unassembled_bwa.fasta


# Search with BLAT for assembled contigs
blat -noHead -minIdentity=90 -minScore=65 microbial_all_cds.fasta \
mouse1_contigs_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_contigs_bwa_blat.blatout

# Search with BLAT for unassembled contigs
blat -noHead -minIdentity=90 -minScore=65 microbial_all_cds.fasta \
mouse1_unassembled_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_unassembled_bwa_blat.blatout

# Convert to SAM format using blast2sam.pl
#CHECK AGAIN - could not test because I did not have the files from the annotations

python blast_to_SAM.py mouse1_contigs_bwa_blat.blatout mouse1_contigs_annotation_bwa_blat.sam mouse1_contigs_bwa.fasta

python blast_to_SAM.py mouse1_unassembled_bwa_blat.blatout mouse1_unassembled_annotation_bwa_blat.sam mouse1_unassembled_bwa.fasta

# Run small python script to extract high confidence alignments
###########################################################################################

# Activate permission
chmod +x 5_BWA_Gene_Map_DB.py

# Command input: 6_BWA_Gene_Map.py <Gene_database> <Contig_Map> <Output_File_For_Gene_Map> 
# <Output_File_For_Gene_sequences> <Contigs_File> <Contig_BWA_SAM> <Unmapped_Contigs> 
# <Unassembled_Reads_File> <Unassembled_Reads_BWA_SAM> <Unmapped_Unassembled_Reads>

# mouse1_unassembled.fastq is from 08_assembling_reads

./5_BWA_Gene_Map_DB.py microbial_all_cds.fasta mouse1_contigs_map.tsv mouse1_genes_map.tsv \
mouse1_genes.fasta mouse1_contigs.fasta mouse1_contigs_annotation_bwa_blat.sam mouse1_contigs_unmapped.fasta \
mouse1_unassembled.fastq mouse1_unassembled_annotation_bwa_blat.sam mouse1_unassembled_unmapped.fasta


# DIAMOND for annotation based on peptides
# SUBMIT TO EULER
###########################################################################################

# DIAMOND is less prone to sequence-changes between strains for annotation
# Uses much larger DB based on non-redundant (NR) proteins

# Get the database
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# uncurl
curl -O ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

# Build NR database
# -p 8: This specifies the number of threads to use (in this case, 8 threads).
# --in nr: This specifies the input FASTA file (nr).
# -d nr: This specifies the name of the DIAMOND database output (nr).
diamond makedb -p 8 --in nr -d nr

# PACKAGE: DIAMOND (https://github.com/bbuchfink/diamond)
# -p: Number of threads to use in the search is 4.
# -q: Input file name.
# -d: Database name.
# -e: Expectation value (E) threshold for saving hits.
# -k: Maximum number of aligned sequences to keep is 10.
# -t: Temporary folder.
# -o: Output file name.
# -f: Output file is in a tabular format.

# Make new temporary directory
mkdir -p dmnd_tmp

# Annotate all assembled & unmapped fasta from BWA & BLAT
diamond blastx -p 4 -d nr -q mouse1_contigs_unmapped.fasta -o mouse1_contigs.dmdout \
-f 6 -t dmnd_tmp -k 10 --id 85 --query-cover 65 --min-score 60

# Annotate all unassembled & unmapped fasta from BWA & BLAT
diamond blastx -p 4 -d nr -q mouse1_unassembled_unmapped.fasta -o mouse1_unassembled.diamondout \
-f 6 -t dmnd_tmp -k 10 --id 85 --query-cover 65 --min-score 60

# Run small python script to extract top matched proteins
###########################################################################################

# Activate permission
chmod +x 6_Diamond_Protein_Map_DB.py

# Command input: 6_Diamond_Protein_Map_DB.py <Protein_database> <Contig_Map> <Gene_Map> \
# <Output_File_For_Protein_sequences> <Unmapped_Contigs_File> <Contig_Diamond_Output> \
# <Output_For_Unannotated_Contigs> <Unmapped_ Unassembled_Reads_File> <Unassembled_Reads_Diamond_Output> \
# <Output_For_Unannotated_Unassembled_Reads>

./6_Diamond_Protein_Map_DB.py nr mouse1_contigs_map.tsv mouse1_genes_map.tsv \
mouse1_proteins.fasta mouse1_contigs_unmapped.fasta mouse1_contigs.dmdout \
mouse1_contigs_unannotated.fasta mouse1_unassembled_unmapped.fasta mouse1_unassembled.dmdout \
mouse1_unassembled_unannotated.fasta
