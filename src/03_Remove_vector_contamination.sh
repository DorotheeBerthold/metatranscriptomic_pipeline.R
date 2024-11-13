### Remove vector contamination

# Generate database for vector, links, adapters & primers
###########################################################################################

# Get NCBI Univec Database
wget ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core

# PACKAGE: BWA
# Generate index for sequences from databases in FASTA format
#-a: algorithm for construction of index (either is, bwtsw or rb2)
bwa index -a bwtsw UniVec_Core

# Make new BLAST DB
#samtools faidx: uses input reference fastq
samtools faidx UniVec_Core
makeblastdb -in UniVec_Core -dbtype nucl

# Align sequence with database using BWA
###########################################################################################

# Align with BWA
# bwa mem: Generates alignments of reads to database (maximal exact matches)
# -t: numer of threads
# Here, we use output from 02 instead of optional step 03
# generates .sam output
bwa mem -t 4 UniVec_Core mouse1_qual.fastq > mouse1_univec_bwa.sam

# Convert .sam into .bam files
samtools view -bS mouse1_univec_bwa.sam > mouse1_univec_bwa.bam

# Generate fastq file for all reads mapping to vector DB (-F 4)
samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam

# Generate fastq file for all reads that did NOT map to vector DB (-f 4)
samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam

# Additional alignment with database using BLAT (often finds reads that BWA misses)
###########################################################################################

# PACKAGE: BLAT
# Only accepts files in fasta format
# --fastq_filter without criteria just uses everything for output
vsearch --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta


# -noHead: suppresses .psl header (so it becomes a tab-separated file)
# -minIdentity: sets minimal identity with DB
# -minScore: Sets minimum score to 65 - matches minus mismatches - gap penalty
# -fine: for high-quality mRNAs
# -q: Query type is RNA sequence
# -t: DB is DNA sequence
blat -noHead -minIdentity=90 -minScore=65  UniVec_Core mouse1_univec_bwa.fasta -fine \
-q=rna -t=dna -out=blast8 mouse1_univec.blatout

# Run small python script to filter reads that BLAT does not confidently align to DB
###########################################################################################

# Command input: 1_BLAT_Filter.py <Input_Reads.fq> <BLAT_Output_File> 
# <Unmapped_Reads_output> <Mapped_Reads_Output>
./1_BLAT_Filter_DB.py mouse1_univec_bwa.fastq mouse1_univec.blatout \
mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq
