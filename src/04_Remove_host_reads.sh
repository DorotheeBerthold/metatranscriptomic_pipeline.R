### Remove host reads

# Generate database for host (mouse DNA) sequences
###########################################################################################

# Optional: DB from 04 & 05 could be combined with cat UniVec_Core mouse_cds.fa > contaminants.fa
# --> however then hard to distinguish how many reads came from host

# Here, we use Ensembl as DB for mouse DNA sequences
wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz
gzip -d Mus_musculus.GRCm39.cds.all.fa.gz

# Rename file
mv Mus_musculus.GRCm39.cds.all.fa mouse_cds.fa

# Generate index with BWA
bwa index -a bwtsw mouse_cds.fa
samtools faidx mouse_cds.fa
makeblastdb -in mouse_cds.fa -dbtype nucl

# Align sequence with database using BWA
###########################################################################################

bwa mem -t 4 mouse_cds.fa mouse1_univec_blat.fastq > mouse1_mouse_bwa.sam
samtools view -bS mouse1_mouse_bwa.sam > mouse1_mouse_bwa.bam

# Generate the two output files
samtools fastq -n -F 4 -0 mouse1_mouse_bwa_contaminats.fastq mouse1_mouse_bwa.bam
samtools fastq -n -f 4 -0 mouse1_mouse_bwa.fastq mouse1_mouse_bwa.bam

# Additional alignment with database using BLAT
###########################################################################################

# Convert to FASTA
vsearch --fastq_filter mouse1_mouse_bwa.fastq --fastaout mouse1_mouse_bwa.fasta

# Search with BLAT
blat -noHead -minIdentity=90 -minScore=65  mouse_cds.fa mouse1_mouse_bwa.fasta -fine \
-q=rna -t=dna -out=blast8 mouse1_mouse.blatout

# Run Python script
./1_BLAT_Filter_DB.py mouse1_mouse_bwa.fastq mouse1_mouse.blatout \
mouse1_mouse_blat.fastq mouse1_mouse_blat_contaminats.fastq
