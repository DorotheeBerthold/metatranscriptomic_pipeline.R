### Process the reads

# Trim adaptor sequences & low quality reads with sliding window
###########################################################################################
# PACAKGE: Trimmomatic
# Old filename: mouse1.fastq
# New filename: mouse1_trim.fastq
# Remove adaptors
# Trim bases below quality score of 3 at beginning & end
# Slide with window of size 4 with local quality of 15 & trim if quality is below
# Delete sequences with length <50

java -jar /Users/dobertho/anaconda3/share/trimmomatic-0.39-2/trimmomatic.jar SE \
mouse1.fastq \
mouse1_trim.fastq \
ILLUMINACLIP:Adapters:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:50

# Check updated QC
fastqc mouse1_trim.fastq

# Add another quality treshold - this time not only with sliding window but overall quality
###########################################################################################
# PACKAGE: vsearch
# fastq_filter: use quality filter algorithm
# fastq_maxee: expected error treshold (set at 1)
# fastq_out output file name
vsearch --fastq_filter mouse1_trim.fastq --fastq_maxee 2.0 --fastqout mouse1_qual.fastq

# Check updated QC
fastqc mouse1_qual.fastq
