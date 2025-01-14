from Bio import SeqIO

def convert_blat_to_sam(blat_file, fasta_file, output_sam):
    # Read sequences from the FASTA file (query sequences)
    fasta_sequences = {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}
    
    # Open BLAT output and SAM file for writing
    with open(blat_file, 'r') as blast, open(output_sam, 'w') as sam:
        # Write SAM header
        sam.write("@HD\tVN:1.0\tSO:unsorted\n")
        
        subject_names = set()
        for line in blast:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            qname, tname = cols[0], cols[1]
            subject_names.add(tname)
        
        # Write @SQ header lines for each subject (transcript IDs from BLAT output)
        for ref_name in subject_names:
            ref_len = 1000  # Placeholder length, adjust as needed
            sam.write(f"@SQ\tSN:{ref_name}\tLN:{ref_len}\n")
        
        blast.seek(0)
        for line in blast:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            qname, tname, identity, length, mismatch, gap, qstart, qend, tstart, tend, evalue, bitscore = cols
            
            # Assuming all reads are mapped with flag 0 (this can be adjusted as needed)
            flag = 0  # Adjust based on alignment status
            cigar = f"{length}M"  # Simple alignment: no gaps or mismatches
            mapq = 255  # Placeholder for mapping quality
            seq = fasta_sequences.get(qname, "*")  # Retrieve sequence for query read
            qual = "*"  # Placeholder for quality score
            
            # Write the SAM alignment record
            sam.write(f"{qname}\t{flag}\t{tname}\t{tstart}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")
