### Functional annotation of enzymes

# Match annotated genes to enzymes in the KEGG pathways
###########################################################################################

# PACKAGE: DIAMOND

# Annotation of microbial genes (identified with BWA & BLAT)
# database with older version of diamond and is incompatible (rerun 09)
diamond blastx -p 4 -d swiss_db -q mouse1_genes.fasta -o mouse1_genes.diamondout -f 6 -t dmnd_tmp -e 10 -k 1

# Annotation of proteins identified with DIAMOND
diamond blastp -p 4 -d swiss_db -q mouse1_proteins.fasta -o mouse1_proteins.diamondout -f 6 -t dmnd_tmp -e 10 -k 1

# Run small python script to generate mapping file listing our gene/protein & 
# enzyme commission (EC)
###########################################################################################

# Activate permission
chmod +x 7_Gene_EC_Map_DB.py

# Command input: 7_Gene_EC_Map_DB.py <SWISS-PROT_EC_Mappings> <Diamond_Output_For_Genes> \
# <Diamond_Output_For_Proteins> <Output_EC_Mapping_File>

./7_Gene_EC_Map_DB.py swiss_map.tsv mouse1_genes.diamondout \
mouse1_proteins.diamondout mouse1_EC_map.tsv
