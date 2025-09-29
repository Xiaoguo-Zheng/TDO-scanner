##Download reference  
#Mus_musculus reference  
https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz  

#Homo_sapiens reference  
https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  

#Mus_musculus gtf  
https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.chr.gtf.gz  

#Homo_sapiens gtf  
https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gtf.gz
  
##install gffread  
conda install bioconda::gffread
  
##Extract transcripts fasta  
#mouse  
gffread -w 0.mm39__matureRNA_seq.fa -g Mus_musculus.GRCm39.dna.primary_assembly.fa -F -W Mus_musculus.GRCm39.113.chr.gtf  
  
#human  
gffread -w 0.hg38_matureRNA_seq.fa -g Homo_sapiens.GRCh38.dna.primary_assembly.fa -F -W Homo_sapiens.GRCh38.113.chr.gtf  
  
##get all candidate prePAM_20bp from transcripts  
perl get_all_candidate_PAM_from_reference.pl Homo_sapiens.GRCh38.dna.primary_assembly.fa > Homo_sapiens.GRCh38_all_candidate_PAM.txt  
