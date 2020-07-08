# KneadData Databases #


#### Bowtie2 Databases

## [human_genome](http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz) 
Version: 0.1
Source of fasta: 


## [human_transcriptome](http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz) 
Version: 0.1
Source of fasta: 


## [mouse_C57BL](http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz)
Version: 0.1
Source of fasta: 

Steps to generate **human_genome,human_transcriptome, and mouse_C57BL bowtie2 database**:
```
# These are the commands run to create the databases and the synthetic reads for human_genome,human_transcriptome, and mouse_C57BL bowtie2 database.

$ bowtie2-inspect /genome_demo | head -n 1000 > genome_demo.fasta

$ ./art_illumina -i genome_demo.fasta -l 100 -o genome_demo_synthetic_reads -sam -ss HS20 -f 0.3 --id genome_demo

$ bowtie2-build genome_demo.fasta genome_demo
```


## [ribosomal_RNA](http://huttenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2.tar.gz) 
Version: 0.2
Source of fasta: https://www.arb-silva.de/no_cache/download/archive/release_128/Exports/

Steps to generate ribosomal_RNA bowtie2 database:
```
# These are the commands to run to create the databases and the synthetic reads.

$ python -u modify_RNA_to_DNA.py ribosomal_RNA_initial.fasta  ribosomal_RNA.fasta

$ bowtie2-inspect /ribosomal_RNA | head -n 1000 > ribosomal_RNA.fasta

$ ./art_illumina -i ribosomal_RNA.fasta -l 100 -o ribosomal_RNA_synthetic_reads -sam -ss HS20 -f 0.3 --id ribosomal_RNA

$ bowtie2-build ribosomal_RNA.fasta ribosomal_RNA
```