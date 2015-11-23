# Settings for the knead_data script. 

# File endings for BMTagger's required database files
bowtie2_db_endings =    [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", 
                        ".rev.2.bt2"]
bmtagger_db_endings =   [".bitmask", ".srprism.amp", 
                        ".srprism.idx", ".srprism.imp",
                        ".srprism.map", ".srprism.pmp", 
                        ".srprism.rmp", ".srprism.ss",
                        ".srprism.ssa", ".srprism.ssd", 
                        ".nhr", ".nin", ".nsq"]
    

# Trimmomatic file endings for single end and paired end, respectively
trimomatic_se_ending = ".trimmed.fastq"

trimomatic_pe_endings =   [".trimmed.1.fastq", 
                    ".trimmed.2.fastq", 
                    ".trimmed.single.1.fastq", 
                    ".trimmed.single.2.fastq"]

# BMTagger file endings if you choose to remove the contaminant reads. For
# single end and paired end reads, respectively.
bmtagger_se_ending = ".fastq"
bmtagger_pe_endings =   ["_1.fastq",
                        "_2.fastq"]
