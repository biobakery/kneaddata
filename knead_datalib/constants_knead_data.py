# Constants for the knead_data script. 

# File endings for BMTagger's required database files
BOWTIE2_DB_ENDINGS =    [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", 
                        ".rev.2.bt2"]
BMTAGGER_DB_ENDINGS =   [".bitmask", ".srprism.amp", 
                        ".srprism.idx", ".srprism.imp",
                        ".srprism.map", ".srprism.pmp", 
                        ".srprism.rmp", ".srprism.ss",
                        ".srprism.ssa", ".srprism.ssd", 
                        ".nhr", ".nin", ".nsq"]
    

# Trimmomatic file endings for single end and paired end, respectively
TRIM_SE_ENDING = ".trimmed.fastq"

TRIM_PE_ENDINGS =   [".trimmed.1.fastq", 
                    ".trimmed.2.fastq", 
                    ".trimmed.single.1.fastq", 
                    ".trimmed.single.2.fastq"]

# BMTagger file endings if you choose to remove the contaminant reads. For
# single end and paired end reads, respectively.
BMTAGGER_SE_ENDING = ".fastq"
BMTAGGER_PE_ENDINGS =   ["_1.fastq",
                        "_2.fastq"]
BMTAGGER_OUTPUT_ENDINGS = ["_pe.out", "_se_1.out", "_se_2.out"]
