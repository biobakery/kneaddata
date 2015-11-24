# Settings for the knead_data script. 

# Default settings for command line arguments
threads=1

strategy_choices=["memory","storage"]
strategy=strategy_choices[1]

log_level_choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
log_level=log_level_choices[0]

bmtagger_exe="bmtagger.sh"

trimmomatic_jar="trimmomatic-0.33.jar"
trimmomatic_memory="500m"
trimmomatic_options=["SLIDINGWINDOW:4:20", "MINLEN:60"]

bowtie2_exe="bowtie2"
bowtie2_options=["--very-sensitive"]

trf_exe="trf"
trf_match=2
trf_mismatch=7
trf_delta=7
trf_match_probability=80
trf_pi=10
trf_minscore=50
trf_maxperiod=500

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
