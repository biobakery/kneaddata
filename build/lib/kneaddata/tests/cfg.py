
import os

data_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")

fastq_file=os.path.join(data_folder,"demo.fastq")
fastq_pair_file=os.path.join(data_folder,"demo2.fastq")
fastq_file_gzipped=os.path.join(data_folder,"demo.fastq.gz")
file_sam=os.path.join(data_folder,"demo.sam")
file_bam=os.path.join(data_folder,"demo.bam")
fastq_file_matches_sam_and_bam=os.path.join(data_folder,"demo_matches_sam_and_bam.fastq")
bowtie2_db_folder=os.path.join(data_folder,"demo_bowtie2_db")
bowtie2_db_file=os.path.join(bowtie2_db_folder,"demo_db.1.bt2")
bowtie2_db_index=os.path.join(bowtie2_db_folder,"demo_db")

merge_files=[os.path.join(data_folder,"merge1.fastq"),
             os.path.join(data_folder,"merge2.fastq"),
             os.path.join(data_folder,"merge3.fastq")]

merge_files_sequences_intersect=["all_files_1","all_files_2"]
merge_files_1_sequences=["all_files_1","file1_1","file1_and_file2_1","all_files_2"]

# kneaddata default file names
log_extension="_kneaddata.log"
single_trim_extension="_kneaddata.trimmed.fastq"

clean_extension="_clean.fastq"
contaminated_extension="_contam.fastq"
paired_contaminated_extension=["_contam_1.fastq","_contam_2.fastq"]

sam_extension=".sam"

name_delimiter="_kneaddata_"
name_delimiter_paired=name_delimiter+"paired_"
name_delimiter_orphan=name_delimiter+"unmatched_"

final_extension="_kneaddata.fastq"
final_extensions_paired=[name_delimiter_paired+"1.fastq",name_delimiter_paired+"2.fastq"]
final_extensions_discordant=[name_delimiter_orphan+"1.fastq",name_delimiter_orphan+"2.fastq"]

paired_trim_extensions=["_kneaddata.trimmed.1.fastq","_kneaddata.trimmed.2.fastq"]
unjoined_trim_extensions=["_kneaddata.trimmed.single.1.fastq","_kneaddata.trimmed.single.2.fastq"]

repeats_removed_extension="_kneaddata.repeats.removed.fastq"
paired_repeats_removed_extensions=["_kneaddata.repeats.removed.1.fastq","_kneaddata.repeats.removed.2.fastq"]

fastqc_extensions=["_fastqc.html","_fastqc.zip"]
