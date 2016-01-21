
import os

data_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")

fastq_file=os.path.join(data_folder,"demo.fastq")
bowtie2_db_folder=os.path.join(data_folder,"demo_bowtie2_db")

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
final_extension="_kneaddata.fastq"
name_delimiter="_kneaddata_"
name_delimiter_paired=name_delimiter+"paired_"
final_extensions_paired=[name_delimiter_paired+"1.fastq",name_delimiter_paired+"2.fastq"]
paired_trim_extensions=["_kneaddata.trimmed.1.fastq","_kneaddata.trimmed.2.fastq"]
unjoined_trim_extensions=["_kneaddata.trimmed.single.1.fastq","_kneaddata.trimmed.single.2.fastq"]
