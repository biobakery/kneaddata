
import os

data_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")

fastq_file=os.path.join(data_folder,"demo.fastq")
bowtie2_db_folder=os.path.join(data_folder,"demo_bowtie2_db")

# kneaddata default file names
log_extension="_kneaddata.log"
single_trim_extension="_kneaddata.trimmed.fastq"
clean_extension="_clean.fastq"
contaminated_extension="_contam.fastq"
sam_extension=".sam"
final_extension="_kneaddata.fastq"
final_extensions_paired=["_kneaddata_pe_1.fastq","_kneaddata_pe_2.fastq"]
name_delimiter="_kneaddata_"
paired_trim_extensions=["_kneaddata.trimmed.1.fastq","_kneaddata.trimmed.2.fastq"]
unjoined_trim_extensions=["_kneaddata.trimmed.single.1.fastq","_kneaddata.trimmed.single.2.fastq"]
