#!/usr/bin/env python

"""
download_db.py

Helper script to download databases for the knead_data.py pipeline
Requires an internet connection and Python v2.7+.
"""

import sys
    
# Try to load argparse to check for python v2.7+
try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to Python version v2.7")
    
import os
import tarfile
import time

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

# the locations of the current databases to download
current_downloads={
    # genome/database
    "human_genome" : {
        # database build type
        "bowtie2" : "http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg37_and_human_contamination_Bowtie2_v0.1.tar.gz",
        "bmtagger" : "http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_BMTagger_v0.1.tar.gz"
    },
    "human_transcriptome" : {
        "bowtie2" : "http://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz"
    },
    "ribosomal_RNA" : {
        "bowtie2" : "http://huttenhower.sph.harvard.edu/kneadData_databases/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2.tar.gz"
    },
    "mouse_C57BL" : {
        "bowtie2" : "http://huttenhower.sph.harvard.edu/kneadData_databases/mouse_C57BL_6NJ_Bowtie2_v0.1.tar.gz"
    }
}

class ReportHook():
    def __init__(self):
        self.start_time=time.time()
        
    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        
        if blocknum == 0:
            self.start_time=time.time()
            if total_size > 0:
                print("Downloading file of size: " + "{:.2f}".format(total_size / (1024.0**3)) + " GB\n")
        else:
            total_downloaded=blocknum*block_size
            status = "{:3.2f} GB ".format(total_downloaded / (1024.0**3))
                    
            if total_size > 0:
                percent_downloaded=total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                download_rate=total_downloaded/(time.time()-self.start_time)
                estimated_time=(total_size-total_downloaded)/download_rate
                estimated_minutes=int(estimated_time/60.0)
                estimated_seconds=estimated_time-estimated_minutes*60.0
                status +="{:3.2f}".format(percent_downloaded) + " %  " + \
                    "{:5.2f}".format(download_rate / (1024.0**2)) + " MB/sec " + \
                    "{:2.0f}".format(estimated_minutes) + " min " + \
                    "{:2.0f}".format(estimated_seconds) + " sec "
            status+="        \r"
            sys.stdout.write(status)
            

def download_tar_and_extract_with_progress_messages(url, filename, folder):
    """
    Download the file at the url
    """
    # check for local file
    local_file = False
    if os.path.isfile(url):
        local_file = True   
 
    if not local_file:
        print("Download URL: " + url) 

    try:
        if not local_file:
            url_handle = urlretrieve(url, filename, reporthook=ReportHook().report)
        else:
            filename = url    
        print("\nExtracting: " + filename)
        tarfile_handle=tarfile.open(filename)
        tarfile_handle.extractall(path=folder)
    except (EnvironmentError, tarfile.ReadError):
        if local_file:
            sys.exit("CRITICAL ERROR: Unable to extract from local file: " + url)
        else:
            sys.exit("CRITICAL ERROR: Unable to download and extract from URL: " + url)

def check_user_database(original, user):
    """ Check that the user database is of the expected version """

    if original.split("/")[-1] == user.split("/")[-1]:
        return True
    else:
        sys.exit("The user database selected does not match that expected: "+original.split("/")[-1])

def download_database(database, build, location, database_location):
    """
    Download and decompress the selected database
    """
    
    if database in current_downloads:
        if build in current_downloads[database]:
            # download the database
            downloaded_file=os.path.join(location,current_downloads[database][build].split('/')[-1])
            if database_location:
                check_user_database(current_downloads[database][build],database_location)
                download_tar_and_extract_with_progress_messages(database_location,
                    downloaded_file, location)
            else:
                download_tar_and_extract_with_progress_messages(current_downloads[database][build], 
                    downloaded_file, location)
            
            # remove the download (if not a local file provided by the user)
            if not database_location or (database_location and not os.path.isfile(database_location)):
                try:
                    os.unlink(downloaded_file)
                except EnvironmentError:
                    print("Unable to remove file: " + downloaded_file)
            
            print("Database installed: " + location + "\n")
        else:
            sys.exit("ERROR: Please select an available build.")
    else:
        sys.exit("ERROR: Please select an available database.")
        
    return location

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "KneadData Databases\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--available", 
        action="store_true",
        help="print the available databases\n")
    parser.add_argument(
        "--download", 
        nargs=3,
        metavar=("<database>","<build>","<install_location>"),
        help="download the selected database to the install location\n")
    parser.add_argument(
        "--database-location",
        help="location (local or remote) to pull the database")
    
    return parser.parse_args()

def main():
    # Parse arguments from the command line
    args=parse_arguments(sys.argv)
    
    if args.download:
        # download the database
        database=args.download[0]
        build=args.download[1]
        location=os.path.abspath(args.download[2])
        
        # create the install location if it does not already exist
        if not os.path.isdir(location):
            try:
                print("Creating directory to install database: " + location)
                os.mkdir(location)
            except EnvironmentError:
                sys.exit("CRITICAL ERROR: Unable to create directory: " + location)
        
        install_location=download_database(database,build,location,args.database_location)
    
    if args.available or not args.download:
        # print the available databases
        print("KneadData Databases ( database : build = location )")
        for database in current_downloads:
            for build, location in current_downloads[database].items():
                print(database+" : "+build+" = "+location)
                
if __name__ == '__main__':
    main()

