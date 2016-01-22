import os
import subprocess
import shutil

import cfg

def run_kneaddata(command):
    """ Run the kneaddata command """
    
    try:
        stdout=subprocess.check_output(command)
    except (EnvironmentError, subprocess.CalledProcessError):
        print("Warning: Unable to execute kneaddata in test.\n"+" ".join(command))    
        
def remove_temp_folder(tempdir):
    """ Remove the temp folder """
    
    try:
        shutil.rmtree(tempdir)
    except EnvironmentError:
        print("Warning: Unable to remove temp directory: " + tempdir)

def remove_temp_file(file):
    """ Try to remove the temp file """
    
    try:
        os.remove(file)
    except EnvironmentError:
        sys.exit("Warning: Unable to remove temp file" + file) 
        
def check_output(output_files_expected,output_folder):
    """ Check the output folder has the expected file and they are all non-zero """
    
    for file in output_files_expected:
        expected_file = os.path.join(output_folder,file)
        # check the file exists
        yield (os.path.isfile(os.path.join(expected_file)), "File does not exist: " + file)
        
        # check the file is not empty
        yield (os.stat(expected_file).st_size > 0, "File is empty: " + file)
    
def file_basename(file):
    """ Get the basename for the file without the extension """
    
    return os.path.splitext(os.path.basename(file))[0]

def database_basename(folder):
    """ Get the basename for the database in the folder """
    
    # get a single database file from the folder
    database_file=os.listdir(folder)[0]
    
    # some files have multiple extensions, so use split instead of splitext
    return database_file.split(".")[0]

def get_filtered_file_basename(basename,db_folder,software, paired=None):
    """ Get the basename of the files created by bowtie2 and bmtagger """
    
    database=database_basename(db_folder)
    
    if paired:
        name_delimiter=cfg.name_delimiter_paired
    else:
        name_delimiter=cfg.name_delimiter
    
    return basename+name_delimiter+database+"_"+software

    