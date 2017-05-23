
import sys
import os

# required python versions (2.7+ or 3.0+)
required_python_version_major = [2,3]
required_python_version_minor = [7,0]

# check for either of the required versions
pass_check=False
try:
    for major, minor in zip(required_python_version_major, required_python_version_minor):
        if (sys.version_info[0] == major and sys.version_info[1] >= minor):
            pass_check=True
except (AttributeError,IndexError):
    sys.exit("CRITICAL ERROR: The python version found (version 1) " +
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")

if not pass_check:
    sys.exit("CRITICAL ERROR: The python version found (version "+
        str(sys.version_info[0])+"."+str(sys.version_info[1])+") "+
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")


try:
    import setuptools
except ImportError:
    sys.exit("Please install setuptools.")
    
import tarfile
import shutil
import zipfile
import tempfile
import re
import time

from setuptools.command.install import install as _install

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

VERSION="0.6.0"
AUTHOR = "KneadData Development Team"
AUTHOR_EMAIL = "kneaddata-users@googlegroups.com"

COUNTER_URL="http://bitbucket.org/biobakery/kneaddata/downloads/counter.txt"

setup_directory = os.path.abspath(os.path.dirname(__file__))

def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """
    
    return byte / (1024.0**2)

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
                print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
        else:
            total_downloaded=blocknum*block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))
                    
            if total_size > 0:
                percent_downloaded=total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                download_rate=total_downloaded/(time.time()-self.start_time)
                estimated_time=(total_size-total_downloaded)/download_rate
                estimated_minutes=int(estimated_time/60.0)
                estimated_seconds=estimated_time-estimated_minutes*60.0
                status +="{:3.2f}".format(percent_downloaded) + " %  " + \
                    "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
                    "{:2.0f}".format(estimated_minutes) + " min " + \
                    "{:2.0f}".format(estimated_seconds) + " sec "
            status+="        \r"
            sys.stdout.write(status)

def download(url, download_file):
    """
    Download a file from a url
    """

    try:
        print("Downloading "+url)
        file, headers = urlretrieve(url,download_file,reporthook=ReportHook().report)
    except EnvironmentError:
        print("Warning: Unable to download "+url)
    

def download_unpack_tar(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("Warning: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file)
    
    error_during_extract=False
    
    try:
        tarfile_handle=tarfile.open(download_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except EnvironmentError:
        print("Warning: Unable to extract "+software_name+".")
        error_during_extract=True
        
    if not error_during_extract:
        try:
            os.unlink(download_file)
        except EnvironmentError:
            print("Warning: Unable to remove the temp download: " + download_file)
        
def download_unpack_zip(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("Warning: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file)
    
    error_during_extract=False
    
    try:
        zipfile_handle=zipfile.ZipFile(download_file)
        zipfile_handle.extractall(path=folder)
        zipfile_handle.close()
    except EnvironmentError:
        print("Warning: Unable to extract "+software_name+".")
        error_during_extract=True
        
    if not error_during_extract:
        try:
            os.unlink(download_file)
        except EnvironmentError:
            print("Warning: Unable to remove the temp download: " + download_file)
            
def find_exe_in_path(exe, bypass_permissions_check=None):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if bypass_permissions_check or os.access(fullexe,os.X_OK):
                return path
    return None
            
def install_bowtie2(final_install_folder, mac_os, replace_install=None):
    """ 
    Download and install the bowtie2 software if not already installed
    """
    
    # Check if bowtie2 is already installed
    bowtie2_installed=find_exe_in_path("bowtie2")
    
    if not bowtie2_installed or replace_install:
        bowtie2_exe="bowtie2"
        bowtie2_file="bowtie2-2.2.3-linux-x86_64.zip"
        bowtie2_url="http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip"

        # if this is a MAC OS, select a different binary download
        if mac_os:
            bowtie2_file="bowtie2-2.2.3-macos-x86_64.zip"
            bowtie2_url="http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-macos-x86_64.zip"
            
        bowtie2_folder="bowtie2-2.2.3"
    
        tempfolder=tempfile.mkdtemp(prefix="bowtie2_download_",dir=setup_directory)

        # install the bowtie2 software
        print("Installing bowtie2.")
        error_during_install=False
        download_unpack_zip(bowtie2_url, bowtie2_file, tempfolder, bowtie2_exe)
        
        # copy the installed software to the final bin location
        # copy all bowtie2* executables
        fullpath_bowtie2_exe=os.path.join(tempfolder, bowtie2_folder)
        
        files=[]
        try:
            files=os.listdir(fullpath_bowtie2_exe)
        except EnvironmentError:
            print("Warning: Bowtie2 files not found.")
            error_during_install=True
        
        for file in files:  
            # check if this file is one of the bowtie2* executables      
            if re.match(bowtie2_exe,file):  
                try:   
                    # copy to the install folder
                    shutil.copy(os.path.join(fullpath_bowtie2_exe,file), final_install_folder)
                    # add executable permissions
                    os.chmod(os.path.join(final_install_folder,file), 0o755)
                except (EnvironmentError, shutil.Error):
                    error_during_install=True
            
        # remove the local bowtie2 install
        try:
            shutil.rmtree(tempfolder)
        except EnvironmentError:
            print("Warning: Unable to remove temp install folder.")
        
        if error_during_install:
            print("Warning: Unable to install bowtie2. Please install bowtie2.")
        else:
            print("Installed bowtie2 at "+final_install_folder)
    else:
        print("Found bowtie2 install at "+bowtie2_installed)
        
def install_trimmomatic(final_install_folder, mac_os, replace_install=None):
    """ Download and install Trimmomatic if not already installed
    """
    
    trimmomatic_jar="trimmomatic-0.33.jar"
    
    # check if trimmomatic is already installed
    trimmomatic_installed=find_exe_in_path(trimmomatic_jar, bypass_permissions_check=True)

    if not trimmomatic_installed or replace_install:
        trimmomatic_file="Trimmomatic-0.33.zip"
        trimmomatic_url="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip"
        
        tempfolder=tempfile.mkdtemp(prefix="trimmomatic_download_",dir=setup_directory)

        # install the trimmomatic software
        print("Installing Trimmomatic.")
        error_during_install=False
        download_unpack_zip(trimmomatic_url, trimmomatic_file, tempfolder, trimmomatic_jar)
        
        trimmomatic_jar_full_path=os.path.join(tempfolder, "Trimmomatic-0.33", trimmomatic_jar)
            
        # copy the installed software to the final bin location
        try:
            # copy to the install folder
            shutil.copy(trimmomatic_jar_full_path, final_install_folder)
            # add executable permissions
            os.chmod(os.path.join(final_install_folder,trimmomatic_jar), 0o755)
        except (EnvironmentError, shutil.Error):
            error_during_install=True
            
        # remove the local Trimmomatic install
        try:
            shutil.rmtree(tempfolder)
        except EnvironmentError:
            print("Warning: Unable to remove temp install folder.")
        
        if error_during_install:
            print("Warning: Unable to install Trimmomatic. Please install Trimmomatic.")
        else:
            print("Installed Trimmomatic at "+final_install_folder)
        
    else:
        print("Found Trimmomatic install at "+trimmomatic_installed)

class Install(_install):
    """
    Custom setuptools install command
    """
    
    _install.user_options=_install.user_options+[('bypass-dependencies-install', 
        None, 'bypass install of dependencies')]
    
    def initialize_options(self):
        self.bypass_dependencies_install=False
        _install.initialize_options(self)
    
    def finalize_options(self):
        _install.finalize_options(self)
    
    def run(self):
        # try to download the bitbucket counter file to count downloads
        # this has been added since PyPI has turned off the download stats
        # this will be removed when PyPI Warehouse is production as it
        # will have download stats
        counter_file="counter.txt"
        if not os.path.isfile(counter_file):
            print("Downloading counter file to track kneadddata downloads"+
                  " since the global PyPI download stats are currently turned off.")
            download(COUNTER_URL,counter_file)
        _install.run(self)
                    
        # find out the platform
        mac_os=False
        if sys.platform in ["darwin","os2","os2emx"]:
            mac_os=True
        
        # install dependencies if not already installed
        if not self.bypass_dependencies_install:
            install_trimmomatic(self.install_scripts,mac_os,replace_install=False)
            install_bowtie2(self.install_scripts,mac_os,replace_install=False)
        else:
            print("Bypass install of dependencies.")

setuptools.setup(
    name='kneaddata',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    long_description="KneadData is a tool designed to perform quality control on metagenomic " + \
        "sequencing data, especially data from microbiome experiments. In these experiments, " + \
        "samples are typically taken from a host in hopes of learning something about the " + \
        "microbial community on the host. However, metagenomic sequencing data from such " + \
        "experiments will often contain a high ratio of host to bacterial reads. This " + \
        "tool aims to perform principled in silico separation of bacterial reads from " + \
        "these \"contaminant\" reads, be they from the host, from bacterial 16S " + \
        "sequences, or other user-defined sources.",
    url="http://huttenhower.sph.harvard.edu/kneaddata",
    keywords=['microbial','microbiome','bioinformatics','microbiology','metagenomic','metatranscriptomic','kneaddata'],
    platforms=['Linux','MacOS'],
    packages=setuptools.find_packages(),
    package_data={
        'kneaddata' : [
            'tests/data/*.*',
            'tests/data/demo_bowtie2_db/*'
        ]},
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    cmdclass={'install': Install},
    entry_points = {
        "console_scripts": [
            "kneaddata = kneaddata.knead_data:main",
            "kneaddata_bowtie2_discordant_pairs = kneaddata.bowtie2_discordant_pairs:main",
            "kneaddata_database = kneaddata.download_db:main",
            "kneaddata_build_database = kneaddata.generate_db:main",
            "kneaddata_read_count_table = kneaddata.read_count_table:main",
            "kneaddata_test = kneaddata.tests.kneaddata_test:main"
        ]
    }
)

