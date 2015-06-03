import os
import sys
import urllib
from zipfile import ZipFile
from setuptools import Command

here = os.path.abspath(os.path.dirname(__file__))


class DownloadTrimmomaticCommand(Command):
    description = "Download and unpack Trimmomatic"
    user_options    = [ ('to=', 't', "Download databases to this directory") ]


    download_config = {
        "trimmomatic.zip": 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip'
    }

    def initialize_options(self):
        self.to = None
        self.download_dir = os.path.join(here, '..')


    def finalize_options(self):
        if self.to:
            self.download_dir = os.path.abspath(self.to)

        if not os.path.exists(self.download_dir):
            os.mkdir(self.download_dir)
                

    def download(self, name):
        url = self.download_config[name]
        fname = os.path.join(self.download_dir, name)

        # report every 10 MB worth of reads or 1280 blocks, since the
        # default block size for urllib is 8KB:
        # https://github.com/python/cpython/blob/master/Lib/urllib/request.py#L204
        report_interval = (1024*1024*10) / (1024*8)
        def _reporthook(blocknum, bs, size):
            if blocknum % report_interval == 0:
                percent_complete = (100 * blocknum * bs) / float(size)
                msg = "\rDownloading. %.1f%% complete"
                sys.stderr.write( msg%(percent_complete) )

        print "Downloading "+name
        urllib.urlretrieve(url, fname, _reporthook)
        sys.stderr.write("\rDownload complete.\n")

        with ZipFile(fname, 'r') as thezip:
            print >> sys.stderr, "Extracting "+name
            thezip.extractall(self.download_dir)
        print >> sys.stderr, "Extract complete"
        
        os.unlink(name)

    def run(self):
        self.download("trimmomatic.zip")
        print >> sys.stderr, "Thanks!"


    help_options    = [ ]
    boolean_options = [ ]
    negative_opt    = { }
    default_format  = { }

