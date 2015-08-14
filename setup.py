import os
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

VERSION="0.4.2"

setup(
    name='knead-datalib',
    version=VERSION,
    description='',
    packages=['knead_datalib'],
    zip_safe=False,
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
    scripts=[
        os.path.join(here, "knead_data.py"), 
        os.path.join(here, "download_db.py"),
        os.path.join(here, "generate_db.py")
        ],
    entry_points = {
        "distutils.commands": [
            "trimmomatic = knead_datalib.util:DownloadTrimmomaticCommand"
        ],
        "console_scripts": [
            "kneaddata = knead_data:main",
            "kneaddata_database = download_db:main"
        ]
    }
)

