import os
import sys
import logging
from math import floor
from multiprocessing import cpu_count

def divvy_threads(args):
    avail_cpus = args.threads or cpu_count()-1
    n_consumers = len(args.reference_db)
    trim_threads = 1
    align_threads = max(1, floor(avail_cpus/float(n_consumers)))
    return int(trim_threads), int(align_threads)
    

def try_create_dir(d):
    if not os.path.exists(d):
        logging.warning("Directory `%s' doesn't exist. Creating.", d)
        try:
            os.makedirs(d)
        except Exception as e:
            logging.crit("Unable to create directory `%s': %s", d, str(e))
            sys.exit(1)
