import os
import sys
import shlex
import logging
import tempfile
from math import floor
from functools import partial
from contextlib import contextmanager
from multiprocessing import cpu_count

def divvy_threads(args):
    avail_cpus = args.threads or cpu_count()-1
    n_consumers = len(args.reference_db)
    trim_threads = 1
    if n_consumers > 0:
        align_threads = max(1, floor(avail_cpus/float(n_consumers)))
    else:
        align_threads = 1
    return int(trim_threads), int(align_threads)
    

def try_create_dir(d):
    if not os.path.exists(d):
        logging.warning("Directory `%s' doesn't exist. Creating.", d)
        try:
            os.makedirs(d)
        except Exception as e:
            logging.crit("Unable to create directory `%s': %s", d, str(e))
            sys.exit(1)


@contextmanager
def mktempfifo(names=("a",)):
    tmpdir = tempfile.mkdtemp()
    names = map(partial(os.path.join, tmpdir), names)
    map(os.mkfifo, names)
    try:
        yield names
    finally:
        # still perform cleanup even if there were exceptions/errors in the
        # "with" block
        map(os.remove, names)
        os.rmdir(tmpdir)


@contextmanager
def mkfifo_here(names=("a",), mode=0600):
    for n in names:
        os.mkfifo(n, mode)
    try:
        yield names
    finally:
        for n in names:
            os.remove(n)


def process_return(name, retcode, stdout, stderr):
    if name:
        logging.debug("Finished running %s!" %name)
    if retcode:
        log = logging.critical
        log("%s exited with exit status %d", name, retcode)
    else:
        log = logging.debug
    if stdout:
        log("%s stdout:\n%s", name, stdout)
    if stderr:
        log("%s stderr:\n%s", name, stderr)
    if retcode:
        sys.exit(retcode)


def parse_positive_int(string):
    try:
        val = int(string)
    except ValueError:
        raise argparse.ArgumentTypeError("Unable to parse %s to int" %string) 
    if val <= 0:
        raise argparse.ArgumentTypeError("%s is not a positive integer" %string)
    return val


def _get_bowtie2_args(bowtie2_args):
    for arg in map(shlex.split, bowtie2_args):
        for a in arg:
            yield a
