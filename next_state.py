import argparse
from glob import glob
import shutil
import random

parser = argparse.ArgumentParser()
parser.add_argument('resultsdir')
args = parser.parse_args()

resultsdir = args.resultsdir

seeds = glob('/home/server/server2/projects/17805/SEEDS/*.xml.bz2')

if seeds:
    selection = random.randint(0, len(seeds)-1)
    shutil.copy(seeds[selection], '%s/checkpointState.xml.bz2' % resultsdir)
