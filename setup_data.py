#!/bin/env python

import argparse
import os
import os.path

parser = argparse.ArgumentParser(description='setup soft links for data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--coverage', '-c', choices=['14','40'], required=True, help='Select photocoverage')
parser.add_argument('--time', '-t', choices=['1ms','100mus','10mus'], required=True, help='Select dataset duration')
parser.add_argument('--events', '-e', default=1, type=int, help='Number of events')

args = parser.parse_args()

#detector
if os.path.islink('detector.txt'):
    os.unlink('detector.txt')
os.symlink('hits/detector_' + args.coverage + '.txt', 'detector.txt')

#pmts
if os.path.islink('all_pmts.txt'):
    os.unlink('all_pmts.txt')
os.symlink('hits/all_pmts_' + args.coverage + 'per.txt', 'all_pmts.txt')

#data
i = 0
while True:
    if os.path.islink('all_hits_' + str(i+1) + '.txt'):
        os.unlink('all_hits_' + str(i+1) + '.txt')
    else:
        break
    i += 1
for i in xrange(args.events):
    os.symlink('hits/all_hits_' + args.coverage + 'per_' + args.time + '.txt', 'all_hits_' + str(i+1) + '.txt')
