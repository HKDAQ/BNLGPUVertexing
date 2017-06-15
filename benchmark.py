#!/bin/env python

from subprocess import call
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import OrderedDict
from itertools import product
from sys import stderr
from os.path import expandvars

parser = ArgumentParser(description='setup soft links for data', formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('--tag', '-t', default='gpu', help='Select a filename tag. For example, GPU type')
args = parser.parse_args()


###############
# SETUP DATASETS (i.e. the commands for setup_data.py)
# permutations
permutationDict = OrderedDict()
permutationDict['coverage'] = [x for x in ['14','40']]
permutationDict['time']     = [x for x in ['1ms','100mus','10mus']]
permutationDict['events']   = [x for x in ['1']]

# create a list of dictionaries for each permutation of the parameter values
permutationDictList = [ OrderedDict(zip(permutationDict, v)) for v in product(*permutationDict.values()) ]
datasets = []
for pDict in permutationDictList:
    command = ['python','setup_data.py']
    stub = ''
    #get the options
    for k,v in pDict.iteritems():
        #print k, v
        command.append('--' + k)
        command.append(v)
        stub += '_' + k[0] + v
        #print command
    datasets.append([command, stub])
#print datasets

###############
# SETUP BUILD OPTIONS (i.e. the commands for setup_data.py)
# permutations
permutationDict = OrderedDict()
permutationDict['histogram_type']   = [x for x in ['__HISTOGRAM_UCHAR__','__HISTOGRAM_USHORT__','__HISTOGRAM_UINT__']]
permutationDict['tof_table_type']   = [x for x in ['__TOF_TABLE_USHORT__','__TOF_TABLE_UINT__','__TOF_TABLE_FLOAT__']]
permutationDict['time_offset_type'] = [x for x in ['__TIME_OFFSET_USHORT__','__TIME_OFFSET_UINT__','__TIME_OFFSET_FLOAT__']]
permutationDict['v3_type']          = [x for x in ['__V3_UINT__','__V3_FLOAT__']]
permutationDict['pmt_id_type']      = [x for x in ['__PMT_ID_USHORT__','__PMT_ID_UINT__']]
permutationDict['l1cache_size']     = [x for x in ['__INCREASED_L1_CACHE__','']]

# create a list of dictionaries for each permutation of the parameter values
permutationDictList = [ OrderedDict(zip(permutationDict, v)) for v in product(*permutationDict.values()) ]
buildfiles = []
for pDict in permutationDictList:
    stub = ''
    buildfile = '#ifndef __BUILD_H__ \n' \
      '#define __BUILD_H__ \n'
    #get the options
    for k,v in pDict.iteritems():
        #print k, v
        if v:
            buildfile += '#define ' + v + '\n'
            stub += '_' + v.replace('_','')
    buildfile += '#endif //__BUILD_H__ \n'
    #print buildfile
    buildfiles.append([buildfile, stub])
#print buildfiles



def subprocess_call(command, stderr=None, stdout=None):
    try:
        retcode = call(command, stderr=stderr, stdout=stdout, shell=False)
        if retcode < 0:
            print >>stderr, "Child was terminated by signal", -retcode
            return True
        elif retcode > 0:
            print >>stderr, "Child returned", retcode
            return True
        else:
            return False
    except OSError as e:
        print >>stderr, "Execution failed:", e
        return True

#loop over datasets
for id, (dataset, datasetstub) in enumerate(datasets, 1):
    print 'Running over dataset', id, 'of', len(datasets)
    #print dataset
    if subprocess_call(dataset):
        continue
    #loop over build options
    for ib, (buildfile, buildfilestub) in enumerate(buildfiles, 1):
        print 'Running build', ib, 'of', len(buildfiles)
        with open('build.h', 'w') as f:
            f.write(buildfile)
            f.close()
        #compile the code
        if subprocess_call(['make','clean']):
            break #continue
        if subprocess_call(['make']):
            break #continue
        #run the code
        filestub = 'trig' + datasetstub + buildfilestub
        print 'STUB', filestub
        outfile = open(filestub+'.out', 'w')
        errfile = open(filestub+'.err', 'w')
        if subprocess_call(['nvprof','-o',filestub+'.nvpp','./daq_code'], stdout=outfile, stderr=errfile):
            break #continue
        outfile.close()
        errfile.close()
        #run the code again, with more metrics
        if subprocess_call(['nvprof','--kernels','kernel_correct_times_and_get_histo_per_vertex_shared','--metrics','all','-o',filestub+'_met.nvpp','./daq_code']):
            break #continue
        break
        
