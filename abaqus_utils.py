'''
------------------------------------------------------------------------
Import Python-related libraries
------------------------------------------------------------------------
'''
import numpy as np
from math import *
import os
import sys
#todo: find out if I comment out the math part what breaks

'''
------------------------------------------------------------------------
Import Abaqus-related libraries
------------------------------------------------------------------------
'''
from abaqus import *
from odbAccess import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


#utility fxns
def printAB(string_):
    '''prints to the abaqus command line'''
    print >> sys.__stdout__, string_

def delete_extra_files(jname, additional_ext = None):
    '''
    deletes extra files:
    * jname: job name to delete files from
    * additional_ext: optional, pass a list of additional ext to remove jname.ext
    '''
    extra_file_extensions = ['.com', '.prt', '.sim', '.dmp', '.ipm', '.mdl', '.simlog', '.stt',
        '.pac', '.cid', '.sel', '.res', '.023', '.SMAFocus']
    if additional_ext is not None:
        extra_file_extensions.extend(additional_ext)
    for i in range(len(extra_file_extensions)):
        if os.path.exists(jname + extra_file_extensions[i]):
            os.remove(jname + extra_file_extensions[i])

def run_inp(jname, num_threads = 4):
    '''
    runs abaqus job
    * jname: runs jname.inp
    * num_threads: default 4
    '''
    if os.path.exists(jname + '.inp'):
        printAB("running " + jname)
        j = mdb.JobFromInputFile(activateLoadBalancing=False, atTime=None,
            explicitPrecision=DOUBLE_PLUS_PACK, inputFileName=jname,
            memory=90, memoryUnits=PERCENTAGE, multiprocessingMode=DEFAULT, name=
            jname, nodalOutputPrecision=SINGLE, numCpus=num_threads, numDomains=num_threads,
            parallelizationMethodExplicit=DOMAIN, queue=None, resultsFormat=ODB,
            scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

        j.submit(consistencyChecking = OFF)
        j.waitForCompletion()
        printAB("finished")

