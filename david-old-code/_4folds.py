'''
============================================================================
Python function to simulate hyper-elastic cylinder under deflation
NB: The units are N and mm (MPa).
Code: run_cylinder.py
Authors: David Melancon, Antonio Elia Forte, Ahmad Zareeeei
Harvard University
email: davidmelancon@g.harvard.edu

Tested on Abaqus 2018
============================================================================
'''
'''
------------------------------------------------------------------------
Import Python-related libraries
------------------------------------------------------------------------
'''
import numpy as np;
from math import *
from linBuck import *
from postBuck import *
from postProcess import *

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
import time
import os.path

def printAB(string_):
    print >> sys.__stdout__, string_
#end

#Variable parameters
#[[thickness, imper %, E_0, damping]]
#from yi: R=10mm, H=18mm, t=0.25 ~ 0.3 mm
params = [[0.54,0.004,1.0,0.0001]]
idx_try = 104

#Fix parameters- geometric param
H     = 18.
R     = 10. 
w     = 5.
theta = 180 * np.pi/180.

#Simulation parameters
ncores   = int(12)
runLin   = True
runPost  = True
hElement = 0.25 * 1.5/50.*H

for i in range(len(params)):
    t = params[i][0]
    imper = params[i][1]*t
    E1    = params[i][2]
    bdamp = params[i][3]
    t2    = t
    t1    = t   
    E2 = E1
    E3 = E1
    project = '4fold-imperfection-'+str(idx_try)
    printAB("t/H = " + str(t/H))
    printAB("R/H = " + str(R/H))
    #Run linear buckling
    printAB("version " + str(idx_try))
    printAB("running linear buckling")
    linBuck(project,H,R,t1,t2,theta,w,E1,E2,E3,hElement,runLin)

    #Get the first positive eigenvalue
    f = open(project+'_lin_buckling_eigenvalues.txt','r')
    lines = f.readlines()
    eig_idx = 1
    for i in range(len(lines)-1):
        eig = float(lines[i+1].split()[-1])
        if eig>0:
            eig_idx = i + 1
            break
    f.close()
    #Run nonlinear buckling

    printAB("running nonlinear buckling")
    postBuck(project,H,R,t1,t2,theta,w,E1,E2,E3,hElement,bdamp,imper,eig_idx,ncores,run_sim = runPost)

    #Run post-process
    postProcess(project+'_post_buckling')

#v100: just running david's code
#v101: running on 1st eigenvalue to lower volume
#v102: new geometry (h = 25, r = 10, t = 0.4 iirc)
#v103: new geometry from yi- too thin had 5 folds
#v104: using david geometry {h = 18, r = 10, t = 0.52-0.56 using mid value of 0.54}