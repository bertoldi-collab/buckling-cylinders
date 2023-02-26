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
import numpy as np
from math import *
from run_cylinder_lin_buckling_3d import *
from run_cylinder_nonlin_buckling_3d import *
from post_process_nonlin_buckling_3d import *

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

#Parameters
H = 20. #3folds
# H       = 12. #4folds
# RoH     = np.linspace(0.1,0.6,10)
# toH     = np.linspace(0.02,0.15,10)
E1      = 1.
E2      = 1.
E3      = 1.
w       = 5. #3folds
# w       = 3. #4folds
i = 5
j = 0
R = 8. #3folds
# R = 6.5 #4folds
t = 0.5 #3 folds
# t = 0.35 #4folds
idx_try = 100
project = '3folds-test-3d-' + str(idx_try)
ncores = int(8)
run_lin = True
run_nonlin = True

#garbage
ridge   = False
ridge_folds = 0
ridge_num   = 0
ridge_w = 2*t
ridge_h = 5*t
#Run linear buckling
run_cylinder_lin_buckling(project,H,R,t,w,E1,E2,E3,ridge,ridge_folds,ridge_num,ridge_w,ridge_h,run_lin)

#Get the first positive eigenvalue
f = open(project+'_lin_buckling_eigenvalues.txt','r')
lines = f.readlines()
eig_idx = 1
for i in range(len(lines)-1):
	eig = float(lines[i+1].split()[-1])
	if eig>0:
		eig_idx = i + 1
		break

#Run nonlinear buckling
imper = 0.3*t
run_cylinder_nonlin_buckling(project,H,R,t,w,E1,E2,E3,imper,eig_idx,ridge,ridge_folds,ridge_num,ridge_w,ridge_h,ncores,run_nonlin)

#Run post-process
post_process(project+'_nonlin_buckling')
