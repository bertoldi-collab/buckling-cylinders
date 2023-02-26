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
from run_cylinder_lin_buckling_3d import *
from run_cylinder_nonlin_buckling_3d import *
from post_process_nonlin_buckling_3d import *

'''
------------------------------------------------------------------------
Import Abaqus-related libraries
------------------------------------------------------------------------
'''
from abaqus import *;
from odbAccess import *;
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
H       = 100.
RoH     = np.linspace(0.1,0.6,10)
toH     = np.linspace(0.02,0.15,10)
E1      = 1.
E2      = 1.
w       = 15.
i = [3,4,5,6,7]
j = [0,1,2,3,4]
ncores = int(24)
run_lin = True
run_nonlin = True

for ii in i:
	for jj in j:
		project = 'cyl_RoH_'+str(ii)+'_toH_'+str(jj)
		R = RoH[ii] * H
		t = toH[jj] * H
		run_cylinder_lin_buckling(project,H,R,t,w,E1,E2,run_lin)
		#Get the first positive eigenvalue
		f = open(project+'_lin_buckling_eigenvalues.txt','r')
		lines = f.readlines()
		eig_idx = 1
		for zz in range(len(lines)-1):
			eig = float(lines[zz+1].split()[-1])
			if eig>0:
				eig_idx = zz + 1
				break

		#Run nonlinear buckling
		imper = 0.25*t
		run_cylinder_nonlin_buckling(project,H,R,t,w,E1,E2,imper,eig_idx,ncores,run_nonlin)
