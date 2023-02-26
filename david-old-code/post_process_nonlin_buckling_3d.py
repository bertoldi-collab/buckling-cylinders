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

def post_process(project):
	odb = openOdb(path=project+'.odb')
	#Get pressure and volume
	out_name_1 = odb.steps['Step-1'].historyRegions.keys(0)[1];
	out_his_1 = odb.steps['Step-1'].historyRegions[out_name_1].historyOutputs['CVOL'].data;
	out_his_2 = odb.steps['Step-1'].historyRegions[out_name_1].historyOutputs['PCAV'].data;
	if out_his_1[-1][0] == out_his_1[-2][0]:
		n_end = len(out_his_1)-1
	else:
		n_end = len(out_his_1)

	time      = [out_his_1[ii][0] for ii in range(n_end)];
	volume    = [out_his_1[ii][1] for ii in range(n_end)];
	pressure  = [out_his_2[ii][1] for ii in range(n_end)];

	#Get the nodes from top and bottom surfaces
	top_nodes = odb.rootAssembly.instances['MERGED-1'].nodeSets['TOP_CAP_FACE']
	bot_nodes = odb.rootAssembly.instances['MERGED-1'].nodeSets['BOTTOM_CAP_FACE']

	#Get height of the cylinder
	H = top_nodes.nodes[0].coordinates[2]-bot_nodes.nodes[0].coordinates[2]
	centroid_bottom = np.array([0.,0.,bot_nodes.nodes[0].coordinates[2]])
	# #Get displacements and rotations
	# out_name_2 = odb.steps['Step-1'].historyRegions.keys(0)[2];
	# out_his_3 = odb.steps['Step-1'].historyRegions[out_name_2].historyOutputs['U1'].data;
	# out_his_4 = odb.steps['Step-1'].historyRegions[out_name_2].historyOutputs['U2'].data;
	# out_his_5 = odb.steps['Step-1'].historyRegions[out_name_2].historyOutputs['U3'].data;
	# out_his_6 = odb.steps['Step-1'].historyRegions[out_name_2].historyOutputs['UR1'].data;
	# out_his_7 = odb.steps['Step-1'].historyRegions[out_name_2].historyOutputs['UR2'].data;
	# out_his_8 = odb.steps['Step-1'].historyRegions[out_name_2].historyOutputs['UR3'].data;

	# #Calculate extension
	# extension = [1./H*np.sqrt(out_his_3[ii][1]*out_his_3[ii][1] + out_his_4[ii][1]*out_his_4[ii][1]+
	# 			(out_his_5[ii][1]+H)*(out_his_5[ii][1]+H)) for ii in range(n_end)]

	# #Calculate bending angle (from centroids)
	# v1 = [np.sqrt(out_his_3[ii][1]*out_his_3[ii][1] + out_his_4[ii][1]*out_his_4[ii][1]+
	# 			(out_his_5[ii][1]+H)*(out_his_5[ii][1]+H)) for ii in range(n_end)]
	# v2 = [np.sqrt(out_his_3[ii][1]*out_his_3[ii][1] + out_his_4[ii][1]*out_his_4[ii][1]) for ii in range(n_end)]
	# angle = [asin(v2[ii]/v1[ii])*180/pi for ii in range(n_end)]

	#Get strain energy
	out_name_3 = odb.steps['Step-1'].historyRegions.keys(0)[0];
	out_his_9 = odb.steps['Step-1'].historyRegions[out_name_3].historyOutputs['ALLSE'].data;
	SE = [out_his_9[ii][1] for ii in range(n_end)];

	#Get bending angle and twist angle from plane normals (RIGID BODY MOTION)
	idx = 1
	X_ref = np.array([top_nodes.nodes[ii].coordinates[0] for ii in range(len(top_nodes.nodes))])
	Y_ref = np.array([top_nodes.nodes[ii].coordinates[1] for ii in range(len(top_nodes.nodes))])
	Z_ref = np.array([top_nodes.nodes[ii].coordinates[2] for ii in range(len(top_nodes.nodes))])
	A_ref = np.array([[top_nodes.nodes[ii].coordinates[0],top_nodes.nodes[ii].coordinates[1],1] for ii in range(len(top_nodes.nodes))])
	v1 = np.array([X_ref[idx],Y_ref[idx],115.])-np.array([X_ref[0],Y_ref[0],115.])
	v1 = v1/np.linalg.norm(v1)
	n_ref = np.array([0,0,1])
	v2 = np.cross(n_ref,v1)
	ba = np.zeros((n_end,1))
	ext = np.zeros((n_end,1))
	ext[0] = 1.
	twist = np.zeros((n_end,1))
	for i in range(1,n_end):
		frame = odb.steps['Step-1'].frames[i]
		disp  = frame.fieldOutputs['U']
		top_disp  = disp.getSubset(region=top_nodes)
		Z_temp = Z_ref + np.array([top_disp.values[ii].data[2] for ii in range(len(top_disp.values))])
		XX_temp = X_ref + np.array([top_disp.values[ii].data[0] for ii in range(len(top_disp.values))])
		YY_temp = Y_ref + np.array([top_disp.values[ii].data[1] for ii in range(len(top_disp.values))])
		centroid_temp = np.array([np.mean(XX_temp),np.mean(YY_temp),np.mean(Z_temp)])
		ext[i] = np.linalg.norm(centroid_temp-centroid_bottom)/H
		A_temp = A_ref + np.array([[top_disp.values[ii].data[0],top_disp.values[ii].data[1],1] for ii in range(len(top_disp.values))])
		x_temp = np.dot(np.dot(np.linalg.inv(np.dot(A_temp.T,A_temp)),A_temp.T),Z_temp) #this is solving the plane equation using the normal equations
		normal_temp = np.array([-x_temp[0],-x_temp[1],1])/np.linalg.norm([-x_temp[0],-x_temp[1],1])
		ba_temp = np.dot(n_ref,normal_temp.T) #angle between the two vectors
		if ba_temp >1.:
			ba[i] = 0.
		else:
			ba[i] = np.arccos(ba_temp)
		linenodes = np.cross(n_ref,normal_temp) #Axis of rotation between the two planes
		normlinenodes = np.linalg.norm(linenodes)
		if normlinenodes > 0.0:
			linenodes = linenodes/normlinenodes
		ux = linenodes[0]
		uy = linenodes[1]
		uz = linenodes[2]
		#Rotation matrix depending on axis of rotation + bending angle (ba)
		rotmat = np.array([[np.cos(ba[i])+ux**2*(1-np.cos(ba[i])),ux*uy*(1-np.cos(ba[i]))-uz*np.sin(ba[i]), ux*uz*(1-np.cos(ba[i]))+uy*np.sin(ba[i])],
						[uy*ux*(1-np.cos(ba[i]))+uz*np.sin(ba[i]),np.cos(ba[i])+uy**2*(1-np.cos(ba[i])),uy*uz*(1-np.cos(ba[i]))-ux*np.sin(ba[i])],
						[uz*ux*(1-np.cos(ba[i]))-uy*np.sin(ba[i]),uz*uy*(1-np.cos(ba[i]))+ux*np.sin(ba[i]),np.cos(ba[i])+uz**2*(1-np.cos(ba[i]))]])
		v3 = np.array([XX_temp[idx],YY_temp[idx],x_temp[0]*XX_temp[idx]+x_temp[1]*YY_temp[idx]+x_temp[2]])-np.array([XX_temp[0],YY_temp[0],x_temp[0]*XX_temp[0]+x_temp[1]*YY_temp[0]+x_temp[2]])
		v4 = np.cross(normal_temp,v3)
		v2_rot = np.dot(rotmat[:,:,0],v2)
		temp = np.dot(v4,v2_rot)/(np.linalg.norm(v4)*np.linalg.norm(v2_rot))
		if temp > 1.:
				twist[i] = 0.
		else:
				twist_temp = np.arccos(temp)
				twist[i] = twist_temp
	odb.close()

	#Write results
	f=open(project+'_results.txt','w')
	#Write the results in a text file
	f.write('time')
	f.write(' ')
	f.write('Volume')
	f.write(' ')
	f.write('Pressure')
	f.write(' ')
	f.write('Energy')
	f.write(' ')
	f.write('Extension')
	f.write(' ')
	# f.write('Bending-centroid')
	# f.write(' ')
	f.write('Bending-normals')
	f.write(' ')
	f.write('Twist')
	f.write("\n")

	for i in range (len(time)):
		f.write(str(time[i]))
		f.write(' ')
		f.write(str(volume[i]))
		f.write(' ')
		f.write(str(pressure[i]))
		f.write(' ')
		f.write(str(SE[i]))
		f.write(' ')
		f.write(str(ext[i][0]))
		f.write(' ')
		# f.write(str(angle[i]))
		# f.write(' ')
		f.write(str(ba[i][0]))
		f.write(' ')
		f.write(str(twist[i][0]))
		f.write("\n")

	f.close()

	return