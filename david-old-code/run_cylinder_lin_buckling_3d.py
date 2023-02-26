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

def run_cylinder_lin_buckling(project,L,R,t,w,E1,E2,E3,ridge,ridge_folds,ridge_num,ridge_w,ridge_h,run_sim):
	''' 
	**Create Abaqus model, run FEA, and output results file
	'''
	'''
	------------------------------------------------------------------------
	Simulation options
	------------------------------------------------------------------------
		'''
	try:
		os.remove(project+'.lck')
		os.remove(project+'.odb')
	except OSError:
		pass

	#Given material properties
	nu      = 0.5 #Poisson's ratio to account for incompressibility
	mu1  	= E1/(2.*(1+nu)) #Initial shear modulus
	K1 		= 1e9 #Bulk modulus very high K = E/3*(1-2nu)
	mu2  	= E2/(2.*(1+nu)) #Initial shear modulus
	K2 		= 1e9 #Bulk modulus very high K = E/3*(1-2nu)
	mu3  	= E3/(2.*(1+nu)) #Initial shear modulus
	K3 		= 1e9 #Bulk modulus very high K = E/3*(1-2nu)

	#Abaqus related names
	model_name    = 'model'
	part_name     = 'part'
	assembly_name = 'instance'
	sketch_name   = 'sketch'
	cae_file      = project+'.cae'
	section_name  = 'solid_section'

	#Abaqus strain energy for neo-hokean form is
	# U = C_10 * (I_1 - 3) + 1/D_1(J^el-1)^2 (right part is zero for incomp.)
	nh_C10_1  	    = 0.5*mu1
	nh_D1_1		    = 0.
	nh_C10_2  	    = 0.5*mu2
	nh_D1_2		    = 0.
	nh_C10_3  	    = 0.5*mu3
	nh_D1_3		    = 0.
	material_name_1 = 'NH-1'
	material_name_2 = 'NH-2'
	material_name_3 = 'NH-3'

	#Element selection 
	quad_element = 'S4R'
	h_element    =  0.5*t#scales with thickness

	'''
	------------------------------------------------------------------------
	Abaqus Pre-processing
	------------------------------------------------------------------------
	'''

	'''MODEL'''
	Mdb()
	m = mdb.models['Model-1']

	'''SKETCHES'''
	# try:
	# 	del Model.sketches[sketch_name]
	# except:
	# 	pass
	# return
	s1 = m.ConstrainedSketch(name='s1', sheetSize=200.0)
	s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, R-0.5*t))
	s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, R+0.5*t))
	s2 = m.ConstrainedSketch(name='s2', sheetSize=200.0)
	s2.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, R+0.5*t))
	if ridge:
		s3 = m.ConstrainedSketch(name='s3', sheetSize=200.0)
		ridge_param = np.sqrt((R+0.5*t)*(R+0.5*t)-(0.5*ridge_w)*(0.5*ridge_w))
		s3.Line(point1=(ridge_param,0.5*ridge_w),point2=(ridge_param+ridge_h,0.5*ridge_w))
		s3.Line(point1=(ridge_param+ridge_h,0.5*ridge_w),point2=(ridge_param+ridge_h,-0.5*ridge_w))
		s3.Line(point1=(ridge_param+ridge_h,-0.5*ridge_w),point2=(ridge_param,-0.5*ridge_w))
		s3.Line(point1=(ridge_param,-0.5*ridge_w),point2=(ridge_param,0.5*ridge_w))

	'''PARTS'''
	p1 = m.Part(name=part_name+'_1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
	p1.BaseSolidExtrude(depth=L, sketch=s1)
	p2 = m.Part(name=part_name+'_2', dimensionality=THREE_D, type=DEFORMABLE_BODY)
	p2.BaseSolidExtrude(depth=w, sketch=s2)
	if ridge:
		p_ridge = m.Part(name=part_name+'_ridge', dimensionality=THREE_D, type=DEFORMABLE_BODY)
		p_ridge.BaseSolidExtrude(depth=L, sketch=s3)

	'''ASSEMBLY'''
	a         = m.rootAssembly
	i1 = a.Instance(dependent=ON, name=assembly_name+'_1',part=p1)
	i2 = a.Instance(dependent=ON, name=assembly_name+'_2',part=p2)
	i3 = a.Instance(dependent=ON, name=assembly_name+'_3',part=p2)
	if ridge:
		i4 = a.Instance(dependent=ON, name=assembly_name+'_ridge',part=p_ridge)
	a.translate(instanceList=(assembly_name+'_2', ), 
	vector=(0.0, 0.0, -w))
	a.translate(instanceList=(assembly_name+'_3', ), 
	vector=(0.0, 0.0, L))
	if ridge:
		a.RadialInstancePattern(axis=(0.0, 0.0, 1.0), 
		instanceList=(assembly_name+'_ridge', ), number=ridge_folds, point=(0.0, 0.0,0.0), totalAngle=360.)
	a.InstanceFromBooleanMerge(domain=GEOMETRY,instances=tuple([a.instances[a.instances.keys()[xxx]] for xxx in range(len(a.instances.keys()))]),
	keepIntersections=ON, name='Merged', originalInstances=SUPPRESS)
	p3 = m.parts['Merged']
	f1 = p3.faces.findAt((0.,0.,-w))
	f2 = p3.faces.findAt((0.,0.,L+w))
	p3.Set(name='bottom_cap_face',faces=(p3.faces[f1.index:f1.index+1], ))
	p3.Set(name='top_cap_face',faces=(p3.faces[f2.index:f2.index+1], ))
	c1 = p3.cells.findAt((0.,0.,-w))
	c2 = p3.cells.findAt((0.,0.,L+w))
	c3 = p3.cells.findAt((0.,R,0.5*L))
	cell_norms = [np.linalg.norm([np.array(p3.cells[xx].pointOn[0])[0],np.array(p3.cells[xx].pointOn[0])[1]]) for xx in range(len(p3.cells))]
	idx_ridge= np.argsort(cell_norms)[-ridge_num:]
	cridge = [p3.cells[xx] for xx in idx_ridge]
	p3.Set(name='bottom_cap_face',faces=(p3.faces[f1.index:f1.index+1], ))
	p3.Set(name='top_cap_face',faces=(p3.faces[f2.index:f2.index+1], ))
	p3.Set(name='bottom_cap_cell',cells=(p3.cells[c1.index:c1.index+1], ))
	p3.Set(name='top_cap_cell',cells=(p3.cells[c2.index:c2.index+1], ))
	p3.Set(name='body',cells=(p3.cells[c3.index:c3.index+1], ))
	if ridge:
		ridge_set = []
		for i in range(len(cridge)):
			ridge_set.append(p3.cells[cridge[i].index:cridge[i].index+1])

		p3.Set(name='ridge',cells=ridge_set)
		#Find connecting cells between body and ridges
		volume_cells = [float(p3.cells[xxx].getSize()) for xxx in range(len(p3.cells))]
		idx_sorted   = np.argsort(volume_cells)
		cell_connected = [p3.cells[idx_sorted[xx]] for xx in range(ridge_num)]
		connection_set = []
		for i in range(len(cridge)):
			connection_set.append(p3.cells[cell_connected[i].index:cell_connected[i].index+1])
		
		p3.Set(name='ridge_connection',cells=connection_set)


	'''SETS'''
	rp1=a.ReferencePoint(point=(0.0, 0.0, -w))
	rp2=a.ReferencePoint(point=(0.0, 0.0, 0.5*L))
	rp3=a.ReferencePoint(point=(0.0, 0.0, L+w))
	a.Set(name='rp1', referencePoints=(a.referencePoints[rp1.id], ))
	a.Set(name='rp2', referencePoints=(a.referencePoints[rp2.id], ))
	a.Set(name='rp3', referencePoints=(a.referencePoints[rp3.id], ))



	'''MATERIALS'''
	#Material
	rho1 = 1e-09
	mat1 = m.Material(name=material_name_1);
	mat1.Density(table=((rho1,),))
	mat1.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_1,nh_D1_1), ), testData=OFF,\
						 type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)
	rho2 = 1e-09
	mat2 = m.Material(name=material_name_2);
	mat2.Density(table=((rho2,),))
	mat2.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_2,nh_D1_2), ), testData=OFF,\
						 type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

	rho3 = 1e-09
	mat2 = m.Material(name=material_name_3);
	mat2.Density(table=((rho3,),))
	mat2.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_3,nh_D1_3), ), testData=OFF,\
						 type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

	'''SECTIONS'''

	m.HomogeneousSolidSection(material='NH-1', name='Section-1', thickness=None)
	m.HomogeneousSolidSection(material='NH-2', name='Section-2', thickness=None)
	m.HomogeneousSolidSection(material='NH-3', name='Section-3', thickness=None)

	p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
	p3.sets['body'], sectionName='Section-1', thicknessAssignment=FROM_SECTION)
	p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
	p3.sets['top_cap_cell'], sectionName='Section-2', thicknessAssignment=FROM_SECTION)
	p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
	p3.sets['bottom_cap_cell'], sectionName='Section-2', thicknessAssignment=FROM_SECTION)
	if ridge:
		p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
		p3.sets['ridge'], sectionName='Section-3', thicknessAssignment=FROM_SECTION)
		p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
		p3.sets['ridge_connection'], sectionName='Section-3', thicknessAssignment=FROM_SECTION)
	
	'''MESH'''
	p3.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=h_element)
	p3.setMeshControls(elemShape=TET,regions=p3.cells, technique=FREE)
	p3.setElementType(elemTypes=(ElemType(
	    elemCode=C3D8RH, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, 
	    hourglassControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
	    ElemType(elemCode=C3D4H, elemLibrary=STANDARD)), regions=(p3.cells, ))
	p3.generateMesh()

	'''CENTRAL NODES'''
	top_nodes = p3.sets['top_cap_face'].nodes
	top_norm = []
	ref_top = [rp3.xValue,rp3.yValue,rp3.zValue]
	for i in range(len(top_nodes)):
			top_norm.append(np.linalg.norm(ref_top-np.array(p3.sets['top_cap_face'].nodes[i].coordinates)))

	top_center_node = np.argmin(top_norm)
	p3.Set(name='top_middle_node',nodes=(p3.sets['top_cap_face'].nodes[top_center_node:top_center_node+1], ))

	bottom_nodes = p3.sets['bottom_cap_face'].nodes
	bottom_norm = []
	ref_bottom = [rp1.xValue,rp1.yValue,rp1.zValue]
	for i in range(len(bottom_nodes)):
			bottom_norm.append(np.linalg.norm(ref_bottom-np.array(p3.sets['bottom_cap_face'].nodes[i].coordinates)))

	bottom_center_node = np.argmin(bottom_norm)
	p3.Set(name='bottom_middle_node',nodes=(p3.sets['bottom_cap_face'].nodes[bottom_center_node:bottom_center_node+1], ))

	'''FLUID-CAVITY INTERACTION'''
	f3 = p3.faces.findAt((0.,0.,0.))
	f4 = p3.faces.findAt((0.,0.,L))
	f5 = p3.faces.findAt((R-0.5*t,0.,0.5*L))
	surfs = []
	surfs.append(a.Surface(name='Surf-1', side1Faces=a.instances['Merged-1'].faces[f3.index:f3.index+1]))
	surfs.append(a.Surface(name='Surf-2', side1Faces=a.instances['Merged-1'].faces[f4.index:f4.index+1]))
	surfs.append(a.Surface(name='Surf-3', side1Faces=a.instances['Merged-1'].faces[f5.index:f5.index+1]))


	a.SurfaceByBoolean(name='Surf-main', surfaces=(surfs))
	# a.Set(name='centernodes', nodes=a.instances['Merged-1'].nodes.getByBoundingBox(-2*R,-2*R,0.5*L-0.5*h_element,2*R,2*R,0.5*L+0.5*h_element))
	p3.Set(name='centernodes', nodes=p3.nodes.getByBoundingBox(-2*R,-2*R,0.5*L-0.9*h_element,2*R,2*R,0.5*L+0.1*h_element))

	m.FluidCavityProperty(bulkModulusTable=((2000.0, ), ), 
	    expansionTable=((1.0, ), ), fluidDensity=1e-09, name='IntProp-1', 
	    useBulkModulus=True, useExpansion=True)

	m.FluidCavity(cavityPoint=a.sets['rp2'], cavitySurface=a.surfaces['Surf-main'],
		createStepName='Initial', interactionProperty='IntProp-1', name='Int-1')

	'''AMPLITUDE'''
	m.SmoothStepAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name='Amp-1', timeSpan=STEP)

	'''STEP'''
	m.BuckleStep(maxEigen=None, name='Step-1', numEigen=5, previous='Initial', vectors=18,maxIterations=500)

	'''BCS'''
	m.DisplacementBC(amplitude=UNSET, buckleCase=
	    PERTURBATION_AND_BUCKLING, createStepName='Step-1', distributionType=
	    UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-1', region=
	    a.sets['Merged-1.bottom_cap_face'], u1=0.0, u2=0.0, u3=0.0, 
	    ur1=0.0, ur2=0.0, ur3=0.0)

	m.FluidCavityPressureBC(amplitude=UNSET, createStepName=
	    'Step-1', fixed=OFF, fluidCavity='Int-1', magnitude=-1.0, name='BC-2')
	'''WRITE INP FILE'''
	# dir_path = os.path.dirname(os.path.realpath(__file__))
	m.keywordBlock.synchVersions(storeNodesAndElements=False)
	nlinesInput = len(m.keywordBlock.sieBlocks)
	m.keywordBlock.replace(nlinesInput-2, 
	    '\n*Output, field, variable=PRESELECT\n*NODE FILE\nU,')
	j=mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
	    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
	    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
	    multiprocessingMode=DEFAULT, name=project+'_lin_buckling', nodalOutputPrecision=SINGLE, 
	    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
	    ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
	mdb.saveAs(cae_file)
	if run_sim:
		j.submit()
		j.waitForCompletion()
		odb = openOdb(path=project+'_lin_buckling.odb')
		step=odb.steps['Step-1']
		frames=step.frames
		initialframe = frames[1]
		f1=open(project+'_lin_buckling_eigenvalues.txt','w')
		f2=open(project+'_lin_buckling_displacement.txt','w')
		#Get eigenvalues
		for frame in frames:
			eig = frame.description
			f1.write(eig)
			f1.write('\n')
		#Get displacement of center nodes during first mode
		displacement = initialframe.fieldOutputs['U']
		center = odb.rootAssembly.instances['MERGED-1'].nodeSets['CENTERNODES']
		centerDisplacement = displacement.getSubset(region=center)
		for v in centerDisplacement.values:
			f2.write(str(v.magnitude))
			f2.write('\n')
		odb.close()
		f1.close()
		f2.close()
	return