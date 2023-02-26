'''
============================================================================
Python function to simulate hyper-elastic cylinder under deflation
NB: The units are N and mm (MPa).
Code: run_cylinder_lin_buclking_shell.py
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
import os
#import shutil

def linBuck(project,H,R,t1,t2,theta,w,E1,E2,E_cap,h_element,run_sim):
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
    mu1     = E1/(2.*(1+nu)) #Initial shear modulus
    K1      = 1e9 #Bulk modulus very high K = E/3*(1-2nu)
    mu2     = E2/(2.*(1+nu)) #Initial shear modulus
    K2      = 1e9 #Bulk modulus very high K = E/3*(1-2nu)
    mu_cap  = E_cap/(2.*(1+nu)) #Initial shear modulus
    K_cap   = 1e9 #Bulk modulus very high K = E/3*(1-2nu)

    #Abaqus related names
    model_name    = 'model'
    part_name     = 'part'
    assembly_name = 'instance'
    sketch_name   = 'sketch'
    cae_file      = project+'.cae'
    section_name  = 'solid_section'

    #Abaqus strain energy for neo-hokean form is
    # U = C_10 * (I_1 - 3) + 1/D_1(J^el-1)^2 (right part is zero for incomp.)
    nh_C10_1        = 0.5*mu1
    nh_D1_1         = 0.
    nh_C10_2        = 0.5*mu2
    nh_D1_2         = 0.
    nh_C10_cap      = 0.5*mu_cap
    nh_D1_cap       = 0.
    material_name_1 = 'NH-1'
    material_name_2 = 'NH-2'
    material_name_cap = 'NH-cap'

    #Element selection 
    quad_element = 'S4R'
    #h_element    =  3./50.*H#scales with thickness

    '''
    ------------------------------------------------------------------------
    Abaqus Pre-processing
    ------------------------------------------------------------------------
    '''

    '''MODEL'''
    #futuredir = os.getcwd()
    #os.chdir('/temp')
    #olddir   = os.getcwd()
    Mdb()
    m = mdb.models['Model-1']

    '''SKETCHES'''
    s1 = m.ConstrainedSketch(name='s1', sheetSize=200.0)
    s1.ArcByCenterEnds(center=(0.0, 0.0), direction=COUNTERCLOCKWISE,
        point1=(R*cos(0.5*pi-0.5*theta), R*sin(0.5*pi-0.5*theta)),
        point2=(0., R))
    s2 = m.ConstrainedSketch(name='s2', sheetSize=200.0)
    s2.ArcByCenterEnds(center=(0.0, 0.0), direction=COUNTERCLOCKWISE,
        point1=(0., R),point2=(-R*cos(0.5*pi-0.5*theta), R*sin(0.5*pi-0.5*theta)))
    s3 = m.ConstrainedSketch(name='s3', sheetSize=200.0)
    s3.ArcByCenterEnds(center=(0.0, 0.0), direction=COUNTERCLOCKWISE,
        point1=(-R*cos(0.5*pi-0.5*theta), R*sin(0.5*pi-0.5*theta)),
        point2=(R*cos(0.5*pi-0.5*theta), R*sin(0.5*pi-0.5*theta)))
    s4 = m.ConstrainedSketch(name='s4', sheetSize=200.0)
    s4.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.,R))

    '''PARTS'''
    p1 = m.Part(name=part_name+'_1', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p1.BaseShellExtrude(depth=H, sketch=s1)
    p2 = m.Part(name=part_name+'_2', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p2.BaseShellExtrude(depth=H, sketch=s2)
    p3 = m.Part(name=part_name+'_3', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p3.BaseShellExtrude(depth=H, sketch=s3)
    p4 = m.Part(name=part_name+'_4', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p4.BaseShell(sketch=s4)


    '''ASSEMBLY'''
    a = m.rootAssembly
    i1 = a.Instance(dependent=ON, name=assembly_name+'_1',part=p1)
    i2 = a.Instance(dependent=ON, name=assembly_name+'_2',part=p2)
    i3 = a.Instance(dependent=ON, name=assembly_name+'_3',part=p3)
    i4 = a.Instance(dependent=ON, name=assembly_name+'_4',part=p4)
    a.translate(instanceList=(assembly_name+'_4', ), 
    vector=(0.0, 0.0, H))
    a.InstanceFromBooleanMerge(domain=GEOMETRY,instances=tuple([a.instances[a.instances.keys()[xxx]] for xxx in range(len(a.instances.keys()))]),
        keepIntersections=ON, name='Merged', originalInstances=SUPPRESS)
    p4 = m.parts['Merged']

    f1 = p4.faces.findAt((0.,0.,H))
    f2 = p4.faces.findAt((R*cos(0.5*pi-0.25*theta),R*sin(0.5*pi-0.25*theta),0.5*H))
    f3 = p4.faces.findAt((-R*cos(0.5*pi-0.25*theta),R*sin(0.5*pi-0.25*theta),0.5*H))
    f4 = p4.faces.findAt((0.,-R,0.5*H))
    e1 = p4.edges.findAt((R*cos(0.5*pi-0.25*theta),R*sin(0.5*pi-0.25*theta),0.))
    e2 = p4.edges.findAt((-R*cos(0.5*pi-0.25*theta),R*sin(0.5*pi-0.25*theta),0.))
    e3 = p4.edges.findAt((0.,-R,0.))
    e4 = p4.edges.findAt((0.,R,0.5*H))
    p4.Set(name='cap_face',faces=(p4.faces[f1.index:f1.index+1], ))
    p4.Set(name='body-1',faces=(p4.faces[f2.index:f2.index+1],p4.faces[f3.index:f3.index+1], ))
    p4.Set(name='body-2',faces=(p4.faces[f4.index:f4.index+1], ))
    p4.Set(name='ring',edges=(p4.edges[e1.index:e1.index+1],p4.edges[e2.index:e2.index+1],p4.edges[e3.index:e3.index+1], ))
    p4.Set(name='curvature',edges=(p4.edges[e4.index:e4.index+1],))


    '''SETS'''
    rp1=a.ReferencePoint(point=(0.0, 0.0, 0.0))
    rp2=a.ReferencePoint(point=(0.0, 0.0, 0.5*H))
    rp3=a.ReferencePoint(point=(0.0, 0.0, H))
    a.Set(name='rp1', referencePoints=(a.referencePoints[rp1.id], ))
    a.Set(name='rp2', referencePoints=(a.referencePoints[rp2.id], ))
    a.Set(name='rp3', referencePoints=(a.referencePoints[rp3.id], ))

    '''MATERIALS'''
    #Material
    rho1 = 1e-09
    mat1 = m.Material(name=material_name_1);
    mat1.Density(table=((rho1,),))
    mat1.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_1,nh_D1_1), ), testData=OFF,type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)
    rho2 = 1e-09
    mat2 = m.Material(name=material_name_2);
    mat2.Density(table=((rho2,),))
    mat2.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_2,nh_D1_2), ), testData=OFF,type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

    rho_cap = 1e-09
    mat_cap = m.Material(name=material_name_cap);
    mat_cap.Density(table=((rho_cap,),))
    mat_cap.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_cap,nh_D1_cap), ), testData=OFF,type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

    '''SECTIONS'''
    m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='NH-1', name='Section-1', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=t1, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='NH-2', name='Section-2', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=t2, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material='NH-cap', name='Section-cap', 
        nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
        preIntegrate=OFF, temperature=GRADIENT, thickness=w, thicknessField='', 
        thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

    p4.SectionAssignment(offset=0.0, 
        offsetField='', offsetType=BOTTOM_SURFACE, region=
        p4.sets['body-1'], sectionName=
        'Section-1', thicknessAssignment=FROM_SECTION)
    p4.SectionAssignment(offset=0.0, 
        offsetField='', offsetType=BOTTOM_SURFACE, region=
        p4.sets['body-2'], sectionName=
        'Section-2', thicknessAssignment=FROM_SECTION)
    p4.SectionAssignment(offset=0.0, 
        offsetField='', offsetType=BOTTOM_SURFACE, region=
        p4.sets['cap_face'], sectionName=
        'Section-cap', thicknessAssignment=FROM_SECTION)

    '''MESH'''
    p4.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=h_element)
    #p3.setMeshControls(elemShape=H,regions=p3.cells, technique=FREE)
    p4.setElementType(elemTypes=(ElemType(
        elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
        hourglassControl=DEFAULT), ElemType(elemCode=S3, elemLibrary=STANDARD)), 
         regions=(p3.faces, ))
    p4.generateMesh()

    '''CENTRAL NODES'''
    top_nodes = p4.sets['cap_face'].nodes
    top_norm = []
    ref_top = [rp3.xValue,rp3.yValue,rp3.zValue]
    for i in range(len(top_nodes)):
            top_norm.append(np.linalg.norm(ref_top-np.array(p4.sets['cap_face'].nodes[i].coordinates)))

    top_center_node = np.argmin(top_norm)
    p4.Set(name='top_middle_node',nodes=(p4.sets['cap_face'].nodes[top_center_node:top_center_node+1], ))

    bottom_nodes = p4.sets['ring'].nodes
    bottom_norm = []
    ref_bottom = [rp1.xValue,rp1.yValue,rp1.zValue]
    for i in range(len(bottom_nodes)):
            bottom_norm.append(np.linalg.norm(ref_bottom-np.array(p4.sets['ring'].nodes[i].coordinates)))

    bottom_center_node = np.argmin(bottom_norm)
    p4.Set(name='bottom_middle_node',nodes=(p4.sets['ring'].nodes[bottom_center_node:bottom_center_node+1], ))

    '''FLUID-CAVITY INTERACTION'''
    surfs = []
    surfs.append(a.Surface(name='Surf-1', side2Faces=a.instances['Merged-1'].faces[f1.index:f1.index+1]))
    surfs.append(a.Surface(name='Surf-2', side2Faces=a.instances['Merged-1'].faces[f2.index:f2.index+1]))
    surfs.append(a.Surface(name='Surf-3', side2Faces=a.instances['Merged-1'].faces[f3.index:f3.index+1]))
    surfs.append(a.Surface(name='Surf-4', side2Faces=a.instances['Merged-1'].faces[f4.index:f4.index+1]))


    a.SurfaceByBoolean(name='Surf-main', surfaces=(surfs))
    # a.Set(name='centernodes', nodes=a.instances['Merged-1'].nodes.getByBoundingBox(-2*R,-2*R,0.5*H-0.5*h_element,2*R,2*R,0.5*H+0.5*h_element))
    p4.Set(name='centernodes', nodes=p4.nodes.getByBoundingBox(-2*R,-2*R,0.5*H-0.9*h_element,2*R,2*R,0.5*H+0.1*h_element))

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
        a.sets['Merged-1.ring'], u1=0.0, u2=0.0, u3=0.0, 
        ur1=0.0, ur2=0.0, ur3=0.0)

    m.FluidCavityPressureBC(amplitude=UNSET, createStepName=
        'Step-1', fixed=OFF, fluidCavity='Int-1', magnitude=-1.0, name='BC-2')

    '''WRITE INP FILE'''
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    m.keywordBlock.synchVersions(storeNodesAndElements=False)
    nlinesInput = len(m.keywordBlock.sieBlocks)
    m.keywordBlock.replace(nlinesInput-2, 
        '\n*Output, field, variable=PRESELECT\n*NODE FILE\nU,')
    j = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name=project+'_lin_buckling', nodalOutputPrecision=SINGLE, 
        numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.saveAs(project+'_lin_buckling')
    if run_sim:
        j.submit()
        j.waitForCompletion()
        odb = openOdb(path=project+'_lin_buckling.odb')
        step=odb.steps['Step-1']
        frames=step.frames
        initialframe = frames[1]
        f1=open(project+'_lin_buckling_eigenvalues.txt','w')
        f2=open(project+'_lin_buckling_displacement.txt','w')
        f3=open(project+'_lin_buckling_parameters.txt','w')
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

        f3.write('H, R, t1, t2, theta, w, E1, E2, E_cap, h_ele\n')
        f3.write(str(H)+' '+str(R)+' '+str(t1)+' '+str(t2)+' '+str(theta)+' '+str(w)+' '+str(E1)+' '+str(E2)+' '+str(E_cap)+' '+str(h_element))
        f1.close()
        f2.close()
        f3.close()
    return