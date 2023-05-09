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

#porting to object oriented code

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
import sys

def printAB(string_):
    print >> sys.__stdout__, string_

def delete_extra_files(jname, additional_ext = None):
    #delete extra files
    extra_file_extensions = ['.com', '.prt', '.sim', '.dmp', '.ipm', '.mdl', '.simlog', '.stt',
        '.pac', '.cid', '.sel', '.res', '.023', '.SMAFocus']
    if additional_ext is not None:
        extra_file_extensions.extend(additional_ext)
    for i in range(len(extra_file_extensions)):
        if os.path.exists(jname + extra_file_extensions[i]):
            os.remove(jname + extra_file_extensions[i])

def run_inp(jname, num_threads = 4):
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

class full_shell:
    def __init__(self, project, imperfection = None, fullProps = None, simpProps = None):
        self.project = project
        self.post_instability = True

        if fullProps is not None:
            self.R = float(fullProps.R)
            self.H = float(fullProps.H)
            self.w = float(fullProps.w)
            self.theta = float(fullProps.theta)

            self.t1 = float(fullProps.t1)
            self.t2 = float(fullProps.t2)

            self.E1 = float(fullProps.E1)
            self.E2 = float(fullProps.E2)
            self.E_cap = float(fullProps.E_cap)
        elif simpProps is not None:
            self.R = float(simpProps.R)
            self.H = float(simpProps.H)
            self.w = float(simpProps.w)
            self.theta = np.pi

            self.t1 = float(simpProps.t)
            self.t2 = self.t1

            self.E1 = float(simpProps.E)
            self.E2 = self.E1
            self.E_cap = self.E1
        else:
            raise ValueError("pls input something thx!!!")

        if imperfection is not None:
            self.imperfection = imperfection * self.t1
        else:
            #set default value of 0.004
            self.imperfection = 0.004 * self.t1


        #finish any other initialization
        # original mult: 0.25 * 1.5/50
        self.h_element = 1.5/50 * self.H
        self.temp_set = -0.332 #removes whole volume
        self.nu_shell = 0.5

        self.mesh_shape = 'quad' #pass in 'tri' for triangular elements
        self.mesh_order = 'lin' #pass in 'quadratic' for 2nd order


    def make_job(self,extra_str,num_threads = 1):
        m = self.model
        m.rootAssembly.regenerate()
        jname = self.project + extra_str
        j = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=jname, nodalOutputPrecision=SINGLE, 
            numCpus=num_threads, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        j.writeInput()
        return jname

    def run_linear_model(self):
        jname = self.make_linear_model()
        run_inp(jname)
        self.post_process_lin_buckle(eig_name = '_lin_buckling')
        return jname

    def make_nonlin_multi_buckle(self,bdamp, max_temp_mult = 0.6, num_steps = 10, eig_idx = None):
        eig_name = '_lin_buckling'
        self.bdamp = bdamp

        max_temp = max_temp_mult * -0.332
        temp_all = np.linspace(0, max_temp, num_steps + 1)[1:]

        self.make_geometry(nonlinear_model = True)
        if eig_idx is None:
            eig_idx = self.get_eig_idx(eig_name = eig_name)
        self.finish_nonlinear_initial()

        self.finish_nonlinear_steps(temp_list = temp_all)
        self.add_imperfection(eig_idx, eig_name = eig_name)
        jname = self.save_cae_write_job('_multi_buckling')

        return jname
        
        #uhh so I want to basically have a normal nonlin model
        #but keep adding more and more steps changing temp_set each time

    def make_riks_model(self,bdamp, temp_mult = None):
        #makes dynamic implicit by default, if is_buckling is true makes a static followed by buckle step
        if temp_mult is not None:
            if temp_mult <= 1 and temp_mult >= 0:
                self.temp_set = -0.332 * temp_mult
            else:
                raise ValueError("temp mult is out of bounds!")

        eig_name = '_lin_buckling'
        extra_str = '_riks_buckling'

        self.bdamp = bdamp
        self.make_geometry(nonlinear_model = True)
        eig_idx = self.get_eig_idx(eig_name = eig_name)
        self.finish_nonlinear_initial()
        self.finish_nonlinear_riks()
        self.add_imperfection(eig_idx, eig_name = eig_name)

        

        jname = self.save_cae_write_job(extra_str)
        return jname

    def make_nonlin_model(self,bdamp, is_buckling = False, temp_set = None, eig_idx = None, eig_name = '_lin_buckling',  alt_name = None):
        '''
        *bdamp: damping value
        *is_buckling: default false- does a dyn imp step; if set to true then does a single [static, freq] step
        *temp_set: default None, allows setting final temp manually
        *eig_idx: default None, allows setting index to pull impefection mode from eig_name manually
        *eig_name: default '_lin_buckling', defines fil file to pull buckling mode from
        *alt_name: default '_post_buckling', defines name of job file
        '''

        if temp_set is not None:
            if temp_set >= -0.332 and temp_set <= 0:
                self.temp_set = temp_set
            else:
                raise ValueError("temp set is out of bounds!")

        if alt_name is not None:
            extra_str = alt_name
            self.post_process_lin_buckle(eig_name = eig_name)
        else:
            extra_str = '_post_buckling'
        
        if eig_idx is None:
            eig_idx = self.get_eig_idx(eig_name = eig_name)

        self.bdamp = bdamp
        self.make_geometry(nonlinear_model = True)
        self.finish_nonlinear_initial()
        self.finish_nonlinear_steps(make_dyn = not is_buckling)
        self.add_imperfection(eig_idx, eig_name = eig_name)

        jname = self.save_cae_write_job(extra_str)
        return jname

    def make_linear_model(self):
        self.make_geometry()
        jname = self.finish_linear_buckle()
        return jname

    def save_cae_write_job(self, extra_str):
        jname = self.make_job(extra_str)
        mdb.saveAs(self.project+extra_str)
        return jname

    def post_process_lin_buckle(self, eig_name):
        odb = openOdb(path=self.project + eig_name + '.odb')
        #todo: can I just grab the last step instead?
        if eig_name is not '_lin_buckling':
            step=odb.steps['Step-1-buckle']
        else:
            step=odb.steps['Step-1']
        frames=step.frames
        initialframe = frames[1]
        #f1 is used
        f1=open(self.project + eig_name + '_eigenvalues.txt','w')
        # f2=open(self.project + eig_name + '_displacement.txt','w')
        # f3=open(self.project + eig_name + '_parameters.txt','w')

        #f1: Get eigenvalues
        for frame in frames:
            eig = frame.description
            f1.write(eig)
            f1.write('\n')
        f1.close()
        
        #f2: Get displacement of center nodes during first mode
        # displacement = initialframe.fieldOutputs['U']
        # center = odb.rootAssembly.instances['MERGED-1'].nodeSets['CENTERNODES']
        # centerDisplacement = displacement.getSubset(region=center)
        # for v in centerDisplacement.values:
        #     f2.write(str(v.magnitude))
        #     f2.write('\n')
        odb.close()

        #f3: write parameters: changing this to write to data-out where it can be used
        # props_all = np.array([self.a, self.b, self.c, self.d, self.alpha_11, self.alpha_22, final_temp])
        props_all = np.array([self.H, self.R, self.t1, self.t2, self.theta, self.w, self.E1, self.E2, self.E_cap, self.h_element])
        np.savetxt("../data_out/" + self.project + '_props.txt', props_all)
        # f3.write('H, R, t1, t2, theta, w, E1, E2, E_cap, h_ele\n')
        # f3.write(str(self.H)+' '+str(self.R)+' '+str(self.t1)+' '+str(self.t2)+' '+
        #     str(self.theta)+' '+str(self.w)+' '+str(self.E1)+' '+str(self.E2)+' '+str(self.E_cap)+' '+str(self.h_element))
        # f2.close()
        # f3.close()

    def get_eig_idx(self, eig_name):
        f = open(self.project + eig_name + '_eigenvalues.txt','r')
        lines = f.readlines()
        #finds first positive eigenvalue and returns the index (1-indexed)
        eig_idx = 1
        for i in range(len(lines)-1):
            eig = float(lines[i+1].split()[-1])
            if eig>0:
                eig_idx = i + 1
                break
        f.close()
        return eig_idx

    def finish_linear_buckle(self):
        m = self.model
        a = self.aa

        '''STEP'''
        m.BuckleStep(maxEigen=None, name='Step-1', numEigen=5, previous='Initial', vectors=18,maxIterations=500)

        '''BCS'''
        m.DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, createStepName='Step-1',
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-1',
            region=a.sets['Merged-1.ring'],u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0)

        m.FluidCavityPressureBC(amplitude=UNSET, createStepName=
            'Step-1', fixed=OFF, fluidCavity='Int-1', magnitude=-1.0, name='BC-2')

        '''WRITE INP FILE'''
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        m.keywordBlock.synchVersions(storeNodesAndElements=False)
        nlinesInput = len(m.keywordBlock.sieBlocks)
        m.keywordBlock.replace(nlinesInput-2, '\n*Output, field, variable=PRESELECT\n*NODE FILE\nU,')

        extra_str = '_lin_buckling'
        jname = self.make_job(extra_str)
        mdb.saveAs(self.project+extra_str)
        return jname

    def add_imperfection(self, eig_idx, eig_name):
        m = self.model
        imper = self.imperfection

        if eig_name is not '_lin_buckling':
            step_num = 2
        else:
            step_num = 1

        '''WRITE IMPERFECTION'''
        m.keywordBlock.synchVersions(storeNodesAndElements=False)
        keyWords = m.keywordBlock.sieBlocks
        keyLine = '** ----------------------------------------------------------------\n** \n** STEP: Step-1\n** '
        idx_change = keyWords.index(keyLine)
        m.keywordBlock.replace(idx_change,
            '** ----------------------------------------------------------------\n*IMPERFECTION,FILE='+
            self.project+eig_name+',STEP='+str(step_num)+'\t\n'+str(eig_idx)+','+str(imper)+'\n** \n** STEP: Step-1\n**')

    def finish_nonlinear_initial(self):
        #missing: steps, requests, temp_set bcs
        m = self.model
        a = self.aa

        set_rp2 = a.sets['rp2']


        '''CONTACT'''
        contact_prop = m.ContactProperty('IntProp-2')
        contact_prop.NormalBehavior(allowSeparation=ON, constraintEnforcementMethod=DEFAULT, pressureOverclosure=HARD)

        int_contact = m.ContactStd(createStepName='Initial', name='Int-2')
        int_contact.includedPairs.setValuesInStep(stepName='Initial', useAllstar=ON)
        int_contact.contactPropertyAssignments.appendInStep(assignments=((GLOBAL, SELF, 'IntProp-2'), ), stepName='Initial')


        '''BCS'''
        m.DisplacementBC(amplitude=UNSET, createStepName='Initial',distributionType=UNIFORM, fieldName='',
                fixed=OFF,localCsys=None, name='BC-ring', region=a.sets['Merged-1.ring'],
                u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
        # m.VelocityBC(amplitude=UNSET, createStepName='Initial',distributionType=UNIFORM, 
        #     fieldName='', localCsys=None, name='BC-2', region=a.sets['Merged-1.ring'],
        #     v1=0.0, v2=0.0, v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0)

        '''TEMPERATURE FIELDS FOR VOLUME CHANGE'''
        #self.temp_set was originally -0.332 for taking out full volume
        m.Temperature(createStepName='Initial', crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
            distributionType=UNIFORM, magnitudes=(0.0, ), name='Predefined Field-1', region=a.sets['rp2'])

    def finish_nonlinear_riks(self):
        m = self.model
        # a = self.aa

        set_rp2 = self.aa.sets['rp2']

        nincre = 5000
        nincre_output = 200 #used to be 5000

        '''AMPLITUDE'''
        m.TabularAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name='Amp-1', timeSpan=STEP)

        m.StaticRiksStep(name='Step-1', nlgeom=ON, previous='Initial', initialArcInc=0.01, maxArcInc=1e+36,
            maxNumInc=nincre, minArcInc=1e-09)

        '''REQUESTS''' 
        #general- field output
        # m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=nincre_output,
        #     variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U','RF', 'CF', 'CSTRESS', 'CDISP', 'COORD'))

        #temp bd condition
        m.Temperature(amplitude='Amp-1', createStepName='Step-1', crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
        distributionType=UNIFORM, magnitudes=(self.temp_set, ), name='Predefined Field-2', region=set_rp2)

    def finish_nonlinear_steps(self, make_dyn = False, make_riks = False, temp_list = None):
        '''
        *default behavior is to make 1 static and 1 buckle step
        *if make_dyn is True, then it makes a dynamic step
        *if make_riks is True, then it makes a riks step
        *if temp_list is provided it makes a series of [static, buckle] steps
        '''

        m = self.model
        # a = self.aa

        set_rp2 = self.aa.sets['rp2']

        nincre = 5000
        nincre_output = 200 #used to be 5000

        if temp_list is None:
            temp_list = [self.temp_set]

        if make_dyn:
            '''AMPLITUDE'''
            m.SmoothStepAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name='Amp-1', timeSpan=STEP)

            m.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
               application=QUASI_STATIC, initialConditions=OFF, initialInc=0.001, maxInc=0.01,
               maxNumInc=nincre, minInc=1e-09, name='Step-1', nlgeom=ON, nohaf=OFF, previous='Initial')


            '''REQUESTS''' 
            #general- field output
            m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=nincre_output, timeMarks=OFF,
                variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U','RF', 'CF', 'CSTRESS', 'CDISP', 'COORD'))

            #temp bd condition
            m.Temperature(amplitude='Amp-1', createStepName='Step-1', crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                distributionType=UNIFORM, magnitudes=(self.temp_set, ), name='Predefined Field-2', region=set_rp2)
        
        else:
            prev_step_name = 'Initial'
            for i in range(len(temp_list)):
                current_temp = temp_list[i]
                idx = i + 1

                #create static step for moving, then buckle step for analysis
                max_inc_step = np.min([0.05*len(temp_list), 1])
                static_step_name = 'Step-' + str(idx)
                buckle_step_name = static_step_name + '-buckle'
                static_step = m.StaticStep(initialInc=0.01*len(temp_list), maxInc=max_inc_step, maxNumInc=nincre, minInc=1e-11,
                    name=static_step_name, nlgeom=ON, previous=prev_step_name)
                # static_step.control.setValues(allowPropagation=OFF,resetDefaultValues=OFF,
                #     timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0,12.0, 25.0, 6.0, 3.0, 50.0))
                #note: the 25.0 above is the number of 1U, 2U,...etc attempts, default is 5
                # m.BuckleStep(maxEigen=None, name=buckle_step_name, numEigen=10,
                #     previous=static_step_name, vectors=18,maxIterations=1000)
                m.FrequencyStep(name=buckle_step_name, numEigen=10, previous=static_step_name)

                #set temp for static
                current_amp_name = 'Amp-temp-' + str(idx)
                m.TabularAmplitude(data=((0.0, float(i)/idx), (1.0, 1.0)), name=current_amp_name, timeSpan=STEP)
                # if i > 0:
                #     m.SmoothStepAmplitude(data=((0.0, float(i)/idx), (1.0, 1.0)), name='Amp-' + str(idx), timeSpan=STEP)

                m.Temperature(amplitude=current_amp_name, createStepName=static_step_name,
                    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=UNIFORM,
                    magnitudes=(current_temp, ), name='Temp-Move-' + str(idx), region=set_rp2)
                # m.FluidCavityPressureBC(amplitude=UNSET, createStepName=buckle_step_name,
                #     fixed=OFF, fluidCavity='Int-1', magnitude=-1.0,  name='BC-buckle-' + str(idx))


                prev_step_name = buckle_step_name

            #add keyword for making fil file- how does this interact w/ multiple buckle steps?
            #currently we're just saving the fil for the last buckle- maybe it's fine bc the imperfection asks for step name?
            m.keywordBlock.synchVersions(storeNodesAndElements=False)
            nlinesInput = len(m.keywordBlock.sieBlocks)
            m.keywordBlock.replace(nlinesInput-2, '\n*Output, field, variable=PRESELECT\n*NODE FILE\nU,')

        #history outputs
        m.historyOutputRequests['H-Output-1'].setValues(numIntervals=nincre_output, timeMarks=OFF, variables=('ALLSE', ))
        m.HistoryOutputRequest(createStepName='Step-1', name='H-Output-1-cavity', numIntervals=nincre_output, rebar=EXCLUDE, 
            region=set_rp2, sectionPoints=DEFAULT, timeMarks=OFF, variables=('PCAV', 'CVOL', 'NT'))

        #field outputs
        m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=nincre_output, timeMarks=OFF,
            variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U','RF', 'CF', 'CSTRESS', 'CDISP', 'COORD'))

    def post_process_multi_buckle(self):
        #Import the relavent data
        project = self.project + '_multi_buckling'
        part_name = 'Merged'
        # Post-processing
        odb = openOdb(project + '.odb')

        # printAB(odb.steps)
        # printAB(len(odb.steps))
        num_freq_steps = len(odb.steps)/2

        data_all = np.empty((num_freq_steps, 11))
        data_all[:] = np.nan

        data_eig_val = np.empty((num_freq_steps, 11))
        data_eig_val[:] = np.nan

        for i in range(num_freq_steps):
            idx = i + 1
            prev_step_name = 'Step-' + str(idx)
            step_name = 'Step-' + str(idx) + '-buckle'
            step_prev = odb.steps[prev_step_name]
            step = odb.steps[step_name]
            # printAB(step)
            # printAB('***')
            # printAB(step_prev.historyRegions)
            his_region_prev = step_prev.historyRegions['Node ASSEMBLY.2']
            # printAB(his_region_prev.historyOutputs)
            # printAB('***')
            his_name = step.historyRegions.keys()[0]
            his_region = step.historyRegions[his_name]
            # printAB(his_region.historyOutputs)
            # printAB(his_region.historyOutputs['EIGFREQ'].data)
            out_his_eig_freq = his_region.historyOutputs['EIGFREQ'].data
            out_his_eig_val = his_region.historyOutputs['EIGVAL'].data
            out_his_nt = his_region_prev.historyOutputs['NT11'].data
            # printAB(out_his_nt)
            if out_his_eig_freq is not None:
                delta_T = out_his_nt[1][-1]
                eig_freq = np.asarray([out_his_eig_freq[j][1] for j in range(len(out_his_eig_freq))])
                eig_val = np.asarray([out_his_eig_val[j][1] for j in range(len(out_his_eig_val))])
                data_all[i,0] = delta_T
                data_all[i,1:] = eig_freq

                data_eig_val[i,0] = delta_T
                data_eig_val[i,1:] = eig_val
            # printAB('---')
        # printAB(data_all)
        odb.close()
        np.savetxt("../data_out/" + self.project + "_eig_freq.txt",data_all)
        np.savetxt("../data_out/" + self.project + "_eig_val.txt",data_eig_val)

    def post_process_pv(self):
        # Import the relavent data
        project = self.project + '_post_buckling'
        part_name = 'Merged'
        # Post-processing
        odb = openOdb(project + '.odb')

        step = odb.steps['Step-1']
        his_name = step.historyRegions.keys()[1]
        his_region = step.historyRegions[his_name]

        out_his_PCAV = his_region.historyOutputs['PCAV'].data
        out_his_CVOL = his_region.historyOutputs['CVOL'].data
        odb.close()

        pcav = np.asarray([out_his_PCAV[idx][1] for idx in range(len(out_his_PCAV))])
        cvol = np.asarray([out_his_CVOL[idx][1] for idx in range(len(out_his_CVOL))])

        data_all = np.array([cvol,pcav]).T
        np.savetxt("../data_out/" + self.project + "_pcav_cvol.txt",data_all)
    
    def post_process_twist(self):
        pass

    def post_process_centernodes(self):
        # Import the relavent data
        project = self.project + '_post_buckling'
        part_name = 'Merged'
        # Post-processing
        odb = openOdb(project + '.odb')
        session.viewports['Viewport: 1'].setValues(displayedObject=odb)
        cc = 0

        #calculate number of frames; defines size of output arrays
        num_frames = len(odb.steps['Step-1'].frames)
        time_all = np.zeros((num_frames))

        #find number of centernodes
        initial_coord_output = odb.steps['Step-1'].frames[0].fieldOutputs['COORD'].getSubset(region=
            odb.rootAssembly.instances[part_name.upper()+'-1'].nodeSets['CENTERNODES'], position=NODAL)
        num_nodes = len(initial_coord_output.values)

        #make data for x,y,z where a col is a timestep
        x_all = np.zeros((num_nodes,num_frames))
        y_all = np.zeros((num_nodes,num_frames))
        z_all = np.zeros((num_nodes,num_frames))

        for frame in odb.steps['Step-1'].frames:
            time_all[cc] = frame.frameValue
            fieldU = frame.fieldOutputs['COORD']
            
            ndFieldU = fieldU.getSubset(region=odb.rootAssembly.instances[part_name.upper()+'-1'].nodeSets['CENTERNODES'],
                position=NODAL)
            #below gives you the [x,y,z] coord of the 0th point
            #printAB(ndFieldU.values[0].data)
            nodal_coord = []
            for i,value in enumerate(ndFieldU.values):
                current_coord = value.data
                x_all[i,cc] = current_coord[0]
                y_all[i,cc] = current_coord[1]
                z_all[i,cc] = current_coord[2]
                #nodal_coord.append(value.data)

            cc = cc+1

        odb.close()
        np.savetxt("../data_out/" + self.project + "_x_all.txt",x_all)
        np.savetxt("../data_out/" + self.project + "_y_all.txt",y_all)
        np.savetxt("../data_out/" + self.project + "_z_all.txt",z_all)
        np.savetxt("../data_out/" + self.project + "_t_all.txt",time_all)
        
        print('Post-processing UC is done. File saved.')

    # def make_extra_buck(self,bdamp):
    #     Mdb()
    #     m = mdb.models['Model-1']
    #     self.model = m
    #     p_all = m.PartFromOdb(frame=22, instance='MERGED-1', name='MERGED-1', 
    #         odb=session.openOdb(self.project + '_post_buckling.odb'), shape=DEFORMED, step=0)

    #     self.add_materials()

    #     '''section bs part 2'''
    #     #todo: perhaps the thickness should come from geometry? rip lmaoooooo
    #     m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    #         integrationRule=SIMPSON, material='NH-1', name='Section-1', 
    #         nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    #         preIntegrate=OFF, temperature=GRADIENT, thickness=self.t1, thicknessField='', 
    #         thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    #     m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    #         integrationRule=SIMPSON, material='NH-2', name='Section-2', 
    #         nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    #         preIntegrate=OFF, temperature=GRADIENT, thickness=self.t2, thicknessField='', 
    #         thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
    #     m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    #         integrationRule=SIMPSON, material='NH-cap', name='Section-cap', 
    #         nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    #         preIntegrate=OFF, temperature=GRADIENT, thickness=self.w, thicknessField='', 
    #         thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

    #     p_all.SectionAssignment(offset=0.0, 
    #         offsetField='', offsetType=MIDDLE_SURFACE, region=
    #         p_all.sets['CAP_FACE'], sectionName=
    #         'Section-cap', thicknessAssignment=FROM_SECTION)
    #     p_all.SectionAssignment(offset=0.0, 
    #         offsetField='', offsetType=BOTTOM_SURFACE, region=
    #         p_all.sets['BODY-1'], sectionName=
    #         'Section-1', thicknessAssignment=FROM_SECTION)
    #     p_all.SectionAssignment(offset=0.0, 
    #         offsetField='', offsetType=BOTTOM_SURFACE, region=
    #         p_all.sets['BODY-2'], sectionName=
    #         'Section-2', thicknessAssignment=FROM_SECTION)

    #     printAB(p_all.sets['CAP_FACE'].elements)

    #     '''make assembly'''
    #     a = m.rootAssembly
    #     self.aa = a

    #     rp1=a.ReferencePoint(point=(0.0, 0.0, 0.0))
    #     rp2=a.ReferencePoint(point=(0.0, 0.0, 0.5*self.H))
    #     rp3=a.ReferencePoint(point=(0.0, 0.0, self.H))
    #     set_rp1 = a.Set(name='rp1', referencePoints=(a.referencePoints[rp1.id], ))
    #     set_rp2 = a.Set(name='rp2', referencePoints=(a.referencePoints[rp2.id], ))
    #     set_rp3 = a.Set(name='rp3', referencePoints=(a.referencePoints[rp3.id], ))

    #     #instance
    #     i_all = a.Instance(dependent=ON, name='MERGED-1-1', part=p_all)

    #     a.Surface(name='Surf-main', side12Elements=i_all.elements.getByBoundingCylinder(center1 = (0,0,-2*self.t1),
    #         center2 = (0,0,self.H + 2*self.t1), radius = self.R + 3*self.t1))

    #     m.FluidCavityProperty(bulkModulusTable=((2000.0, ), ), expansionTable=((1.0, ), ),
    #         fluidDensity=1e-09, name='IntProp-1', useBulkModulus=True, useExpansion=True)

    #     m.FluidCavity(cavityPoint=set_rp2, cavitySurface=a.surfaces['Surf-main'],
    #         createStepName='Initial', interactionProperty='IntProp-1', name='Int-1')

    #     '''AMPLITUDE'''
    #     m.SmoothStepAmplitude(data=((0.0, 0.0), (1.0, 1.0)), name='Amp-1', timeSpan=STEP)

    #     '''STEP'''
    #     m.BuckleStep(maxEigen=None, name='Step-1', numEigen=10, previous='Initial', vectors=18,maxIterations=500)

    #     '''BCS'''
    #     m.DisplacementBC(amplitude=UNSET, buckleCase=
    #         PERTURBATION_AND_BUCKLING, createStepName='Step-1', distributionType=
    #         UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-1', region=
    #         a.sets['MERGED-1-1.RING'], u1=0.0, u2=0.0, u3=0.0, 
    #         ur1=0.0, ur2=0.0, ur3=0.0)

    #     m.FluidCavityPressureBC(amplitude=UNSET, createStepName=
    #         'Step-1', fixed=OFF, fluidCavity='Int-3', magnitude=-1.0, name='BC-2')

    #     '''WRITE INP FILE'''
    #     # dir_path = os.path.dirname(os.path.realpath(__file__))
    #     m.keywordBlock.synchVersions(storeNodesAndElements=False)
    #     nlinesInput = len(m.keywordBlock.sieBlocks)
    #     m.keywordBlock.replace(nlinesInput-2, '\n*Output, field, variable=PRESELECT\n*NODE FILE\nU,')

    #     extra_str = '_lin_buckling_deformed'
    #     jname = self.make_job(extra_str)
    #     mdb.saveAs(self.project+extra_str)
    #     return jname

    def add_materials(self, nonlinear_model = False):
        m = self.model
        #Given material properties
        nu      = 0.5 #Poisson's ratio to account for incompressibility
        mu1     = self.E1/(2.*(1+nu)) #Initial shear modulus
        K1      = 1e9 #Bulk modulus very high K = E/3*(1-2nu)
        mu2     = self.E2/(2.*(1+nu)) #Initial shear modulus
        K2      = 1e9 #Bulk modulus very high K = E/3*(1-2nu)
        mu_cap  = self.E_cap/(2.*(1+nu)) #Initial shear modulus
        K_cap   = 1e9 #Bulk modulus very high K = E/3*(1-2nu)

        #Abaqus strain energy for neo-hokean form is
        # U = C_10 * (I_1 - 3) + 1/D_1(J^el-1)^2 (right part is zero for incomp.)
        # nh_d1_all = 0.1*0.0833 #(smol) david suggested this on 1/20/23

        nh_d1_all = 0.0 #(smol) david suggested this on 1/20/23

        nh_C10_1        = 0.5*mu1
        nh_D1_1         = nh_d1_all
        nh_C10_2        = 0.5*mu2
        nh_D1_2         = nh_d1_all
        nh_C10_cap      = 0.5*mu_cap
        nh_D1_cap       = nh_d1_all
        material_name_1 = 'NH-1'
        material_name_2 = 'NH-2'
        material_name_cap = 'NH-cap'

        '''MATERIALS'''
        #Material
        rho1 = 1e-09
        mat1 = m.Material(name=material_name_1);
        mat1.Density(table=((rho1,),))
        mat1.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_1,nh_D1_1), ), testData=OFF,
            type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

        rho2 = 1e-09
        mat2 = m.Material(name=material_name_2);
        mat2.Density(table=((rho2,),))
        mat2.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_2,nh_D1_2), ), testData=OFF,
            type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

        rho_cap = 1e-09
        mat_cap = m.Material(name=material_name_cap);
        mat_cap.Density(table=((rho_cap,),))
        mat_cap.Hyperelastic(materialType=ISOTROPIC, table=((nh_C10_cap,nh_D1_cap), ), testData=OFF,
            type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)

        self.mat_all = [mat1,mat2,mat_cap]

        if nonlinear_model:
            for i in range(len(self.mat_all)):
                self.mat_all[i].Damping(alpha=0., beta=self.bdamp)

    def mesh_part(self):
        p4 = self.part
        p4.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=self.h_element)

        if self.mesh_shape == 'tri':
            faces_all = p4.faces
            face_cap = FaceArray((p4.faces.findAt((0, 0, self.H)), ) )

            face_shell_only = [faces_all[i] for i in range(len(faces_all)) if faces_all[i] not in face_cap]
            face_shell_only = FaceArray(face_shell_only)

            p4.setMeshControls(elemShape=TRI, regions=face_shell_only, technique=STRUCTURED)
        else:
            p4.setMeshControls(elemShape=QUAD_DOMINATED, regions = p4.faces, technique=FREE)

        if self.mesh_order == 'quadratic':
            p4.setElementType(elemTypes=(ElemType(elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=STRI65, 
                elemLibrary=STANDARD)),  regions=(self.p3.faces, ))
        else:
            p4.setElementType(elemTypes=(ElemType(elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
                hourglassControl=DEFAULT), ElemType(elemCode=S3, elemLibrary=STANDARD)), regions=(self.p3.faces, ))
        p4.generateMesh()



        # if mesh_type == 'quad_lin':
        #         #todo: regions used to be p3.faces not p4 check if this is ok

        #     #p3.setMeshControls(elemShape=H,regions=p3.cells, technique=FREE)
        #     p4.setElementType(elemTypes=(ElemType(elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
        #         hourglassControl=DEFAULT), ElemType(elemCode=S3, elemLibrary=STANDARD)), regions=(self.p3.faces, ))
        # elif mesh_type == 'quad_quad':
        #     pass
        # elif mesh_type == 'tri_lin':
        #     # faces_all = p4.Set(name = 'faces_all', faces = p4.faces)
        #     # face_cap = p4.Set(name = 'face_cap', faces = FaceArray((p4.faces.findAt((0, 0, self.H)), ) ))

        #     faces_all = p4.faces
        #     face_cap = FaceArray((p4.faces.findAt((0, 0, self.H)), ) )
        #     printAB(p4.faces.findAt((0, 0, self.H)))
        #     # nodes_edges_only = [edge_nodes[i] for i in range(len(edge_nodes)) if edge_nodes[i] not in corner_nodes]

        #     face_shell_only = [faces_all[i] for i in range(len(faces_all)) if faces_all[i] not in face_cap]
        #     face_shell_only = FaceArray(face_shell_only)
        #     printAB(face_shell_only)

        #     p4.setMeshControls(elemShape=TRI, regions=face_shell_only, technique=STRUCTURED)

        # p4.generateMesh()



    def make_geometry(self, nonlinear_model = False):
        try:
            os.remove(self.project+'.lck')
            os.remove(self.project+'.odb')
        except OSError:
            pass


        #Abaqus related names
        model_name    = 'model'
        part_name     = 'part'
        assembly_name = 'instance'
        sketch_name   = 'sketch'
        cae_file      = self.project+'.cae'
        section_name  = 'solid_section'

        #Element selection 
        quad_element = 'S4R'#if shells are used
        h_element    =  self.h_element #scales with thickness

        #assign local var bc I copied this from david's code
        theta = self.theta
        H = self.H
        R = self.R
        t1 = self.t1
        t2 = self.t2
        w = self.w




        '''
        ------------------------------------------------------------------------
        Abaqus Pre-processing
        ------------------------------------------------------------------------
        '''

        '''MODEL'''
        Mdb()
        m = mdb.models['Model-1']
        self.model = m

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
        self.aa = a

        i1 = a.Instance(dependent=ON, name=assembly_name+'_1',part=p1)
        i2 = a.Instance(dependent=ON, name=assembly_name+'_2',part=p2)
        i3 = a.Instance(dependent=ON, name=assembly_name+'_3',part=p3)
        i4 = a.Instance(dependent=ON, name=assembly_name+'_4',part=p4)
        a.translate(instanceList=(assembly_name+'_4', ), 
        vector=(0.0, 0.0, H))
        a.InstanceFromBooleanMerge(domain=GEOMETRY,instances=tuple([a.instances[a.instances.keys()[xxx]] 
            for xxx in range(len(a.instances.keys()))]),
        keepIntersections=ON, name='Merged', originalInstances=SUPPRESS)
        p4 = m.parts['Merged']

        self.part = p4
        self.p3 = p3

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
        p4.Set(name='ring',edges=(p4.edges[e1.index:e1.index+1],
            p4.edges[e2.index:e2.index+1],p4.edges[e3.index:e3.index+1], ))
        p4.Set(name='curvature',edges=(p4.edges[e4.index:e4.index+1],))


        '''SECTIONS'''
        m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
            integrationRule=SIMPSON, material='NH-1', name='Section-1', 
            nodalThicknessField='', numIntPts=5, poissonDefinition=VALUE, poisson = self.nu_shell, 
            preIntegrate=OFF, temperature=GRADIENT, thickness=t1, thicknessField='', 
            thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
        m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
            integrationRule=SIMPSON, material='NH-2', name='Section-2', 
            nodalThicknessField='', numIntPts=5, poissonDefinition=VALUE, poisson = self.nu_shell, 
            preIntegrate=OFF, temperature=GRADIENT, thickness=t2, thicknessField='', 
            thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
        m.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
            integrationRule=SIMPSON, material='NH-cap', name='Section-cap', 
            nodalThicknessField='', numIntPts=5, poissonDefinition=VALUE, poisson = self.nu_shell, 
            preIntegrate=OFF, temperature=GRADIENT, thickness=w, thicknessField='', 
            thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

        p4.SectionAssignment(offset=0.0, offsetField='', offsetType=BOTTOM_SURFACE, region=p4.sets['body-1'],
            sectionName='Section-1', thicknessAssignment=FROM_SECTION)
        p4.SectionAssignment(offset=0.0, offsetField='', offsetType=BOTTOM_SURFACE, region=p4.sets['body-2'],
            sectionName='Section-2', thicknessAssignment=FROM_SECTION)
        p4.SectionAssignment(offset=0.0, offsetField='', offsetType=BOTTOM_SURFACE, region=p4.sets['cap_face'],
            sectionName='Section-cap', thicknessAssignment=FROM_SECTION)

        '''MESH'''
        self.mesh_part()

        '''SETS'''
        rp1=a.ReferencePoint(point=(0.0, 0.0, 0.0))
        rp2=a.ReferencePoint(point=(0.0, 0.0, 0.5*H))
        rp3=a.ReferencePoint(point=(0.0, 0.0, H))
        set_rp1 = a.Set(name='rp1', referencePoints=(a.referencePoints[rp1.id], ))
        set_rp2 = a.Set(name='rp2', referencePoints=(a.referencePoints[rp2.id], ))
        set_rp3 = a.Set(name='rp3', referencePoints=(a.referencePoints[rp3.id], ))

        '''MATERIALS'''
        self.add_materials(nonlinear_model)
        
        
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
        p4.Set(name='centernodes', nodes=p4.nodes.getByBoundingBox(-2*R,-2*R,0.5*H-0.9*h_element,
            2*R,2*R,0.5*H+0.1*h_element))

        m.FluidCavityProperty(bulkModulusTable=((2000.0, ), ), expansionTable=((1.0, ), ),
            fluidDensity=1e-09, name='IntProp-1', useBulkModulus=True, useExpansion=True)

        m.FluidCavity(cavityPoint=set_rp2, cavitySurface=a.surfaces['Surf-main'],
            createStepName='Initial', interactionProperty='IntProp-1', name='Int-1')

        


class beam_model:
    def __init__(self, project, imperfection = None, fullProps = None, simpProps = None, beamProps = None):
        self.project = project
        self.post_instability = True

        if fullProps is not None:
            self.R = float(fullProps.R)
            self.H = float(fullProps.H)
            self.w = float(fullProps.w)
            self.theta = float(fullProps.theta)

            self.t1 = float(fullProps.t1)
            self.t2 = float(fullProps.t2)

            self.E1 = float(fullProps.E1)
            self.E2 = float(fullProps.E2)
            self.E_cap = float(fullProps.E_cap)
        elif simpProps is not None:
            self.R = float(simpProps.R)
            self.H = float(simpProps.H)
            self.w = float(simpProps.w)
            self.theta = float(simpProps.theta)

            self.t1 = float(simpProps.t)
            self.t2 = self.t1

            self.E1 = float(simpProps.E)
            self.E2 = self.E1
            self.E_cap = self.E1
        else:
            raise ValueError("pls input something thx!!!")

        if beamProps is not None:
            self.base = float(beamProps.base)
            self.height = float(beamProps.height)

        if imperfection is not None:
            self.imperfection = imperfection * self.t1
        else:
            #set default value of 0.004
            self.imperfection = 0.004 * self.t1

        #finish any other initialization
        # original mult: 0.25 * 1.5/50
        self.h_element = 1.5/50 * self.H

    def run_linear_model(self):
        jname = self.make_linear_model()
        run_inp(jname)

    def make_linear_model(self, extra_springs = False):
        self.make_geometry(extra_springs = extra_springs)
        m = self.model

        #step
        m.BuckleStep(name='Step-1', numEigen=10, previous='Initial', vectors=18)
        m.ConcentratedForce(cf3=-1.0, createStepName='Step-1', distributionType=UNIFORM, field='',
            localCsys=None, name='Load-1', region=self.set_rp)

        '''WRITE INP FILE'''
        # dir_path = os.path.dirname(os.path.realpath(__file__))
        m.keywordBlock.synchVersions(storeNodesAndElements=False)
        nlinesInput = len(m.keywordBlock.sieBlocks)
        m.keywordBlock.replace(nlinesInput-2, '\n*Output, field, variable=PRESELECT\n*NODE FILE\nU,')

        extra_str = '_lin_buckling'
        jname = self.make_job(extra_str = extra_str)
        mdb.saveAs(self.project + extra_str)
        return jname

    def make_nonlin_model(self, bdamp, eig_idx):
        self.bdamp = bdamp
        self.make_geometry(nonlinear_model = True)
        imper = self.imperfection
        m = self.model
        aa = self.aa

        #step
        nincre = 5000
        m.ImplicitDynamicsStep(alpha=DEFAULT, amplitude=RAMP, 
           application=QUASI_STATIC, initialConditions=OFF, initialInc=0.001, maxInc=0.05,
           maxNumInc=nincre, minInc=1e-09, name='Step-1', nlgeom=ON, nohaf=OFF, previous='Initial')

        #bc force
        m.ConcentratedForce(cf3=-7.307, createStepName='Step-1', distributionType=UNIFORM, field='',
            localCsys=None, name='Load-1', region=self.set_rp)

        '''WRITE IMPERFECTION'''
        # asdsad
        m.keywordBlock.synchVersions(storeNodesAndElements=False)
        keyWords = m.keywordBlock.sieBlocks
        keyLine = '** ----------------------------------------------------------------\n** \n** STEP: Step-1\n** '
        idx_change = keyWords.index(keyLine)
        m.keywordBlock.replace(idx_change,
            '** ----------------------------------------------------------------\n*IMPERFECTION,FILE='+
            self.project+'_lin_buckling'+',STEP=1\t\n'+str(eig_idx)+','+str(imper)+'\n** \n** STEP: Step-1\n**')

        extra_str = '_post_buckling'
        jname = self.make_job(extra_str)
        mdb.saveAs(self.project+extra_str)
        return jname


    def make_geometry(self, nonlinear_model = False, extra_springs = False):
        Mdb()
        m = mdb.models['Model-1']
        self.model = m

        #discrete rigid top
        s1 = m.ConstrainedSketch(name='sketch_circ', sheetSize=200.0)
        s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(self.R, 0.0))
        p1 = m.Part(dimensionality=THREE_D, name='Part-cap', type=DISCRETE_RIGID_SURFACE)
        p1.BaseSolidExtrude(depth=self.w, sketch=s1)
        p1.RemoveCells(cellList=p1.cells)
        #note: extrudes in pos z direction

        #beam part
        p2 = m.Part(dimensionality=THREE_D, name='Part-beam', type=DEFORMABLE_BODY)
        yz_plane = p2.DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=YZPLANE)
        y_axis = p2.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
        transform_beam = p2.MakeSketchTransform(sketchPlane=p2.datums[yz_plane.id], sketchPlaneSide=SIDE1,
            sketchUpEdge=p2.datums[y_axis.id], sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
        s2 = m.ConstrainedSketch(name='sketch_beam', sheetSize=200.0, transform=transform_beam)
        s2.Line(point1=(0.0, 0.0), point2=(self.H, 0.0))
        
        p2.Wire(sketch=s2, sketchOrientation=RIGHT, sketchPlane=p2.datums[yz_plane.id],
            sketchPlaneSide=SIDE1, sketchUpEdge=p2.datums[y_axis.id])

        if extra_springs:
            p3 = m.Part(dimensionality=THREE_D, name='Part-spring', type=DEFORMABLE_BODY)
            s3 = m.ConstrainedSketch(name='sketch_bar', sheetSize=200.0)
            length_spring = self.R * np.sqrt(2) #only for 4 ridges obv
            s3.Line(point1=(0.0, 0.0), point2=(length_spring, 0.0))
            p3.BaseWire(sketch = s3)
        

        #assembly
        aa = m.rootAssembly
        self.aa = aa
        aa.DatumCsysByDefault(CARTESIAN)
        i_cap = aa.Instance(dependent=ON, name='Part-cap-1', part=p1)

        aa.Instance(dependent=ON, name='Part-beam-1', part=p2)
        aa.Instance(dependent=ON, name='Part-beam-2', part=p2)
        aa.Instance(dependent=ON, name='Part-beam-3', part=p2)
        aa.Instance(dependent=ON, name='Part-beam-4', part=p2)
        
        aa.translate(instanceList=('Part-beam-1', ), vector=(0.0, self.R, 0.0))
        aa.rotate(angle=90.0, axisDirection=(0.0, 0.0, -self.H), axisPoint=(0.0, self.R, 0.0),
            instanceList=('Part-beam-1', ))
        aa.translate(instanceList=('Part-beam-2', ), vector=(self.R, 0.0, 0.0))
        # aa.rotate(angle=90.0, axisDirection=(0.0, 0.0, -self.H), axisPoint=(self.R, 0.0, 0.0),
        #     instanceList=('Part-beam-2', ))
        aa.translate(instanceList=('Part-beam-3', ), vector=(0.0, -self.R, 0.0))
        aa.rotate(angle=90.0, axisDirection=(0.0, 0.0, -self.H), axisPoint=(0.0, -self.R, 0.0),
            instanceList=('Part-beam-3', ))
        aa.translate(instanceList=('Part-beam-4', ), vector=(-self.R, 0.0, 0.0))
        # aa.rotate(angle=90.0, axisDirection=(0.0, 0.0, -self.H), axisPoint=(-self.R, 0.0, 0.0),
        #     instanceList=('Part-beam-4', ))

        all_coord_top = [(0.0, self.R, 0.), (self.R, 0.0, 0.), (0.0, -self.R, 0.), (-self.R, 0.0, 0.)]

        if extra_springs:
            angle_rot = [45, 135, -135, -45]
            for i in range(4):
                idx = i + 1
                aa.Instance(dependent = ON, name = 'Part-spring-' + str(idx), part = p3)
                current_coord = (all_coord_top[i][0], all_coord_top[i][1], all_coord_top[i][2] - self.H/2)
                aa.translate(instanceList = ('Part-spring-' + str(idx), ), vector = current_coord)
                aa.rotate(angle = angle_rot[i], axisDirection = (0.0, 0.0, -self.H), axisPoint = current_coord,
                    instanceList = ('Part-spring-' + str(idx), ))


        #beam profile
        m.RectangularProfile(a=self.base, b=self.height, name='Profile-1')

        #material
        rho = 1e-9
        nu = 0.5
        mu = self.E1/(2.*(1 + nu))
        mat_nh = m.Material(name='Material-1')
        mat_nh.Hyperelastic(materialType=
            ISOTROPIC, table=((0.5*mu, 0.0), ), testData=OFF, type=NEO_HOOKE, 
            volumetricResponse=VOLUMETRIC_DATA)
        mat_nh.Density(table=((rho,),))

        if nonlinear_model:
            mat_nh.Damping(alpha = 0., beta = self.bdamp)


        m.BeamSection(consistentMassMatrix=False, integration=DURING_ANALYSIS, material='Material-1',
            name='Section-beam', poissonRatio=0.5, profile='Profile-1', temperatureVar=LINEAR)
        set_beam = p2.Set(edges=p2.edges, name='Set-beam')
        p2.assignBeamSectionOrientation(method=N1_COSINES, n1=(0.0, 1.0, 0.0), region=set_beam)
        p2.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=
            set_beam, sectionName='Section-beam', thicknessAssignment=FROM_SECTION)

        if extra_springs:
            set_spring = p3.Set(edges = p3.edges, name = 'Set-spring')
            area_spring = self.H * self.t1
            m.TrussSection(area = area_spring, material = 'Material-1', name='Section-spring')
            p3.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=set_spring,
            sectionName='Section-spring', thicknessAssignment=FROM_SECTION)

        #ref pt
        rp_cap = p1.ReferencePoint(point=(0.0, 0.0, self.w/2.))
        set_rp = p1.Set(name='set-rp', referencePoints=(p1.referencePoints[rp_cap.id], ))
        set_rp_inst = i_cap.sets['set-rp']
        self.set_rp = set_rp_inst

        #problems: regions for aa need to do reference that version of the set

        rho = 1e-9
        mass_cylinder = np.pi * self.R * self.R * self.w * rho
        I_zz = 0.5 * mass_cylinder * self.R * self.R
        I_xx = mass_cylinder * (0.25 * self.R * self.R + 1/3. * self.w * self.w)
        aa.engineeringFeatures.PointMassInertia(alpha=0.0, composite=0.0, mass=mass_cylinder,
            i11=I_xx, i12=0.0, i13=0.0, i22=I_xx, i23=0.0, i33=I_zz, name='Inertia-1', region=set_rp_inst)

        #mesh
        p1.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=2.8)
        p1.generateMesh()

        p2.setElementType(elemTypes=(ElemType(elemCode=B31H, elemLibrary=STANDARD), ), regions=set_beam)
        p2.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=self.H/10.)
        p2.generateMesh()

        if extra_springs:
            p3.seedPart(deviationFactor = 0.1, minSizeFactor = 0.1, size = self.R * np.sqrt(2))
            p3.setElementType(elemTypes=(ElemType(elemCode=T3D2, elemLibrary=STANDARD), ), regions=set_spring)
            p3.generateMesh()

        
        p2.Set(name = 'Set-center-node', nodes = p2.nodes.getByBoundingSphere(center = (0, 0, -self.H/2),
            radius = self.H/20))

        for i in range(4):
            idx = i + 1
            bar_1_idx = (i - 1) % 4 + 1
            bar_2_idx = i + 1

            set_name = 'Set-springs-node-' + str(idx)
            current_coord = (all_coord_top[i][0], all_coord_top[i][1], all_coord_top[i][2] - self.H/2)
            vertex_1 = aa.instances['Part-spring-' + str(bar_1_idx)].vertices.findAt(current_coord)
            vertex_2 = aa.instances['Part-spring-' + str(bar_2_idx)].vertices.findAt(current_coord)


            set_cur = aa.Set(name = set_name, vertices = VertexArray([vertex_1, vertex_2]))
            set_master = aa.instances['Part-beam-' + str(idx)].sets['Set-center-node']

            constraint_name = 'Constraint-spring-' + str(idx)
            m.Tie(adjust=OFF, master=set_master, name=constraint_name, positionTolerance=self.H/20,
                positionToleranceMethod=SPECIFIED, slave=set_cur, thickness=ON, tieRotations=ON)


        aa.regenerate()

        #BCs
        vertex_all = []
        vertex_all.append(aa.instances['Part-beam-1'].vertices.findAt((0.0, self.R, -self.H)))
        vertex_all.append(aa.instances['Part-beam-2'].vertices.findAt((self.R, 0.0, -self.H)))
        vertex_all.append(aa.instances['Part-beam-3'].vertices.findAt((0.0, -self.R, -self.H)))
        vertex_all.append(aa.instances['Part-beam-4'].vertices.findAt((-self.R, 0.0, -self.H)))
        vertex_all = VertexArray(vertex_all)

        set_beam_end = aa.Set(name='Set-end', vertices=vertex_all)

        

        m.DisplacementBC(amplitude=UNSET, buckleCase=PERTURBATION_AND_BUCKLING, createStepName='Initial',
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-1', region=set_beam_end,
            u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
        
        for i in range(4):
            idx = i + 1
            instance_name = 'Part-beam-' + str(int(idx))
            set_name = 's_Set-' + str(int(idx))
            constraint_name = 'Constraint-' + str(int(idx))
            vert_find = aa.instances[instance_name].vertices.findAt(all_coord_top[i])
            set_cur = aa.Set(name=set_name, vertices=VertexArray([vert_find]))
            m.Tie(adjust=OFF, master=set_rp_inst, name=constraint_name, positionTolerance=1.5 * self.R,
                positionToleranceMethod=SPECIFIED, slave=set_cur, thickness=ON, tieRotations=ON)

    def post_process_lin_buckle(self):
        odb = openOdb(path=self.project+'_lin_buckling.odb')
        step=odb.steps['Step-1']
        frames=step.frames
        initialframe = frames[1]
        f1=open(self.project+'_lin_buckling_eigenvalues.txt','w')
        f2=open(self.project+'_lin_buckling_displacement.txt','w')
        f3=open(self.project+'_lin_buckling_parameters.txt','w')
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
        f3.write(str(self.H)+' '+str(self.R)+' '+str(self.t1)+' '+str(self.t2)+' '+
            str(self.theta)+' '+str(self.w)+' '+str(self.E1)+' '+str(self.E2)+' '+str(self.E_cap)+' '+str(self.h_element))
        f1.close()
        f2.close()
        f3.close()
        
    def make_job(self,extra_str,num_threads = 1):
        m = self.model
        m.rootAssembly.regenerate()
        jname = self.project + extra_str
        j = mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name=jname, nodalOutputPrecision=SINGLE, 
            numCpus=num_threads, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
            ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
        j.writeInput()
        return jname

