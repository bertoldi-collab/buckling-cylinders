import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

num_folds = 2
if num_folds == 2: geo_props_use = geo_prop_two
elif num_folds == 3: geo_props_use = geo_prop_three
elif num_folds == 4: geo_props_use = geo_prop_four
else: raise ValueError('yo')

# sim_type = '2folds'
sim_type = str(num_folds)+'folds'

version = 221

# final_temp_mult = 0.5
elem_size_mult = 1.0

proj_name = 'sim-long-' + sim_type + '-' + str(version)

test = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.002)
test.stabilization_factor = 2e-8
# test.tangential_contact = True
# test.E_cap = 4e3
# test.theta = 4*np.pi/3
# test.static_stable = False

test.h_element = elem_size_mult * test.h_element
jname_lin = test.run_linear_model()
# jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.55, num_steps = 50)
# jname_nonlin = test.make_nonlin_model(temp_mult = final_temp_mult, is_buckling = True, eig_idx = 3)
jname_nonlin = test.make_static_dyn_model(temp_mult_static = 0.31, temp_mult_final = 0.7)
run_inp(jname_nonlin)
# # run_inp(jname_multi)
# # test.post_process_multi_buckle()

test.post_process_multi_pv(step_div_factor = 1, extra_str = '_post_buckling')
test.post_process_multi_contraction_twist(step_div_factor = 1, extra_str = '_post_buckling')
# test.post_process_centernodes()
# test.post_process_lin_centernodes(mode = 3)
delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
# # delete_extra_files(jname_multi)
delete_extra_files(jname_nonlin)

# num_folds = test.post_process_num_folds()
# printAB(num_folds)

#new rp
#v200: 2folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002)
#v201: 3folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002)
#[OLD used figure making] v202: 4folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002)
#[OLD used figure making] v203: 2folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + adding tangential friction
#[OLD used figure making] v204: 3folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + mesh_size/2
#v205: 4folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + mesh_size/2
#v206: 2folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + mesh_size/2 + adding tangential friction + E_cap = 4e3 [MPa] [worse fit to data]
#v207: [FORMED 4 FOLDS] 3folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + mesh_size/2 + E_cap = 4e3 [MPa] + theta = 4*pi/3
#210: 4folds just lin, mesh_mult = 0.25
#211: 3folds just lin, mesh_mult = 0.25
#212: 2folds just lin, mesh_mult = 0.25

#213: 4folds, 1 static step, mesh_mult = 0.5, temp_mult = 0.6, stable_fac = 2e-8, imper = 0.002, extracting pv/cont/twist/centernodes [did not finish]
#214: 3folds, 1 static step, mesh_mult = 0.5, temp_mult = 0.6, stable_fac = 2e-8, imper = 0.002, extracting pv/cont/twist/centernodes, eig_idx = 3
#[used figure ses] 215: 2folds, 1 static step, mesh_mult = 0.5, temp_mult = 0.6, stable_fac = 2e-8, imper = 0.002, extracting pv/cont/twist/centernodes
#[used figure ses] 216: 4folds, 1 static step, mesh_mult = 1.0, temp_mult = 0.5, stable_fac = 2e-8, imper = 0.002, extracting pv/cont/twist/centernodes
#[used figure ses] 217: 3folds, 1 static step, mesh_mult = 1.0, temp_mult = 0.5, stable_fac = 2e-8, imper = 0.005, extracting pv/cont/twist/centernodes
#218: 3folds, 1 static step, mesh_mult = 1.0, temp_mult = 0.5, stable_fac = 2e-8, imper = 0.005, extracting pv/cont/twist/centernodes, eig_idx = 3

#219: 4folds, static/dyn 0.22/0.7 (stable_fac = 2e-8, imper = 0.005)
#220: 3folds, static/dyn 0.31/0.7 (stable fac = 2e-8, imper = 0.002)
#221: 2folds, static/dyn 0.31/0.7 (stable fac = 2e-8, imper = 0.002)

#some minimal tests
#v208: 3folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + mesh_size/2 +  theta = 4*pi/3
#v209: 3folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + mesh_size/2 +  theta = default value (iirc pi)

#########################################
#OLD RP
#FOUR FOLDS
#v100: ran normal, failed at 0.777: looks like we got local hourglassing that caused the sample to bend
#v101: trying 0.75*mesh size
#v102: trying 0.5*mesh size: showed more oscillation post 2nd instab
#v103: bumping bdamp from 0.0001 to 0.0005
#v104: bumping bdamp from 0.0005 to 0.001

#LINEAR BUCKLING ONLY
#v105: just linear buckling 4folds
#v106: just linear buckling 3folds
#v107: just linear buckling 2folds

#long 2/3/4
#v108: 3folds up to 0.75
#v109: 4folds up to 0.75
#v110: 2folds up to 0.75 + multi buckling (up to 0.55)

#long testing 2
#v111: 2folds mutli up to 0.55 w/ static stable

#long 2/3/4 [changed rp location]
#v112: 3folds up to 0.75 (eq of v108)
#v113: 4folds up to 0.75 (eq of v109) (did not finish)
#v114: 2folds up to 0.75 (eq of v110)


