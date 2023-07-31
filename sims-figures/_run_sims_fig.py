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

version = 203

final_temp_mult = 0.75
bdamp = 0.001
# elem_size_mult = 0.5

proj_name = 'sim-long-' + sim_type + '-' + str(version)

test = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.002)
test.tangential_contact = True
# test.static_stable = False

# test.h_element = elem_size_mult * test.h_element
jname_lin = test.run_linear_model()
# jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.55, num_steps = 50)
jname_nonlin = test.make_nonlin_model(bdamp, temp_set = final_temp_mult*-0.332)
run_inp(jname_nonlin,2)
# run_inp(jname_multi)
# test.post_process_multi_buckle()

test.post_process_pv()
test.post_process_contraction_twist()
test.post_process_centernodes()
delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log'])
# delete_extra_files(jname_multi)
delete_extra_files(jname_nonlin)

# num_folds = test.post_process_num_folds()
# printAB(num_folds)

#new rp
#v200: 2folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002)
#v201: 3folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002)
#v202: 4folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002)
#v203: 2folds up to 0.75, extracting pv/contraction/twist/centernodes (imperfection 0.002) + adding tangential friction

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


