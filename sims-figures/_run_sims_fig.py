import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

sim_type = '2folds'
geo_props_use = geo_prop_two
version = 111

final_temp_mult = 0.75
bdamp = 0.001
elem_size_mult = 0.5

proj_name = 'sim-long-' + sim_type + '-' + str(version)

test = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.05)
# test.static_stable = False

# test.h_element = elem_size_mult * test.h_element
jname_lin = test.run_linear_model()
# jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.55, num_steps = 50)
# jname_nonlin = test.make_nonlin_model(bdamp, temp_set = final_temp_mult*-0.332)
# run_inp(jname_nonlin,6)
# run_inp(jname_multi)
# test.post_process_multi_buckle()

# test.post_process_pv()
# test.post_process_contraction_twist()
delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log'])
# delete_extra_files(jname_multi)
# delete_extra_files(jname_nonlin)

# num_folds = test.post_process_num_folds()
# printAB(num_folds)

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


