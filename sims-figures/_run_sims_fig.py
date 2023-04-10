import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

sim_type = '4folds'
geo_props_use = geo_prop_four
version = 104

final_temp_mult = 0.9
bdamp = 0.001
elem_size_mult = 0.5

proj_name = 'sim-long-' + sim_type + '-' + str(version)

test = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.001)

test.h_element = elem_size_mult * test.h_element
jname_lin = test.run_linear_model()
jname_nonlin = test.make_nonlin_model(bdamp, temp_set = final_temp_mult*-0.332)
run_inp(jname_nonlin,6)

test.post_process_pv()
delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log'])
delete_extra_files(jname_nonlin)

#FOUR FOLDS
#v100: ran normal, failed at 0.777: looks like we got local hourglassing that caused the sample to bend
#v101: trying 0.75*mesh size
#v102: trying 0.5*mesh size: showed more oscillation post 2nd instab
#v103: bumping bdamp from 0.0001 to 0.0005
#v104: bumping bdamp from 0.0005 to 0.001