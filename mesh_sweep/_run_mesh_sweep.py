import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 105

num_folds = 4
proj_name = str(num_folds) + 'fold-test_mesh-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')

mesh_mult_all = [1, 0.5, 0.25]

for i, mesh_mult in enumerate(mesh_mult_all):
    test = full_shell(project = proj_name + str(idx_try + i), imperfection = 0.002, simpProps = props_use)
    test.h_element = mesh_mult * test.h_element
    test.stabilization_factor = 2e-6

    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.5)
    # jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.3, num_steps = 200)

    # run_inp(jname_multi)
    run_inp(jname_nonlin)

    test.post_process_pv()

    # test.post_process_multi_buckle()
    # test.post_process_multi_pv()
    
    delete_extra_files(jname_lin)
    delete_extra_files(jname_nonlin)
    # delete_extra_files(jname_multi)

#100-103: multi buckle mesh sweep, mesh_mult_all = [1, 0.5, 0.25, 0.125], imperfection = 0.05
#105-107: dyn imp mesh sweep, mesh_mult_all = [1, 0.5, 0.25], imperfection = 0.002