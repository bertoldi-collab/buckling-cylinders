import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 100

num_folds = 4
proj_name = str(num_folds) + 'fold-test_mesh-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')

mesh_mult_all = [1, 0.5, 0.25, 0.125]

for i, mesh_mult in enumerate(mesh_mult_all):
    test = full_shell(project = proj_name + str(idx_try + i), imperfection = 0.05, simpProps = props_use)
    test.h_element = mesh_mult * test.h_element

    jname_lin = test.run_linear_model()
    jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.3, num_steps = 200)

    run_inp(jname_multi)

    test.post_process_multi_buckle()
    test.post_process_multi_pv()
    