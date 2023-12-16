import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

num_folds = 3

sim_type = str(num_folds)+'folds'

version = 100

E_all = [0.87, 1.2, 2.16]

w_all = 5

#R = 18 mm, H = 45 mm, t = 1.5 mm

for i, E in enumerate(E_all):
    geo_props_use = geoProps(18, 45, w_all, 1.5, E)
    proj_name = 'wrist-' + sim_type + '-' + str(version + i)

    test = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.002)
    test.stabilization_factor = 2e-8


    jname_lin = test.run_linear_model()
    num_folds_lin = test.post_process_num_folds()
    printAB(num_folds_lin)
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.6)
    # jname_nonlin = test.make_nonlin_model(temp_mult = final_temp_mult, is_buckling = True, eig_idx = 3)
    # jname_nonlin = test.make_static_dyn_model(temp_mult_static = 0.31, temp_mult_final = 0.7)
    run_inp(jname_nonlin)

    test.post_process_pv()

    # test.post_process_multi_pv(step_div_factor = 1, extra_str = '_post_buckling')
    # test.post_process_multi_contraction_twist(step_div_factor = 1, extra_str = '_post_buckling')

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
    delete_extra_files(jname_nonlin)