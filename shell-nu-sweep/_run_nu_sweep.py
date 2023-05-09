import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *
from time import time

idx_try = 520
bdamp = 0.0001
# proj_name = '2fold-test_nu-'
proj_name = '3fold-test_nu-'
# proj_name = 'bender-test_nu-'


# E_try = [1.0]
# t_try = [0.5]

nu_try = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
# time_all = [0, 0, 0]

for i in range(len(nu_try)):
    test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_three, imperfection = 0.004)
    test.nu_shell = nu_try[i]
    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.7*-0.332, eig_idx = 3)
    run_inp(jname_nonlin,4)

    test.post_process_pv()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
    delete_extra_files(jname_nonlin)

    idx_try += 1

# for i, nu in enumerate(nu_try):
#     test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_three, imperfection = 0.05)
#     # test = full_shell(project = proj_name+str(idx_try), fullProps = geo_prop_bend, imperfection = 0.05)
#     # test.mesh_shape = 'tri'
#     # test.mesh_order = 'quadratic'
#     test.nu_shell = nu
#     jname_lin = test.run_linear_model()
#     jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.45, num_steps = 40, eig_idx = 3)
#     # jname_riks = test.make_riks_model(bdamp, temp_mult = 0.3)

#     run_inp(jname_multi)
#     # run_inp(jname_riks)

#     test.post_process_multi_buckle()

#     delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
#     # delete_extra_files(jname_riks)
#     delete_extra_files(jname_multi)

#     idx_try += 1

#100s series: multi up to 0.245
#200s series: multi up to 0.26
#210s series: multi up to 0.28
#220s series: multi up to 0.28 and 30 steps
#230s series: multi up to 0.3 and 30 steps
#240s series: multi up to 0.35 and 30 steps
#250s series: multi up to 0.35 and 30 steps w/ tri elems
#[didn't work] 260s series: multi up to 0.35 and 30 steps w/ 2nd order quad elems
#300s series: 2 folds multi up to 0.35 and 30 steps
#400s series: bender up to 0.35 and 30 steps: they all failed early due to 1st eigenmode being the non bending one
#500s series: 3 folds multi up to 0.35 and 30 steps
#510s series: same as 500 but using eig idx of 3 to get 3 folds
#520s series: same as 510 but going up to 0.45 and w/ 40 steps
#530s series: same as 520 but w/ a larger imperfection




