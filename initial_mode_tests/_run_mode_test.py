import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_initial = 100
num_try_t = 25
num_try_R = 35
# num_try_t = 2
# num_try_R = 3

v_R_try = np.linspace(0.1,0.6, num_try_R) #R/H
v_t_try = np.linspace(0.01,0.05, num_try_t) #t/H

H_set = 20 #standardizes the mesh size as a bonus
E_set = 1.2
w_set = 5

n_all = np.zeros((num_try_R, num_try_t), dtype = int)
n_all_dyn = np.zeros((num_try_R, num_try_t), dtype = int)

for i, v_R in enumerate(v_R_try):
    R = H_set * v_R
    for j, v_t in enumerate(v_t_try):
        t = H_set * v_t
        idx_cur = idx_initial + len(v_t_try)*i + j
        props_use = geoProps(R, H_set, w_set, t, E_set)

        proj_name = 'dyn_mode_sweep_v' + str(idx_cur)
        # proj_name = 'initial_mode_sweep_v' + str(idx_cur)
        test = full_shell(project = proj_name, simpProps = props_use, imperfection = 0.002)
        test.max_timestep_vol = 0.005

        jname_lin = test.run_linear_model()
        num_folds = test.post_process_num_folds()

        jname_nonlin = test.make_nonlin_model(temp_mult = 0.14)
        run_inp(jname_nonlin)
        num_folds_dyn = test.post_process_num_folds(extra_str = '_post_buckling', frame_sp = -1)

        n_all[i,j] = num_folds
        n_all_dyn[i,j] = num_folds_dyn
        delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg', '.odb'])
        delete_extra_files(jname_nonlin, ['.sta', '.log', '.dat', '.msg'])

printAB(n_all)
printAB('----------------------')
printAB(n_all_dyn)
# np.savetxt('../data_out/initial_mode_tests.txt', n_all, fmt = '%i')
np.savetxt('../data_out/dyn_mode_tests_lin.txt', n_all, fmt = '%i')
np.savetxt('../data_out/dyn_mode_tests_nonlin.txt', n_all_dyn, fmt = '%i')
# np.savetxt('../data_out/test_initial_10.txt', n_all, fmt = '%i')
# np.savetxt('../data_out/test_dyn_10.txt', n_all_dyn, fmt = '%i')
# np.savetxt('../data_out/test_matrix_10.txt', n_all, fmt='%i')


#v10-15 test sweep: num_try_t = 2, num_try_R = 3
#v100 sweep: num_try_t = 25, num_try_R = 35