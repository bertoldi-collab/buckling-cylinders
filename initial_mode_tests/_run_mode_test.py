import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_initial = 10
# num_try_t = 25
# num_try_R = 35
num_try_t = 2
num_try_R = 3

v_R_try = np.linspace(0.1,0.6, num_try_R) #R/H
v_t_try = np.linspace(0.01,0.05, num_try_t) #t/H

H_set = 20 #standardizes the mesh size as a bonus
E_set = 1.2
w_set = 5

n_all = np.zeros((num_try_R, num_try_t), dtype = int)

for i, v_R in enumerate(v_R_try):
    R = H_set * v_R
    for j, v_t in enumerate(v_t_try):
        t = H_set * v_t
        idx_cur = idx_initial + len(v_t_try)*i + j
        props_use = geoProps(R, H_set, w_set, t, E_set)

        proj_name = 'initial_mode_sweep_v' + str(idx_cur)
        test = full_shell(project = proj_name, simpProps = props_use)

        jname_lin = test.run_linear_model()
        num_folds = test.post_process_num_folds()

        n_all[i,j] = num_folds
        delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg', '.odb'])

printAB(n_all)
np.savetxt('../data_out/initial_mode_tests.txt', n_all, fmt = '%i')
# np.savetxt('../data_out/test_matrix_10.txt', n_all, fmt='%i')


#v10-15 test sweep: num_try_t = 2, num_try_R = 3
#v100 sweep: num_try_t = 25, num_try_R = 35