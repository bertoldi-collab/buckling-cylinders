import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 800
bdamp = 0.0001
proj_name = '3fold-geo-test-'


#reduced multi buckle sweep
R_try = [10]
H_try = np.linspace(20,25,5)
t_try = np.linspace(1.1,1.2,5)



num_folds_all = np.zeros((len(R_try), len(H_try), len(t_try)))

for i in range(len(R_try)):
    R_cur = R_try[i]
    for j in range(len(H_try)):
        H_cur = H_try[j]
        for k in range(len(t_try)):
            printAB(idx_try)
            t_cur = t_try[k]
            geo_prop_cur = geoProps(R_cur, H_cur, 5, t_cur, 1.4) #all else is taken from geo_prop_three

            test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur, imperfection = 0.05)
            jname_lin = test.run_linear_model()
            jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.45, num_steps = 50)

            run_inp(jname_multi)
            test.post_process_multi_buckle()
            num_folds = test.post_process_num_folds()
            num_folds_all[i,j,k] = num_folds

            delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])
            delete_extra_files(jname_multi, ['.log', '.dat', '.msg'])

            idx_try += 1
# printAB(num_folds_all)
# np.savetxt('_3folds_num_folds_mode_1.txt', np.reshape(num_folds_all, (len(H_try), len(R_try)*len(t_try))))

#v600-799
# t_try = np.linspace(1.1,1.2,10)
# H_try = np.linspace(20,25,10)
# R_try = [8.8, 10]

#v500s:
# t_try = np.linspace(0.5,0.6,5)
# H_try = np.linspace(30,35,5)
# R_try = [8.8]

#original lin buckle sweep
# t_try = np.linspace(0.5,0.6,20)
# H_try = np.linspace(30,35,10)
# R_try = [8.8, 10]

#v100-v499: lin sweep of parameter space
#v500-524: multi sweep of more limited parameters


