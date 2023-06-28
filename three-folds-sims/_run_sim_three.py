import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 500
bdamp = 0.0001
# proj_name = '3fold-fitting-'
proj_name = '3fold-multi-imper'

imper_second = np.linspace(0.0005, 0.01, 10)

for i, imper in enumerate(imper_second):
    idx_cur = idx_try + i
    test = full_shell(project = proj_name+str(idx_cur), simpProps = geo_prop_three)

    extra_imper = [(1,imper)]

    jname_lin = test.run_linear_model()
    jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.55, num_steps = 50,
                    eig_idx = 3, extra_imper = extra_imper)
    
    run_inp(jname_multi)
    test.post_process_multi_buckle()


#v100+ trial
# E_try = [0.9, 0.95, 1.0, 1.05, 1.1]
# t_try = [0.45, 0.475, 0.5, 0.525, 0.55]

#v200+ trial
# E_try = [1.1, 1.15, 1.2, 1.25, 1.3]
# t_try = [0.55, 0.575, 0.6, 0.625, 0.65]

#v300+ trial
# E_try = [1.4]
# t_try = [0.45, 0.5, 0.55, 0.60, 0.65]

#v310+ trial
# E_try = [1.4]
# t_try = [0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62]

#v320+
# E_try = [1.4]
# t_try = [0.59, 0.60]

#v400+
# E_try = [1.4]
# t_try = [0.54,0.56,0.58,0.60,0.62,0.64]

#v410+
# E_try = [1.4]
# t_try = [0.57,0.58]


# for i in range(len(E_try)):
#     E_cur = E_try[i]
#     for j in range(len(t_try)):
#         # printAB(idx_try)
#         t_cur = t_try[j]
#         geo_prop_cur = geoProps(10, 27, 5, t_cur, E_cur) #all else is taken from geo_prop_three

#         test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur)
        # jname_lin = test.run_linear_model()
#         jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.7*-0.332)
#         run_inp(jname_nonlin,4)

#         test.post_process_pv()

#         delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
#         delete_extra_files(jname_nonlin)

#         idx_try += 1


