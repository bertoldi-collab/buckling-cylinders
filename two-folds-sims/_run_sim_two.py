import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 400
bdamp = 0.0001
proj_name = '2fold-fitting-'

#v100+ trial
# E_try = [1.0, 1.1, 1.15, 1.2, 1.25, 1.3]
# t_try = [1.13, 1.14, 1.15, 1.16, 1.17, 1.18]

#v200+ trial
# E_try = [1.1,1.2,1.3,1.4,1.5]
# t_try = [1.2,1.3,1.4,1.5,1.6]

#v300+ trial
# E_try = [1.4, 1.5]
# t_try = [1.2]

#v400+ trial
E_try = [1.2]
t_try = [1.3]




# E_try = [1.0]
# t_try = [0.5]

for i in range(len(E_try)):
    E_cur = E_try[i]
    for j in range(len(t_try)):
        # printAB(idx_try)
        t_cur = t_try[j]
        geo_prop_cur = geoProps(8.8, 44.75, 5, t_cur, E_cur)

        test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur, imperfection = 0.001)
        jname_lin = test.run_linear_model()
        jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.7*-0.332)
        run_inp(jname_nonlin,4)

        test.post_process_pv()

        delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log'])
        delete_extra_files(jname_nonlin)

        idx_try += 1


