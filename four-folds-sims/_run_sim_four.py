import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 500
bdamp = 0.0001
proj_name = '4fold-fitting-'


#v500+ trial
E_try = [0.795]
t_try = [0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6]


for i in range(len(E_try)):
    E_cur = E_try[i]
    for j in range(len(t_try)):
        # printAB(idx_try)
        t_cur = t_try[j]
        geo_prop_cur = geoProps(10, 18, 5, t_cur, E_cur)

        test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur, imperfection = 0.002)
        test.transverse_shear = True
        jname_lin = test.run_linear_model()

        jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.5*-0.332)
        run_inp(jname_nonlin,4)

        test.post_process_pv()
        delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log'])
        delete_extra_files(jname_nonlin)

        idx_try += 1


#v100+ trial
# E_try = [1.0, 1.1, 1.15, 1.2, 1.25, 1.3]
# t_try = [0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56]

#v200+ trial
# E_try = [1.1, 1.2, 1.3, 1.4, 1.5]
# t_try = [0.50, 0.52, 0.54, 0.56, 0.58, 0.6]

#v300+ trial
# E_try = [1.4, 1.5]
# t_try = [0.52]

#v400+ trial
# E_try = [1.4]
# t_try = [0.56]