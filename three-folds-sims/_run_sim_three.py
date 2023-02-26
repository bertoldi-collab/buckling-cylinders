import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 200
bdamp = 0.0001
proj_name = '3fold-fitting-'

#v100+ trial
# E_try = [0.9, 0.95, 1.0, 1.05, 1.1]
# t_try = [0.45, 0.475, 0.5, 0.525, 0.55]

#v200+ trial
E_try = [1.1, 1.15, 1.2, 1.25, 1.3]
t_try = [0.55, 0.575, 0.6, 0.625, 0.65]

# E_try = [1.0]
# t_try = [0.5]

for i in range(len(E_try)):
    E_cur = E_try[i]
    for j in range(len(t_try)):
        # printAB(idx_try)
        t_cur = t_try[j]
        geo_prop_cur = geoProps(8, 20, 5, t_cur, E_cur) #all else is taken from geo_prop_three

        test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur)
        test.run_linear_model()
        jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.5*-0.332)
        run_inp(jname_nonlin,4)

        test.post_process_pv()

        idx_try += 1


