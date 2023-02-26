import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *
from time import time

idx_try = 100
bdamp = 0.0001
proj_name = '4fold-test_w-'


# E_try = [1.0]
# t_try = [0.5]

w_try = [2.5, 5, 10]
time_all = [0, 0, 0]

for i in range(len(w_try)):
    geo_prop_cur = geoProps(10, 18, w_try[i], 0.54, 1.0)
    test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur, imperfection = 0.001)
    test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.5*-0.332)
    run_inp(jname_nonlin,4)

    test.post_process_pv()

    idx_try += 1



