import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *


idx_try = 200


idx_try = 100
bdamp = 0.0001
proj_name = 'bender-test-H-v'

H_try = [20, 40, 60, 80, 100, 120, 140]


for i in range(len(H_try)):
    geo_prop_cur = geoPropsFull(10, H_try[i], 5, 0.75, 0.25, 1.4, 1.4, 1.4, 3*np.pi/2)

    test = full_shell(project = proj_name+str(idx_try), fullProps = geo_prop_cur, imperfection = 0.001)
    test.h_element = 1.5/50 * H_try[0]
    jname_lin = test.run_linear_model()


    delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])

    idx_try += 1


