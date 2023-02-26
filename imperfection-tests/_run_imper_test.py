import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 105
bdamp = 0.0001

imperfection_test = [0.0, 0.001, 0.002, 0.003, 0.004]

for i in range(len(imperfection_test)):
    test = full_shell(project = 'imper-test-4fold-v'+str(idx_try), simpProps = geo_prop_four,
        imperfection = imperfection_test[i])
    jname_lin = test.run_linear_model()

    # jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.016, num_steps = 10)
    jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.5*-0.332)
    run_inp(jname_nonlin)
    test.post_process_pv()

    delete_extra_files(jname_lin)
    delete_extra_files(jname_nonlin)

    idx_try += 1

#100-104: trying normal

