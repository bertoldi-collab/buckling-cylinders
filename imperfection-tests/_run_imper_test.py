import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 120

# imperfection_test = [0.0, 0.001, 0.002, 0.003, 0.004]
imperfection_test = [0.0, 0.001, 0.005, 0.01, 0.05]

for i in range(len(imperfection_test)):
    test = full_shell(project = 'imper-test-4fold-v'+str(idx_try + i), simpProps = geo_prop_four,
        imperfection = imperfection_test[i])
    test.stabilization_factor = 2e-8
    jname_lin = test.run_linear_model()

    # jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.016, num_steps = 10)
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.3, is_buckling = True)
    run_inp(jname_nonlin)
    test.post_process_pv()

    delete_extra_files(jname_lin)
    delete_extra_files(jname_nonlin)

#100-104: trying normal, imperfection_test = [0.0, 0.001, 0.002, 0.003, 0.004]
#v110+: imperfection_test = [0.001, 0.005, 0.01, 0.05], static
#v120+: imperfection_test = [0.0, 0.001, 0.005, 0.01, 0.05] static, 2e-8, temp_mult = 0.3


