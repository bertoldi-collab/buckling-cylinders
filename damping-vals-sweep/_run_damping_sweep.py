import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 100
bdamp = 0.0001

num_folds = 4
proj_name = str(num_folds) + 'fold-test_damping-'
# proj_name = 'bender-test_nu-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')

alpha_vals = np.linspace(0, 2, 11)

for i, adamp in enumerate(alpha_vals):
    test = full_shell(project = proj_name + str(idx_try + i), imperfection = 0.001, simpProps = props_use)
    test.adamp = adamp

    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.5)

    run_inp(jname_nonlin)

    test.post_process_pv()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log', '.dat', '.msg'])
    delete_extra_files(jname_nonlin)
