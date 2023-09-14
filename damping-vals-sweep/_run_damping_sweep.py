import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 140

num_folds = 4
proj_name = str(num_folds) + 'fold-test_damping-'
# proj_name = 'bender-test_nu-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')

# alpha_vals = np.linspace(0, 2, 11)
# imper_vals = 5*np.logspace(-4, -2, 3)
beta_vals = np.logspace(-8, -4, 5)


for i, bdamp in enumerate(beta_vals):
    test = full_shell(project = proj_name + str(idx_try + i), imperfection = 0.001, simpProps = props_use)
    # test.adamp = adamp
    test.bdamp = bdamp

    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.5)

    run_inp(jname_nonlin)

    test.post_process_pv()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.odb', '.log', '.dat', '.msg'])
    delete_extra_files(jname_nonlin)

#v100-110: damping sweep alpha_vals = np.linspace(0, 2, 11): no change in any of them (temp_mult = 0.5)
#v120: imperfection sweep imper_vals = 5*np.logspace(-4, -2, 3): still overshot for 2nd buckling (temp_mult = 0.5)
#v130: damping sweep beta_vals = np.logspace(-4, -1, 4), large change but none that look like static (temp_mult = 0.5)
#v140: damping sweep beta_vals = np.logspace(-8, -4, 5),