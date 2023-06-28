import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *
from time import time

idx_try = 100
bdamp = 0.0001

num_folds = 4
proj_name = str(num_folds) + 'fold-test_3d-v'
# proj_name = 'bender-test_nu-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')


test = full_3d(project = proj_name + str(idx_try), simpProps = props_use)

# test.make_linear_model()
jname_lin = test.run_linear_model()
jname_nonlin = test.make_nonlin_model(bdamp, temp_set = -0.332*0.7)
run_inp(jname_nonlin)
test.post_process_pv()

delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
delete_extra_files(jname_nonlin)

#100: 4folds (shows 4 folds for lin but 5 folds for 3d)
#101: 3folds (eig idx auto)
#102: 2folds