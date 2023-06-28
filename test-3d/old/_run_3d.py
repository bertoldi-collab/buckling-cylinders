import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 100
bdamp = 0.0001
jname_nonlin = '3fold-brick-elem-'+str(idx_try) + '_post_buckling'
test = full_shell(project = '3fold-brick-elem-'+str(idx_try), simpProps = geo_prop_three, shell_mesh = False)

test.run_linear_model()
# test.make_nonlin_model(bdamp)

#todo: having issues w/ assigning region for meshing 3d elements
