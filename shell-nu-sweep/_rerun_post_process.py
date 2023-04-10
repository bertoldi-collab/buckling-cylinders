import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 100
idx_last = 261



for i in range(idx_try, idx_last):
    odb_name = '4fold-test_nu-' + str(i) + '_multi_buckling.odb'
    proj_name = '4fold-test_nu-' + str(i)
    if os.path.exists(odb_name):
        instance = full_shell(proj_name, simpProps = geo_prop_four)
        instance.post_process_multi_buckle()