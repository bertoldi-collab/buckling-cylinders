import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 300
idx_last = 350



for i in range(idx_try, idx_last):
    odb_name = '2fold-test_nu-' + str(i) + '_multi_buckling.odb'
    proj_name = '2fold-test_nu-' + str(i)
    if os.path.exists(odb_name):
        instance = full_shell(proj_name, simpProps = geo_prop_two)
        instance.post_process_multi_buckle()