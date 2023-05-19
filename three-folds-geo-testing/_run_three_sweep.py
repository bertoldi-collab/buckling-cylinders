import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 100
bdamp = 0.0001
proj_name = '3fold-geo-test-'


t_try = np.linspace(0.5,0.6,20)
H_try = np.linspace(30,35,10)
R_try = [8.8, 10]

for i in range(len(E_try)):
    E_cur = E_try[i]
    for j in range(len(t_try)):
        # printAB(idx_try)
        t_cur = t_try[j]
        geo_prop_cur = geoProps(10, 27, 5, t_cur, E_cur) #all else is taken from geo_prop_three

        test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_cur)
        jname_lin = test.run_linear_model()

        num_folds = test.post_process_num_folds()

        delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])

        idx_try += 1


