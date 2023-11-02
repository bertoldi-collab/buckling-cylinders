import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *
from time import time

idx_try = 520

num_folds = 4
proj_name = str(num_folds) + 'fold-test_static-stable-'
# proj_name = 'bender-test_nu-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')


# num_samp = 5
# damping_sweep = 2*np.logspace(-7, -4, num_samp)

# num_samp = 5
# damping_sweep = 2*np.logspace(-10, -6, num_samp)

# num_samp = 5
# damping_sweep = 2*np.logspace(-8, -4, num_samp)

num_samp = 7
damping_sweep = 2*np.logspace(-10, -4, num_samp)

for i, stab_factor in enumerate(damping_sweep):
    test = full_shell(project = proj_name+str(idx_try + i), simpProps = props_use, imperfection = 0.002)
    test.tangential_contact = True
    # test.theta = 4*np.pi/3
    test.stabilization_factor = stab_factor
    # test.nu_shell = 0.45
    # test = full_shell(project = proj_name+str(idx_try), fullProps = geo_prop_bend, imperfection = 0.05)
    jname_lin = test.run_linear_model()
    # jname_nonlin = test.make_nonlin_model(is_buckling = True, temp_mult = 0.35)
    # run_inp(jname_nonlin)

    jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.35, num_steps = 100)
    run_inp(jname_multi)

    # test.post_process_pv()
    # test.post_process_multi_buckle()
    test.post_process_multi_pv()
    # test.post_process_multi_contraction_twist()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.odb'])
    # delete_extra_files(jname_nonlin)
    delete_extra_files(jname_multi)

#v110s: 2folds up to 0.45, damping_sweep = 2*np.logspace(-7, -4, 5)
#v120s: 2folds up to 0.45, damping_sweep = 2*np.logspace(-10, -6, 5)
#v130s: 2folds up to 0.65, damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 100 (no tangential contact) (failed a little before 0.5)
#v140s: 2folds up to 0.50, damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 100, tangential contact [re-exported]

#v220s: 3folds up to 0.45 (eig_idx passed in as 3), damping_sweep = 2*np.logspace(-7, -4, 5)
#v230s: 3folds up to 0.45 (eig_idx passed in as 3), damping_sweep = 2*np.logspace(-10, -6, 5)
#v240s: 3folds up to 0.45 (eig_idx passed in as 3), damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 200, theta = 4*pi/3
#v250s: 3folds up to 0.45 (eig_idx passed in as 3), damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 200
#v260s: 3folds up to 0.45 (eig_idx passed in as 3), damping_sweep = 2*np.logspace(-8, -4, 5), num_steps = 400 [re-exported]

#v330s: 4folds up to 0.45, damping_sweep = 2*np.logspace(-7, -4, 5)
#v340s: 4folds up to 0.45, damping_sweep = 2*np.logspace(-10, -6, 5)
#v350s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-10, -6, 5)
#v360s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 100
#v370s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 200
#v380s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-7, -4, 5), num_steps = 100
#v390s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-10, -6, 5), num_steps = 400
#v400s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-8, -4, 5), num_steps = 800 [re-exported] (0.25 is wrong right? it's like 0.45 right? + a random number of steps bc conv)

#v500s: 4folds up to 0.25, damping_sweep = 2*np.logspace(-10, -4, 7), static single, imper = 0.002
#v510s: 4folds up to 0.35, damping_sweep = 2*np.logspace(-10, -4, 7), static single, imper = 0.002
#v520s: 4folds up to 0.35, damping_sweep = 2*np.logspace(-10, -4, 7), num_steps = 100, imper = 0.002

#v100s: 2folds up to 0.55
#v200s: 3folds up to 0.55 (eig_idx manually passed in as 3)
#v210s: 3folds up to 0.45 + nu = 0.45 (eig_idx passed in)
#v300s: 4folds up to 0.55
#v310s: 4folds up to 0.30
#v320s: 4folds up to 0.30 + nu = 0.45