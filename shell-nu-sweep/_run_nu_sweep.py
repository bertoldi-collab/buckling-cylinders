import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *
from time import time

idx_try = 710
bdamp = 0.0001

num_folds = 2
# proj_name = '2fold-test_nu-'
proj_name = str(num_folds) + 'fold-test_nu-'
# proj_name = 'bender-test_nu-'

if num_folds == 2: props_use = geo_prop_two
elif num_folds == 3: props_use = geo_prop_three
elif num_folds == 4: props_use = geo_prop_four
else: raise ValueError('yo')



# E_try = [1.0]
# t_try = [0.5]

nu_try = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
# nu_try = [0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49]


for i, nu in enumerate(nu_try):
    test = full_shell(project = proj_name+str(idx_try), simpProps = props_use, imperfection = 0.05)
    # test = full_shell(project = proj_name+str(idx_try), fullProps = geo_prop_bend, imperfection = 0.05)
    test.nu_shell = nu
    # test.static_stable = False
    jname_lin = test.run_linear_model()
    jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.45, num_steps = 50)

    run_inp(jname_multi)

    test.post_process_multi_buckle()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.odb'])
    delete_extra_files(jname_multi)

    idx_try += 1

#NEW RP
#700s series: 2 folds up to 0.45, 50 steps, static_stable = False
#710s series: 4 folds up to 0.45, 50 steps, static_stable = True

#800s series: 3 folds up to 0.45, 50 steps, eig_idx set manually, static_stable = False
#810s series: 4 folds up to 0.45, 50 steps, static_stable = True

#900s series: 4 folds up to 0.45, 50 steps, static_stable = False
#910s series: 4 folds up to 0.45, 50 steps, static_stable = True

#OLD RP
#100s series: multi up to 0.245
#200s series: multi up to 0.26
#210s series: multi up to 0.28
#220s series: multi up to 0.28 and 30 steps
#230s series: multi up to 0.3 and 30 steps
#240s series: multi up to 0.35 and 30 steps
#250s series: multi up to 0.35 and 30 steps w/ tri elems
#260s series: 4 folds up to 0.35 and 50 steps w/ damping factor to static step
#[OLD, deleted files] [didn't work] 260s series: multi up to 0.35 and 30 steps w/ 2nd order quad elems

#300s series: 2 folds multi up to 0.35 and 30 steps

#400s series: bender up to 0.35 and 30 steps: they all failed early due to 1st eigenmode being the non bending one

#500s series: 3 folds multi up to 0.35 and 30 steps
#510s series: same as 500 but using eig idx of 3 to get 3 folds
#520s series: same as 510 but going up to 0.45 and w/ 40 steps

#600s series: 3 folds, new geometric parameters (eig_idx still set manually)
#610s series: 3 folds, 30 steps instead of 40
#620s series: 3 folds, 40 steps, increase imperfection to 0.07
#630s series: 3 folds, 50 steps, imperfection back to 0.05
#640s series: 3 folds, 50 steps, added damping factor to static step (buck_cylinder_obj.py change)
#650s series: 3 folds, 50 steps, increase imperfection to 0.1 [lmao this one is literally the 4 folds]
#660s series: 3 folds, 50 steps, manually set eig_idx to force 3 folds
#670s series: 3 folds, 50 steps, using nu \in linspace(0.40, 0.49)



#measuring pressure as a fxn of nu:
# for i in range(len(nu_try)):
#     test = full_shell(project = proj_name+str(idx_try), simpProps = geo_prop_three, imperfection = 0.004)
#     test.nu_shell = nu_try[i]
#     jname_lin = test.run_linear_model()
#     jname_nonlin = test.make_nonlin_model(bdamp, temp_set = 0.7*-0.332, eig_idx = 3)
#     run_inp(jname_nonlin,4)

#     test.post_process_pv()

#     delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.odb'])
#     delete_extra_files(jname_nonlin)

#     idx_try += 1
