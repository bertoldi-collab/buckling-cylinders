import sys

sys.path.append('../')

# from buck_cylinder_obj import *
from cylinder_obj import *
from geo_prop import *

idx_try = 210
# idx_try = 207
bdamp = 0.0001
# jname_nonlin = '4fold-imperfection-'+str(idx_try) + '_post_buckling'
# test = full_3d(project = '2fold-3d-'+str(idx_try), simpProps = geo_prop_two)
# test.num_elem_thickness = 3
test = full_shell(project = '4fold-imperfection-' + str(idx_try), simpProps = geo_prop_four)
# test.post_process_pv()

# test.make_linear_model()
jname_lin = test.run_linear_model()
jname_riks = test.make_riks_model(bdamp, pressure_set = -2)
run_inp(jname_riks, num_threads = 2)
# jname = test.make_nonlin_model(bdamp, temp_set = -0.332*0.3, extra_imper = [(2,0.001), (3,0.002)])

delete_extra_files(jname_lin)
delete_extra_files(jname_riks)

# run_inp(jname)
# test.post_process_twist()
# jname_post = test.make_nonlin_model(bdamp, is_buckling = True, temp_set = -0.08)
# run_inp(jname_post,4)
# jname_post_2 = test.make_nonlin_model(bdamp, temp_set = -0.332, alt_name = '_post_alt_eig', eig_name = '_post_buckling')
# jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.4)
# run_inp(jname_multi)

# run_inp(jname_nonlin,4)
# jname_post_2 = '4fold-imperfection-'+str(idx_try) + '_post_alt_eig'
# run_inp(jname_post_2,8)

#NOTE TO FUTURE HELEN BEFORE YOU TRY TO RUN 118: WE NEED TO ADD THE STORAGE OF NODE OUTPUT FILES FOR THE POST BUCKLING ANALYSIS TO GRAB
#AND ALSO DO THE POST PROCESSING LINEAR BUCKLE PART FOR THE POST PROCESSSING OF NONLINEAR BUCKLING

# test.post_process_centernodes()

#v109 is attempt at 2 steps, did not work oops I'm overwriting it shrug
#new v109 is my attempt to check if the static step is fine --> it's fine
#v110 did not work uhh I'm lacking any bd conditions lmao
#static step ones used imperfection = 0.05
#v112 was -0.08 as temp_set and static followed by buckle [sim used for the I_x I_y stuff for 242r]
#v113 is trying to recreate take out whole vol: initially didn't work bc it bent for ??? reasons
#v114 is actually unused I'm just testing if my more general fxns for the nonliner static vs dynamic work
#v115 is testing the multi buckle thing; it has a bug where static steps don't inherit from prev static steps
#v116 is running 3 folds as long as it will go
#v117 is running 2 folds as long as it will go
#v118 is (a) remaking v112 (b) reading that eigenvalue in to do post buckling on the weird mode - NOT DONE
#204 was 3d testing
#205 is more 3d testing
#206: 3d 3folds
#207: 3d 2folds
#208: 3d testing partition on ring face (2folds)
#209/210: riks testing (not working I think)
