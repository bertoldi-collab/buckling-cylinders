import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *

idx_try = 232
# test = full_3d(project = '2fold-3d-'+str(idx_try), simpProps = geo_prop_two)
# test.num_elem_thickness = 3

#### stuff for 3folds
test = full_shell(project = '3fold-imperfection-' + str(idx_try), simpProps = geo_prop_three, imperfection = 0.002)
test.stabilization_factor = 2e-8
jname_lin = test.run_linear_model()
jname_nonlin = test.make_nonlin_model(temp_mult = 0.15, is_buckling = True)

run_inp(jname_nonlin)

test.post_process_centernodes()

delete_extra_files(jname_lin)
delete_extra_files(jname_nonlin)

#### stuff for testing 3folds linear mode as a fxn of mesh size and shape ###
# mesh_mult = 0.25
# test = full_shell(project = '3fold-imperfection-' + str(idx_try), simpProps = geo_prop_three, imperfection = 0.002)

# # test.mesh_shape = 'tri'
# test.h_element = mesh_mult * test.h_element

# test.post_process_lin_centernodes(mode = 3)
# jname_lin = test.run_linear_model()
# num_folds = test.post_process_num_folds()
# printAB('num folds: ' + str(num_folds))

# delete_extra_files(jname_lin)
##########


##### stuff for testing dyn imp bonus bump #####
# test = full_shell(project = '4fold-imperfection-' + str(idx_try), simpProps = geo_prop_four, imperfection = 0.002)
# test = full_shell(project = '4fold-imperfection-' + str(idx_try), simpProps = geo_prop_four, imperfection = 0.05)

# test.stabilization_factor = 2e-8
# test.static_stable = False

# test.max_timestep_vol = 0.00002
# # test.adamp = 50

# jname_lin = test.run_linear_model()
# jname_nonlin = test.make_nonlin_model(temp_mult = 0.35, is_buckling = True)
# jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.35, num_steps = 100)

# run_inp(jname_nonlin)
# run_inp(jname_multi)
# test.post_process_pv()
# test.post_process_multi_pv()
# test.post_process_multi_buckle()
###########################

# jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.25, num_steps = 1)
# run_inp(jname_multi)

# test.post_process_multi_pv()
# test.add_cap_bc()
# test.regenerate_job(jname_multi)

# run_inp(jname_multi)
# test.post_process_multi_buckle()
# test.post_process_multi_pv()

# test.make_linear_model()
# jname_lin = test.run_linear_model()
# jname_riks = test.make_riks_model(pressure_set = -2)
# run_inp(jname_riks, num_threads = 2)
# jname_pressure = test.make_force_buckling_model(temp_mult = 0.2, pressure_app = 0.1)
# run_inp(jname_pressure)
# jname = test.make_nonlin_model(temp_set = -0.332*0.3, extra_imper = [(2,0.001), (3,0.002)])

# delete_extra_files(jname_lin)
# delete_extra_files(jname_pressure)
# delete_extra_files(jname_riks)

# run_inp(jname)
# test.post_process_twist()
# jname_post = test.make_nonlin_model(is_buckling = True, temp_set = -0.08)
# run_inp(jname_post,4)
# jname_post_2 = test.make_nonlin_model(temp_set = -0.332, alt_name = '_post_alt_eig', eig_name = '_post_buckling')
# jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.4)
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
#211: testing dyn --> pressure control on top face w/ (0.2, 0.1)
#212: testing multi buckling w/ UR3 = 0 on the cap (pensive emoji)
#213: try linear perturbation instead of freq for a single step [didn't work]
#214: single step 4folds up to 0.25 to compare pv
#215: dyn imp w/ new timestep algorithm to see if it helps 2nd buckle overshoot (imperfection = 0.05 oops)
#216: dyn imp w/ max_timestep_vol = 0.001 in cylinder_obj (imperfection = 0.002) temp_mult = 0.3 
#217: dyn imp w/ max_timestep_vol = 0.0002 in cylinder_obj (imperfection = 0.002) temp_mult = 0.3
#218: dyn imp w/ max_timestep_vol = 0.0002 in cylinder_obj (imperfection = 0.002) temp_mult = 0.3, adamp = 50 (for fun)
#219: dyn imp w/ max_timestep_vol = 0.00002 now class var (imper = 0.002, temp_mult = 0.3, adamp = 0)
#220: 3folds linear buckling quad mesh
#221: 3folds linear buckling tri mesh
#222: 3folds linear buckling quad mesh, mesh_mult = 0.25
#223: 3folds linear buckling tri mesh, mesh_mult = 0.25
#224: 4folds single static step (followed by freq) temp_mult = 0.35, inc_mult = 1e-2 (cylinder_obj), imper = 0.002, stable_fac = 2e-6
#225: 4folds single static step (followed by freq) temp_mult = 0.35, inc_mult = 1e-2 (cylinder_obj), imper = 0.05, stable_fac = 2e-6
#226: 4folds multi buckle, num_steps = 100, max_temp_mult = 0.35, imper = 0.05, stable_fac = 2e-6
#227: 4folds multi buckle, num_steps = 5, max_temp_mult = 0.35, imper = 0.05, stable_fac = 2e-6, manually deleting freq steps
#228: 4folds multi buckle, num_steps = 100, max_temp_mult = 0.35, imper = 0.05, stable_fac = 2e-6, manually deleting freq steps --> got same multi result
#229: 4folds single static step (followed by freq), temp mult = 0.35, imper = 0.05, no stable
#230: 4folds single static step (followed by freq), temp mult = 0.35, imper = 0.05, no stable, auto time max timestep = 0.05
#231: 4folds single static step (followed by freq), temp mult = 0.35, imper = 0.05, stable_fac = 2e-8, inc_mult = 1e-1 --> resolved multi behavior!!
#232: 3folds single static step (""), temp_mult = 0.15, imper = 0.002, stable_fac = 2e-8, extracting centernodes to test ability to resolve num folds