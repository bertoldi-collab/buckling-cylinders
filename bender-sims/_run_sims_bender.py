import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *



def damping_sweep(idx_try):
    num_samp = 5
    damping_sweep = 2*np.logspace(-8, -4, num_samp)

    for i, stab_factor in enumerate(damping_sweep):
        proj_name = 'bender-static-stable-v' + str(idx_try + i)
        test = full_shell(project = proj_name, fullProps = geo_prop_bend, imperfection = 0.05)

        test.stabilization_factor = stab_factor

        # jname_lin = test.run_linear_model()
        # jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.5, num_steps = 100)

        # run_inp(jname_multi)

        # test.post_process_multi_pv()
        test.post_process_multi_buckle()

        # delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])
        # delete_extra_files(jname_multi)

def dyn_imp(idx_try, deg_thin = 90, t_thin = 0.224, t_thick = 0.856):
    proj_name = 'bender-dyn-imp-v' + str(idx_try)

    rad_thick = (360 - deg_thin)*np.pi/180

    props_use = geoPropsFull(10, 18, 5, t_thick, t_thin, 1.2, 1.2, 1.2, rad_thick)
    test = full_shell(project = proj_name, fullProps = props_use, imperfection = 0.002)


    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.5)

    run_inp(jname_nonlin)

    test.post_process_pv()
    test.post_process_bend_angle()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])
    delete_extra_files(jname_nonlin)

def extract_centernodes(idx_try):
    proj_name = 'bender-dyn-imp-v' + str(idx_try)
    props_use = geo_prop_bend
    test = full_shell(project = proj_name, fullProps = props_use, imperfection = 0.002)

    test.post_process_lin_centernodes()




def main():
    # idx_try = 141
    idx_initial = 200
    theta_sweep = np.linspace(30, 150, 7)
    thickness_ratio_sweep = np.linspace(0.1, 0.9, 9)

    t_thick_const = 0.85

    for i, theta in enumerate(theta_sweep):
        for j, t_ratio in enumerate(thickness_ratio_sweep):
            idx_try = idx_initial + i*len(thickness_ratio_sweep) + j

            t_thin_cur = t_ratio * t_thick_const
            dyn_imp(idx_try, deg_thin = theta, t_thin = t_thin_cur, t_thick = t_thick_const)

def main_single():
    idx_try = 143
    t_thick = 0.92
    t_thin = 0.23

    dyn_imp(idx_try, t_thin = t_thin, t_thick = t_thick)
    
def main_damping():
    idx_try = 150
    damping_sweep(idx_try)

main_damping()

# extract_centernodes(141)

# main_damping()

# main()
# main_single()

#v140: dyn imp bender
#v141: dyn imp bender w/ new parameters, imperfection = 0.002
#v142: dyn imp bender w/ same prop as above, but 70 deg thin angle

#v200s: dyn imp sweep


#idk what this stuff below is, but it maybe got overwritten oops (might also be pre fixing the rp tho so)
#v200s: bender multi, damping_sweep = 2*np.logspace(-8, -4, num_samp), steps = 100, max_temp_mult = 0.5, imperfection = 0.05
#v210s: bender multi, damping_sweep = 2*np.logspace(-8, -4, num_samp), steps = 100, max_temp_mult = 0.5, imperfection = 0.1, eig_idx = 2

# idx_try = 100
# bdamp = 0.0001
# proj_name = 'bender-test-H-v'

# H_try = [20, 40, 60, 80, 100, 120, 140]


# for i in range(len(H_try)):
#     geo_prop_cur = geoPropsFull(10, H_try[i], 5, 0.75, 0.25, 1.4, 1.4, 1.4, 3*np.pi/2)

#     test = full_shell(project = proj_name+str(idx_try), fullProps = geo_prop_cur, imperfection = 0.001)
#     test.h_element = 1.5/50 * H_try[0]
#     jname_lin = test.run_linear_model()


#     delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])

#     idx_try += 1


