import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *


def damping_sweep(idx_try):
    num_samp = 5
    damping_sweep = 2*np.logspace(-8, -4, num_samp)

    for i, stab_factor in enumerate(damping_sweep):
        proj_name = 'bender-static-stable-v' + str(idx_try + i)
        test = full_shell(project = proj_name, fullProps = geo_prop_bend, imperfection = 0.1)

        test.stabilization_factor = stab_factor

        jname_lin = test.run_linear_model()
        jname_multi = test.make_nonlin_multi_buckle(max_temp_mult = 0.5, num_steps = 100, eig_idx = 2)

        run_inp(jname_multi)

        test.post_process_multi_pv()
        test.post_process_multi_buckle()

        delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])
        delete_extra_files(jname_multi)

def dyn_imp(idx_try):
    proj_name = 'bender-dyn-imp-v' + str(idx_try)

    deg_thin = 70
    rad_thick = (360 - deg_thin)*np.pi/180
    props_use = geoPropsFull(10, 18, 5, 0.856, 0.224, 1.2, 1.2, 1.2, rad_thick)
    test = full_shell(project = proj_name, fullProps = props_use, imperfection = 0.002)


    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.5)

    run_inp(jname_nonlin)

    test.post_process_pv()
    test.post_process_bend_angle()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])
    delete_extra_files(jname_nonlin)


def main():
    idx_try = 141
    dyn_imp(idx_try)

main()

#v140: dyn imp bender
#v141: dyn imp bender w/ new parameters, imperfection = 0.002
#v142: dyn imp bender w/ same prop as above, but 70 deg thin angle

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


