import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *


def static_single(idx_try, R, H, t):
    proj_name = 'geo_sweep_v' + str(idx_try)

    w = 5.
    props_use = geoProps(R, H, w, t, 1.2)
    test = full_shell(project = proj_name, simpProps = props_use, imperfection = 0.002)
    test.stabilization_factor = 2e-8


    jname_lin = test.run_linear_model()
    jname_nonlin = test.make_nonlin_model(temp_mult = 0.4, is_buckling = True)

    run_inp(jname_nonlin)

    test.post_process_pv()
    test.post_process_contraction_twist()

    delete_extra_files(jname_lin, ['.fil', '.sta', '.log', '.dat', '.msg'])
    delete_extra_files(jname_nonlin)


def main():

    # idx_initial = 100
    idx_initial = 136
    # n1 = 3
    # n2 = 2
    # v1_sweep = np.linspace(0.1, 0.6, n1)
    # v2_sweep = np.linspace(0.01, 0.05, n2)

    n1 = 7
    n2 = 6
    # v1_sweep = np.linspace(0.2, 0.5, n1)
    v2_sweep = np.linspace(0.02, 0.04, n2)
    v1_sweep = [0.5]

    H_const = 20

    for i, v1 in enumerate(v1_sweep):
        for j, v2 in enumerate(v2_sweep):
            idx_try = idx_initial + i*len(v2_sweep) + j

            thickness = v2 * H_const
            R = v1 * H_const
            static_single(idx_try, R = R, H = H_const, t = thickness)


main()


#v10-15 was test sweep w/ v1_sweep = np.linspace(0.1, 0.6, 3), v2_sweep = np.linspace(0.01, 0.05, 2)


