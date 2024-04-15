import sys

sys.path.append('../')

from cylinder_obj import *
from geo_prop import *



def run_wrist(E_all, w_all, version, static_end = None):
    num_folds = 3
    sim_type = str(num_folds)+'folds'
    #R = 18 mm, H = 45 mm, t = 1.5 mm
    if static_end is None: static_end = [0.42]*len(E_all)


    for i, E in enumerate(E_all):
        geo_props_use = geoProps(18, 45, w_all, 1.5, E)
        proj_name = 'wrist-' + sim_type + '-' + str(version + i)

        test = full_shell(project = proj_name, simpProps = geo_props_use, imperfection = 0.002)
        test.stabilization_factor = 2e-8
        test.h_element = test.h_element / 2.


        jname_lin = test.run_linear_model()
        num_folds_lin = test.post_process_num_folds()
        printAB(num_folds_lin)
        jname_nonlin = test.make_nonlin_model(temp_mult = 0.6)
        # jname_nonlin = test.make_nonlin_model(temp_mult = final_temp_mult, is_buckling = True, eig_idx = 3)
        # jname_nonlin = test.make_static_dyn_model(temp_mult_static = static_end[i], temp_mult_final = 0.6)
        run_inp(jname_nonlin)

        test.post_process_pv()
        test.post_process_contraction_twist()

        # test.post_process_multi_pv(step_div_factor = 1, extra_str = '_post_buckling')
        # test.post_process_multi_contraction_twist(step_div_factor = 1, extra_str = '_post_buckling')

        # test.post_process_multi_pv(step_div_factor = 1, extra_str = '_post_buckling')
        # test.post_process_multi_contraction_twist(step_div_factor = 1, extra_str = '_post_buckling')

        delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
        delete_extra_files(jname_nonlin)

def run_finger(E_all, w_all, version):

    for i, E in enumerate(E_all):
        geo_props_use = geoPropsFull(10, 18, w_all, 0.856, 0.224, E, E, E, 3*np.pi/2)

        proj_name = 'finger-bender' + '-' + str(version + i)

        test = full_shell(project = proj_name, fullProps = geo_props_use, imperfection = 0.002)
        test.stabilization_factor = 2e-8


        jname_lin = test.run_linear_model()
        num_folds_lin = test.post_process_num_folds()
        printAB(num_folds_lin)
        jname_nonlin = test.make_nonlin_model(temp_mult = 0.6)
        # jname_nonlin = test.make_nonlin_model(temp_mult = final_temp_mult, is_buckling = True, eig_idx = 3)
        # jname_nonlin = test.make_static_dyn_model(temp_mult_static = 0.31, temp_mult_final = 0.7)
        run_inp(jname_nonlin)

        test.post_process_pv()
        test.post_process_bend_angle()

        # test.post_process_multi_pv(step_div_factor = 1, extra_str = '_post_buckling')
        # test.post_process_multi_contraction_twist(step_div_factor = 1, extra_str = '_post_buckling')

        delete_extra_files(jname_lin, ['.fil', '.sta', '.log'])
        delete_extra_files(jname_nonlin)


def main():
    version = 130

    E_all = [0.87, 1.2, 2.16]

    w_all = 5

    # run_finger(E_all, w_all, version)
    run_wrist(E_all, w_all, version, static_end = [0.37, 0.37, 0.37])

main()

