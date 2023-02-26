import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 134
bdamp = 0.0001
# jname_nonlin = '4fold-imperfection-'+str(idx_try) + '_post_buckling'
test = full_shell(project = '2fold-imperfection-'+str(idx_try), simpProps = geo_prop_two, imperfection = 0.0)
test.run_linear_model()

jname_multi = test.make_nonlin_multi_buckle(bdamp, max_temp_mult = 0.0194, num_steps = 10)
run_inp(jname_multi)
test.post_process_multi_buckle()
delete_extra_files(jname_multi)

#v100: freq steps failed after like 8
#v101: changed max steps to 10 (11 in code bc 1st is deleted)
#v102: changed max_temp_mult from 0.4 --> 0.3; they are all running up to -0.081897
#v103: changed max_temp_mult from 0.3 --> 0.245 should converge? --> works!!
#v104: just trying to get more fine-grained info: changing max steps from 10 (11 in code) --> 15 (16 in code)
#v105: try 1 step see if it works w/ temp_mult = 0.27: did work!
#v106: do temp_mult = 0.27 for 10 steps: did not work!; then I tried for 2 steps and overwrote oops
#v107: testing 2 steps w/ ramp for static instead of smooth step: failed (at normal time I think)
#v108: trying 1 step w/ ramp: worked!! aaahhh
#v109: [BAD] 1 step but with nu = 0 set in section, and making initial inc for static 0.01 instead of 0.001: nu should not be 0
#v110: [BAD] ?????????
#v111: fixed v109 and made nu = 0.5 again: running again to check current state w/ max_temp_mult = 0.28
#v112: same as v111 but w/ num_steps = 1
#v113: same as v112 but w/ bd cond being disp instead of velocity (idk why it was velocity!!)
#v114: trying 0.27 as max_temp_mult again [works?]
#v115: trying 0.275 [works!]
#v116: trying 0.2775 [did not work]
#v117: trying 0.275 w/ 10 steps
#v118: trying 0.275 with david's 25U change
#v119: mult initial and max increment by len(temp_list) to try to scale like 1 step sim
#v120: sanity checking the 1 step thing [same as v115 basically]
#v121: [BAD] changed D_1 values to 0.0833 for all materials per david's suggestion: trying 1 step still
#v122: [BAD] same as v121 but w/ 10 steps
#v123: [BAD] "Maybe try with factor of 10 and 0.1" [factor of 10]
#v124: [BAD] "Maybe try with factor of 10 and 0.1" [factor of .1]
#v125: no imperfection, resolve initial buckle (usually 0.05, max temp has been 0.275)- did not finish
#v126: trying 0.025
#v127: trying 0.02
#v128: trying 0.015: finished but did not resolve freq = 0
#v129: trying 0.017
#v130: trying 0.016

#v131: doing 2 folds and same 0.016 (still called 4 fold oops)
#v132: 0.017
#v133: 0.02- stopped at 9.71 --> go until 0.0194
#v134: 0.0194
