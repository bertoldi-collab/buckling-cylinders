import sys

sys.path.append('../')

from buck_cylinder_obj import *
from geo_prop import *

idx_try = 12
bdamp = 0.0001
# jname= 'beam-test-v'+str(idx_try)
test = beam_model(project = 'beam-model-v'+str(idx_try), imperfection = 0.05,simpProps = geo_prop_four, beamProps = beam_prop_square)

jname_lin = test.make_linear_model(extra_springs = True)
# run_inp(jname_lin)

# jname_nonlin = test.make_nonlin_model(bdamp = bdamp, eig_idx = 3)
# run_inp(jname_nonlin)


#ALL PRE V5 ARE WRONG BC I DID NOT INPUT REAL MATERIAL PROPERTIES
#v2 is wrong even tho it is beautiful sad
#v3 works
#v4 is v3 but I made it so it should output the fil file
#---------
#oh no; I uhh never put in real material properties oh no
#v5 -10 is too large of a force
#v6 post buckling works w/ force -7
#v7: -8 too large
#v8: -7.5
#v9: -7.307
#v10: trying different beam geometry from analysis of 1/8 of model
#v11: try square for fun