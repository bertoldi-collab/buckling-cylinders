from collections import namedtuple
import numpy as np

geoProps = namedtuple('geoProps',['R','H','w','t','E'])
geoPropsFull = namedtuple('geoPropsFull',['R', 'H', 'w', 't1', 't2','E1', 'E2', 'E_cap', 'theta'])

geo_prop_four = geoProps(10, 18, 5, 0.54, 1.2)
geo_prop_two = geoProps(8.8, 44.75, 5, 1.2, 1.2)
geo_prop_three = geoProps(10, 27, 5, 0.54, 1.2)
geo_prop_bend = geoPropsFull(10, 18, 5, 1.0, 0.3, 1.2, 1.2, 1.2, 3*np.pi/2)

#old props from pre fixing RP days
# geo_prop_four = geoProps(10, 18, 5, 0.56, 1.4)
# geo_prop_two = geoProps(8.8, 44.75, 5, 1.2, 1.4)
# geo_prop_three = geoProps(10, 27, 5, 0.58, 1.4)
# geo_prop_bend = geoPropsFull(10, 18, 5, 0.75, 0.25, 1.4, 1.4, 1.4, 3*np.pi/2)


geo_prop_four_old = geoProps(10, 18, 5, 0.54, 1.0)
geo_prop_two_old = geoProps(8.8, 44.75, 5, 1.15, 1.0)
geo_prop_three_old_geo = geoProps(8, 20, 5, 0.59, 1.4) #fit parameters but on old (R,H)
geo_prop_three_old = geoProps(8, 20, 5, 0.5, 1.0) #non fit parameters on old (R,H)

beamProps = namedtuple('beamProps',['base', 'height'])
beam_prop_v112 = beamProps(5.79404699008849, 3.7670449798818)
beam_prop_v112_eighth = beamProps(3.971594980027066, 1.8781385071582015)
beam_prop_square = beamProps(4.0,4.0)