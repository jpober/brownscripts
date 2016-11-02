from astropy.time import Time
import numpy as np
import os.path as op
from astropy.utils import iers
import uvdata

data_path = op.join(uvdata.__path__[0], 'data')
iers_a = iers.IERS_A.open(op.join(data_path, 'finals.all'))

def set_lsts_from_time_array(uvd):
      lsts = []
      curtime = uvd.time_array[0]
      for ind, jd in enumerate(uvd.time_array):
          if ind == 0 or not np.isclose(jd, curtime, atol=1e-6, rtol=1e-12):
#              print 'Curtime/jd: ', curtime, jd
              curtime = jd
              latitude, longitude, altitude = uvd.telescope_location_lat_lon_alt_degrees
#              print 'Loc: ', latitude, longitude, altitude
              t = Time(jd, format='jd', location=(longitude, latitude))
#              print 't: ', t
              t.delta_ut1_utc = iers_a.ut1_utc(t)
#              print 't.delta_ut1_utc: ', t.delta_ut1_utc
          print "LST: ", t.sidereal_time('apparent')
	  lsts.append(t)
      return lsts
