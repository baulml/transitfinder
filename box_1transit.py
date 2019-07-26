import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.timeseries import BoxLeastSquares
from astropy.timeseries import TimeSeries
np.set_printoptions(threshold=sys.maxsize)

filename = "/home/baku/mastDownload/Kepler/kplr011446443_sc_Q113313330333033302/kplr011446443-2011116030358_slc.fits"
filename_tess = "/home/baku/mastDownload/Kepler/kplr011446443_sc_Q113313330333033302/tess2019112060037-s0011-0000000287587961-0143-s_lc.fits"

ts = TimeSeries.read(filename, format="tess.fits")
periodogram = BoxLeastSquares.from_timeseries(ts,'sap_flux')
results = periodogram.autopower(0.1 * u.day)
#print(results.period)

best = np.argmax(results.power)
period = results.period[best]
transit_time = results.transit_time[best]

ts_folded = ts.fold(period=period, midpoint_epoch=transit_time)
plt.plot(ts_folded.time.jd, ts_folded['sap_flux'], 'k.', markersize=1)
plt.xlabel('Time (days)')
plt.ylabel('SAP FLUX (e-/s)')
plt.show()