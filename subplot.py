import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.timeseries import BoxLeastSquares
from astropy.stats import sigma_clipped_stats
from astropy.timeseries import TimeSeries
from astropy.timeseries import aggregate_downsample
from astropy.io import registry, fits
import tkinter as tk
from tkinter import filedialog
np.set_printoptions(threshold=sys.maxsize)


root = tk.Tk()
root.withdraw()
filepath = filedialog.askopenfilename()

keplerID = "555"



ts = TimeSeries.read(filepath,format='kepler.fits')
fig, axs = plt.subplots(2, 2)
fig.suptitle('Kepler Id:011446443', fontsize=15)
fig.suptitle('Kepler Id:' +keplerID, fontsize=15)
fig.set_size_inches(18.5, 10.5)

#full data plot
axs[0, 0].plot(ts.time.jd, ts['sap_flux'], 'k.', markersize=1, label='whole')
axs[0, 0].set_xlabel('Julian Date')
axs[0, 0].set_ylabel('SAP Flux (e-/s)')
axs[0, 0].set_title('full data')


#folded time series
periodogram = BoxLeastSquares.from_timeseries(ts, 'sap_flux')
results = periodogram.autopower(0.1 * u.day)
best = np.argmax(results.power)
period = results.period[best]
transit_time = results.transit_time[best]
ts_folded = ts.fold(period=period, midpoint_epoch=transit_time)
axs[0, 1].plot(ts_folded.time.jd, ts_folded['sap_flux'], 'k.', markersize=1)
axs[0, 1].set_xlabel('Time (days)')
axs[0, 1].set_ylabel('SAP FLUX (e-/s)')
axs[0, 1].set_title('folded time series')

#binned time series
mean, median, stddev = sigma_clipped_stats(ts_folded['sap_flux'])
ts_folded['sap_flux_norm'] = ts_folded['sap_flux'] / median
ts_binned = aggregate_downsample(ts_folded,time_bin_size=0.003 * u.day)
axs[1, 0].plot(ts_folded.time.jd, ts_folded['sap_flux_norm'], 'k.', markersize=1)
axs[1, 0].plot(ts_binned.time_bin_start.jd, ts_binned['sap_flux_norm'], 'r-', drawstyle='default')
axs[1, 0].set_xlabel('Time (days)')
axs[1, 0].set_ylabel('Normalized Flux')
axs[1, 0].set_title('binned time series')

#final transit model
axs[1, 1].plot(ts_folded.time.jd, ts_folded['sap_flux_norm'], 'k.', markersize=1)
axs[1, 1].plot(ts_binned.time_bin_start.jd, ts_binned['sap_flux_norm'], 'r-', drawstyle='default')
axs[1, 1].set_xlabel('Time (days)')
axs[1, 1].set_ylabel('Normalized Flux')
axs[1, 1].set_title('transit model')
plt.show()
