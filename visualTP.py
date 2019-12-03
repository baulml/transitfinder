import sys
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

def main():
    root = tk.Tk()
    root.withdraw()
    filepath = filedialog.askopenfilename()
    hdulist = fits.open(filepath)
    telescope = hdulist[0].header['telescop'].lower()
    if telescope == 'tess':
        hdu = hdulist['LIGHTCURVE']
        print("TESS file recognized")
        tessID = str(hdulist[1].header['TICID'])
        print(tessID)
        ts = TimeSeries.read(filepath, format='tess.fits')
    elif telescope == 'kepler':
        hdu = hdulist[1]
        print("Kepler file recognized")
        keplerID = str(hdulist[0].header['KEPLERID'])
        print(keplerID)
        ts = TimeSeries.read(filepath, format='kepler.fits')
    else:
        raise NotImplementedError("{} is not implemented, only KEPLER or TESS are "
                                  "supported through this reader".format(hdulist[0].header['telescop']))
    fig, axs = plt.subplots(2, 2)
    if telescope == 'tess':
        fig.suptitle('TESS Id:' +tessID, fontsize=15)
    if telescope == 'kepler':
        fig.suptitle('Kepler Id:' + keplerID, fontsize=15)
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
