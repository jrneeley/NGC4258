import numpy as np
import matplotlib.pyplot as mp
import glob
import re
import sys

star_list = glob.glob('lc_data/*.dat')
star_names = np.zeros(len(star_list), dtype=int)
for ii, star in enumerate(star_list):
    star_names[ii] = int(star[8:-4])
star_names = np.sort(star_names)

f = open('results.dat', 'w')
for star in star_names:

    print star
    star_file = 'lc_data/'+str(star)+'.dat'
    dtype = np.dtype([('name', int), ('filter', 'S5'), ('mjd', float),
        ('mag', float), ('err', float)])
    data = np.loadtxt(star_file, dtype=dtype, skiprows=2)

    Jmag = data['mag'][data['filter'] == 'F125W']
    Hmag = data['mag'][data['filter'] == 'F160W']
    Jerr = data['err'][data['filter'] == 'F125W']
    Herr = data['err'][data['filter'] == 'F160W']
    Jmjd = data['mjd'][data['filter'] == 'F125W']
    Hmjd = data['mjd'][data['filter'] == 'F160W']

    # clean up bad data
    J_avg_err = np.median(Jerr)
    J_std_err = np.std(Jerr)
    J_bad = (np.abs(Jerr - J_avg_err) > 3*J_std_err) #| (Jerr > 0.5)
    H_avg_err = np.median(Herr)
    H_std_err = np.std(Herr)
    H_bad = (np.abs(Herr - H_avg_err) > 3*H_std_err) #| (Herr > 0.5)
    Jmag[J_bad] = np.nan
    Hmag[H_bad] = np.nan
    Jerr[J_bad] = np.nan
    Herr[H_bad] = np.nan

    # set up figure
    fig = mp.figure()
    ax1 = mp.subplot(111)
    ax1.set_xlabel('Phase')
    ax1.set_ylabel('mag')
    ax1.set_ylim((np.nanmax(Jmag+0.2), np.nanmin(Hmag-1.0)))
    ax1.set_xlim((-0.5, 2.0))
    ax1.set_title('NGC 4258: Mira '+ str(star))

    # calculate color
    Jmask = np.in1d(np.floor(Jmjd), np.floor(Hmjd)) == True
    Hmask = np.in1d(np.floor(Hmjd), np.floor(Jmjd)) == True
    color = Jmag[Jmask] - Hmag[Hmask]
    color_err = np.sqrt(Jerr[Jmask]**2 + Herr[Hmask]**2)
    color_epochs = Hmjd[Hmask]


    filters = np.unique(data['filter'])
    colors = ['xkcd:ocean blue', 'xkcd:rose']
    for ind, filt in enumerate(filters):
        if filt == 'F160W':
            mag = Hmag
            err = Herr
            mjd = Hmjd
        if filt == 'F125W':
            mag = Jmag
            err = Jerr
            mjd = Jmjd

        mag_max = np.nanmax(mag)
        mag_min = np.nanmin(mag)
        epoch_max = mjd[mag == mag_max][0]
        epoch_min = mjd[mag == mag_min][0]

        quality_flag = 0
        if (epoch_max == mjd[0]) | (epoch_min == mjd[0]): quality_flag=1
        if (epoch_max == mjd[-1]) | (epoch_min == mjd[-1]): quality_flag=1

        # derive light curve parameters
        amplitude = mag_max - mag_min
        average2 = np.nanmean(mag)
        average = (mag_max + mag_min)/2.0
        average_err = amplitude/(12.*np.sqrt(10.))
        period = 2.*np.abs(epoch_max - epoch_min)

        # derive distance modulus
        if filt == filters[1]:
            dm = average + 3.15*np.log10(period) - 19.06 + 18.50
        if filt == filters[0]:
            dm = average + 2.92*np.log10(period) - 19.37 + 18.50

        print '\n{} results: '.format(filt)
        print 'Period = {:5.1f}'.format(period)
        print 'Amplitude = {}'.format(amplitude)
        print '<{}> = {:5.2f}'.format(filt, average)
        print '<{}> err = {:4.2f}'.format(filt, average_err)
        print 'mu = {:4.2f}\n'.format(dm)

        #phase data
        #phase_all = np.mod((mjd_all - epoch_min)/period, 1)
        phase = np.mod((mjd - epoch_min)/period,1)

        # repeat data for plots
        phase2 = np.concatenate((phase, phase+1))
        mag2 = np.tile(mag, 2)
        err2 = np.tile(err, 2)
    #    ax1.errorbar(phase_all, mag_all, yerr=err_all, fmt='o', mfc='none', color=colors[ind])
        ax1.errorbar(phase2, mag2, yerr=err2, fmt='o', color=colors[ind])

        #cosine function
        t = np.linspace(0, 2, num=100)
        y = amplitude/2.*np.cos(2*np.pi*(t-0.5)) + average
        ax1.plot(t, y, color=colors[ind])
        ax1.plot([0,2], [np.nanmean(y), np.nanmean(y)], 'k--')
        ax1.annotate(filt, xy=(-0.4, average))

        # goodness of fit
        yfit = 0.5*amplitude*np.cos(2*np.pi*(phase-0.5)) + average
        residuals = mag - yfit
        N = float(len(residuals))
        gof = np.nansum(residuals**2)/(N-1.0)

        # save results to file
        data_save = np.array(zip([star], [filters[ind]], [period], [average],
            [average_err], [amplitude], [dm], [quality_flag], [gof]),
            dtype=[('name', int), ('filter', 'S5'), ('period', float),
            ('avg', float), ('err', float), ('amp', float), ('dm', float),
            ('qual', int), ('gof', float)])
        np.savetxt(f, data_save, fmt='%6i %5s %5.1f %5.2f %4.2f %4.2f %5.2f %2s %5.3f')

    # add color curve
    color2 = np.tile(color, 2)
    color_err2 = np.tile(color_err, 2)
    color_offset = np.nanmin(Hmag) - 1.2
    color_phase = np.mod((color_epochs - epoch_min)/period, 1)
    color_phase2 = np.concatenate((color_phase, color_phase+1))
    ax1.errorbar(color_phase2, color2 + color_offset, yerr=color_err2, fmt='o', color='xkcd:dark lavender')
    ax1.plot([0,2], [np.nanmean(color)+color_offset, np.nanmean(color)+color_offset], 'k--')
    color_str = '('+filters[0]+'-'+filters[1]+')'
    ax1.annotate(color_str, xy=(-0.4, np.nanmean(color)+color_offset-0.05))

    if quality_flag == 0:
        mp.savefig('good_plots/'+str(data['name'][0])+'.eps', format='eps')
    if quality_flag == 1:
        mp.savefig('bad_plots/'+str(data['name'][0])+'.eps', format='eps')
    mp.close()
f.close()
