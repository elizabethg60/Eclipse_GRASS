import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import timeit
import os
from pathlib import Path
import sys
import time
import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
from math import floor
from random import choice
import glob
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants
from scipy import interpolate
from scipy import optimize
from scipy import stats
from scipy import ndimage
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

def blaze_correct(lamp_spec,test_ord,wave,spec):
    percentile_flux = np.nanpercentile(spec[test_ord,:], 95)
    # 1. Fit with the lamp
    if test_ord == 15: # order 17 if we don't remove the NaN orders
        # for lower
        ratio = spec[test_ord-1,3000:8000]/lamp_spec[test_ord-1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-1,:]*scale_fac
    
        # for upper
        ratio = spec[test_ord+2,3000:8000]/lamp_spec[test_ord+2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+2,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        fitspec = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/fitspec
        wvl_test = wave[test_ord,:]
    elif test_ord == 16: # order 18 if we don't remove the NaN order
        # for lower
        ratio = spec[test_ord-2,3000:8000]/lamp_spec[test_ord-2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-2,:]*scale_fac
    
        # for upper
        ratio = spec[test_ord+1,3000:8000]/lamp_spec[test_ord+1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+1,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        fitspec = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/fitspec
        wvl_test = wave[test_ord,:]
    # elif test_ord > 70:
    #     lamp_norm = lamp_spec[test_ord,1000:8000]/np.nanmax(lamp_spec[test_ord,1000:8000])
    #     lamp_norm = lamp_norm * np.nanpercentile(spec[test_ord,1000:8000],95)
    #     spec_norm = spec[test_ord,1000:8000]
    #     wvl_test = wave[test_ord,1000:8000]
    #     spec_noblaze = spec_norm/lamp_norm
    else:
        lamp_norm = lamp_spec[test_ord,:]/np.nanmax(lamp_spec[test_ord,:])
        lamp_norm = lamp_norm * np.nanpercentile(spec[test_ord,:],95)
        spec_norm = spec[test_ord,:]
        wvl_test = wave[test_ord,:]
        spec_noblaze = spec_norm/lamp_norm

    # ok try to do a various polynomial fit
    if test_ord > 10 and test_ord < 115:
        wave_portions = np.array_split(wvl_test, 9)
        flux_portions = np.array_split(spec_noblaze, 9)
        maxs_flux = []
        maxs_wvl = []
        for i in range(len(flux_portions)):
            portion = flux_portions[i]
            # maxes
            maxs_flux.append(np.nanmax(portion))
            if len(wave_portions[i][np.where(portion == np.nanmax(portion))]) > 1:
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))][0]))
            else:
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))]))
        # add first -- max in the first 
        maxs_flux.append(np.nanmax(spec_noblaze[0:25]))
        maxs_wvl.append(wave_portions[0][0])

        # add last
        maxs_flux.append(np.nanmax(spec_noblaze[-25:-1]))
        maxs_wvl.append(wave_portions[-1][-1])

        # spline
        spline = interpolate.interp1d(maxs_wvl, maxs_flux)
        spline_lines = spline(wvl_test)

        spec_noblaze = spec_noblaze/spline_lines

    return spec_noblaze*percentile_flux

def reverse_blaze(lamp_spec,test_ord,wave,spec):
    percentile_flux = np.nanpercentile(spec[test_ord,:], 95)
    # 1. Fit with the lamp
    if test_ord == 15: # order 17 if we don't remove the NaN orders
        # for lower
        ratio = spec[test_ord-1,3000:8000]/lamp_spec[test_ord-1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-1,:]*scale_fac
    
        # for upper
        ratio = spec[test_ord+2,3000:8000]/lamp_spec[test_ord+2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+2,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        lamp_norm = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/lamp_norm
        wvl_test = wave[test_ord,:]
    elif test_ord == 16: # order 18 if we don't remove the NaN order
        # for lower
        ratio = spec[test_ord-2,3000:8000]/lamp_spec[test_ord-2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-2,:]*scale_fac
    
        # for upper
        ratio = spec[test_ord+1,3000:8000]/lamp_spec[test_ord+1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+1,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        lamp_norm = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/lamp_norm
        wvl_test = wave[test_ord,:]
    else:
        lamp_norm = lamp_spec[test_ord,:]/np.nanmax(lamp_spec[test_ord,:])
        lamp_norm = lamp_norm * np.nanpercentile(spec[test_ord,:],95)
        spec_norm = spec[test_ord,:]
        wvl_test = wave[test_ord,:]
        spec_noblaze = spec_norm/lamp_norm

    # ok try to do a various polynomial fit
    if test_ord > 10 and test_ord < 115:
        wave_portions = np.array_split(wvl_test, 9)
        flux_portions = np.array_split(spec_noblaze, 9)
        maxs_flux = []
        maxs_wvl = []
        for i in range(len(flux_portions)):
            portion = flux_portions[i]
            # maxes
            maxs_flux.append(np.nanmax(portion))
            if len(wave_portions[i][np.where(portion == np.nanmax(portion))]) > 1:
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))][0]))
            else:
                #print(wave_portions[i][np.where(portion == np.nanmax(portion))])
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))]))
        # add first -- max in the first 
        maxs_flux.append(np.nanmax(spec_noblaze))
        maxs_wvl.append(wave_portions[0][0])

        # add last
        maxs_flux.append(np.nanmax(spec_noblaze))
        maxs_wvl.append(wave_portions[-1][-1])

        # spline
        spline = interpolate.interp1d(maxs_wvl, maxs_flux)
        spline_lines = spline(wvl_test)
        #spline_lines = np.pad(spline_lines, (1000,9216-8000), 'empty')
        #print(spline_lines)

        spec_noblaze = spec_noblaze/spline_lines
        return lamp_norm*spline_lines
    return lamp_norm
    #return wvl_test, spec_noblaze*percentile_flux

def wave_vel_shift(wave, vel):
    '''
    Appropriately shifts wavelength(s) 'wave' for velocity 'vel'
    
    Parameters
    ----------
    wave : float or array of floats
        Wavelength(s) to be stretched
    
    vel : float 
        Velocity to shift wave [km/s]
    Returns
    -------
    wave_shifted : float or array of floats
        Wavelength(s) shifted by velocity 'vel'
    '''

    # shift wavelengths to provided velocity
    wave_shifted = wave * (1. + vel / constants.c)

    return wave_shifted

def airtovac(wvl_air):
    s = 104 / wvl_air
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    wvl_vac = wvl_air * n
    
    return wvl_vac

def primed(wave,flux):
    return np.gradient(flux,wave)

def func(X, A, b):
    s_temp,s_prime = X
    return A*s_temp + b*s_prime

def get_coefficients(object_flux, template_wave, template_flux, index=None, window_min=None, window_max=None):
    # object_window = object[index][window_min:window_max]
    # template_window = template_flux[index][window_min:window_max]
    # template_wave = template_wave[index][window_min:window_max]
    # template_prime = primed(template_wave,template_window)
    #----------------------------Passing In Pre-Cut Windows---------------------------------
    object_window = object_flux
    #template_wave *= 10**(-9)
    template_window = template_flux
    template_prime = primed(template_wave, template_window)
    # initial guesses for A,v (will guess 1 otherwise)
    
    while True:
        try:
            popt, pcov = curve_fit(func, xdata=(template_window,template_prime), ydata=object_window)
            #return leastsq(func_lstsq,0)
            return popt, pcov
        except (RuntimeError, TypeError, NameError, ValueError):
            popt, pcov = [],[]
            return popt, pcov

# make sure the units of the wavelength and the b coefficient (look at gradient/wave) are the same
def rv_shift(A, b, wvl):
    #wvl*=10**(-9)
    return ((constants.c)/wvl)*(b/A)

# Adeltav as one v. A and delta v?
def rv_error(A, b, sigmaA,sigmab, wvl):
    #wvl*=10**(-9)
    return ((constants.c)/wvl) * np.sqrt((1/A)**2*sigmab**2 + (-b/A**2)**2*sigmaA**2)

# keep a tally of what is being tossed
def sigma_clip(data, carried_data, a):
    #---------------------------Calculate Median and Sigma----------------------------------
    med = np.median(data) # want this to just take the lower of the two
    sigma = np.std(data, axis=0)
    #-----------------------------Calculate Interval--------------------------------
    lower_lim = med - a*sigma
    upper_lim = med + a*sigma
    where = np.where((data>lower_lim) & (data<upper_lim))[0]
    #-----------------------------Cut--------------------------------
    return data[data.index.isin(where)], carried_data[carried_data.index.isin(where)]
    #return data[where], carried_data[where]

def rv_cleaning(rvs, errors):
    #-----------------------------0. Convert to arrays?--------------------------------
    # errors = np.array(errors)
    # rvs = np.array(rvs)
    #-----------------------------1. 2 6-sigma cuts on the ERROR--------------------------------
    errors, rvs = sigma_clip(errors, rvs, 6)
    #print(len(errors), len(rvs))
    errors, rvs = sigma_clip(errors, rvs, 6)
    #-----------------------------2. 2 6-sigma cuts on the RVS--------------------------------
    rvs, errors = sigma_clip(rvs, errors, 6)
    rvs, errors = sigma_clip(rvs, errors, 6)
    return rvs, errors

def total_rv(rvs, errors):
    numerator = 0.0
    denominator = 0.0
    j=0
    for i in range(len(rvs)):
        error = errors[i]
        rv = rvs[i]
        if not np.isnan(error):
            if not np.isnan(rv):
                j+=1
                numerator += (rv * (1/(error**2)))
                denominator += (1/(error**2))
    return numerator/denominator

def total_err(line_errs):
    sum_square = 0
    for err in line_errs:
        sum_square += (1/(float(err)**2))
    return (1/(sum_square)**(1/2))

def rms(err_col):
    ms_total = err_col - np.median(err_col)
    return np.sqrt(np.mean(ms_total**2))

def calculate_rvs(jd, fits_file, template_wave, template_flux, wavelengths, plot=False, blazecorr=True, quality_cut=None):
    #----------------------------Create Empty Template to Edit---------------------------------
    temp_flux = np.zeros(np.shape(template_flux))
    temp_wave = np.zeros(np.shape(template_wave))
    
    for ind in range(len(template_flux)):
        temp_flux[ind] = template_flux[ind]
        temp_wave[ind] = template_wave[ind]
    #----------------------------Read In Observation---------------------------------
    # read in wave, spec, and var
    # the actual spectrum (in flux units)
    flux_obs = fits.getdata(fits_file, 'SCIFLUX')

    # the variance or noise of the spectrum
    var_obs = fits.getdata(fits_file, 'SCIVAR')

    # the wavelength solution of the spectrum
    wave_obs = fits.getdata(fits_file, 'SCIWAVE') / 10. # nm

    blaze_obs = fits.getdata(fits_file, 'SCIBLAZE')
    # while True:
    #     try:
    #         blaze_obs = fits.getdata(fits_file, 'SCIBLAZE')
    #     except KeyError:
    #         return pd.DataFrame(data=np.zeros((len(wavelengths),5)), columns=['wavelength', 'rv', 'error', 'order', 'date'])
    #---------------------------Clean NaN orders----------------------------------
    # get rid of NaN orders, if any
    # get order mean wavelengths
    nord_obs = flux_obs.shape[0]
    wvls_ord_obs = []

    # get mean order wavelengths
    for ind, order in enumerate(range(nord_obs)):
        wvls_ord_obs.append(np.nanmean(wave_obs[ind,:]))
    wvls_ord_obs = np.asarray(wvls_ord_obs)

    # now cut the nan orders
    nan_ords_obs = np.isnan(wvls_ord_obs)

    #--------------------------------Barycentric Velocity of Observation File-------------------------------
    # red-shift each spectrum to account for the barrycentric velocity
    # shift: lambda_rest = lambda_obs / (1+v/c)
    # pull velocity information from header -- this is from neid_measure_act slash this is how we were shifting for velocity offset there
    hdr_spec = fits.getheader(fits_file)
    ccf_header = fits.getheader(fits_file,'CCFS')


    bc_vel = hdr_spec['SSBRV160'] # km/s # 'Barycentric corr. (km/s) for order 160' why do you take the 160th order? # barrycentric correction for each order WE CAN USE ONE FOR ALL ORDERS BC SUN
    q_vel = hdr_spec['QRV'] # km/s # 'Queue RV for Star'
    #measured_rv = ccf_header['CCFRVMOD'] # 'Bary. corrected RV for weighted orders' -- this is the actual RV

    # shift wavelengths to stellar rest frame
    vel_factor = -1. * (float(bc_vel)) + float(q_vel)
    z_obs = (vel_factor*10**3)/constants.c # putting it into m/s because that is what I am working in
    z_brv_file = hdr_spec['SSBZ060']
    
    #--------------------------------Shift Template------------------------------- # do I want to add or subtract it?
    wave_obs = wave_obs*(1-z_obs) # or is it divded by?
    #wave_obs = wave_obs*(1+z_brv_file)

    #--------------------------------Adjust Total Flux per Order-------------------------------
    #print(flux_obs[65][3000:3016])
    # for ind in range(len(wave_obs)):
    #     if 'CCFWT'+str(ind).rjust(3,'0') in ccf_header:
    #         if ccf_header['CCFWT'+str(ind).rjust(3,'0')] is not None:
    #             if ccf_header['CCFWT'+str(ind).rjust(3,'0')] != 0:
    #                 flux_obs[ind] *= ccf_header['CCFWT'+str(ind).rjust(3,'0')]
    #print(flux_obs[65][3000:3016])
    #------------------------------Interpolate to Match Wavelengths-------------------------------
    # Note: only looking at indices 1000:8000 because worried abt the nan buffers at the begining and end of the 
    for ind in range(len(template_wave)):
        f = interpolate.interp1d(temp_wave[ind][1000:8000], temp_flux[ind][1000:8000], fill_value="extrapolate")# can only use the fill value if your new values are BARELY out of range
        temp_flux[ind][1000:8000] = f(wave_obs[ind][1000:8000])
        temp_wave[ind][1000:8000] = wave_obs[ind][1000:8000]
    #------------------------------Mask-------------------------------
    # if mask[0][0]==1:
    #     mask = 1-mask
    # flux_obs *= mask
    
    #------------------------------Blaze Correct: Template-------------------------------
    plt.rcParams['figure.figsize'] = (12,5)
    # blaze correct each order
    if blazecorr==True:
        for i in range(len(temp_flux)):
            if not nan_ords_obs[i]:
                #fit_ord = reverse_blaze(i,wave_obs,flux_obs) # level 2 data use below method
                # temp_flux_ord = temp_flux[i]*fit_ord
                temp_flux_ord = temp_flux[i] * blaze_obs[i]
                # temp_flux_ord = temp_flux_ord * ndimage.maximum_filter1d(wave_obs[i]/blaze_obs[i], 3, axis=- 1, output=None, mode='reflect', cval=0.0, origin=0)
                temp_flux[i] = temp_flux_ord/np.nanmax(temp_flux_ord)
            else:
                temp_flux[i] = temp_flux[i]/np.nanmax(temp_flux[i])
    # elif blazecorr=='opposite':
    #     for i in range(len(flux_obs)):
    #         if not nan_ords[i]:
    #             flux_obs[i] = blaze_correct(i, wave_obs, flux_obs)
    #------------------------------Plot to Check-------------------------------
    if plot==True:
        plt.rcParams['figure.figsize'] = (12,5)
        plt.plot(template_wave[18], template_flux[18])
        plt.title('Template Spectrum')
        plt.show()
        plt.plot(temp_wave[18],temp_flux[18])
        plt.title('Shifted and Blazed Template Spectrum')
        plt.show()
        plt.plot(wave_obs[18],flux_obs[18])
        plt.title('Sample Spectrum')
        plt.show()
    #------------------------------Initialize Data Frame-------------------------------
    df = pd.DataFrame(data=wavelengths, columns=['wavelength', 'rv', 'error', 'order', 'date'])
    #print(df['wavelength'])
    #------------------------------Loop Through Wavelengths-------------------------------
    for i in range(len(df)):
        df['date'][i] = jd
        wavelength = df['wavelength'][i]/10 # made vacuum before
        #------------------------------Retrive Index-------------------------------
        ind_line = np.nanargmin(np.abs(wvls_ord_obs-wavelength))
        df['order'][i] = ind_line
        #------------------------------Flux normalize order-------------------------------
        if 'CCFWT'+str(173-ind_line).rjust(3,'0') in ccf_header:
            if ccf_header['CCFWT'+str(173-ind_line).rjust(3,'0')] is not None:
                if ccf_header['CCFWT'+str(173-ind_line).rjust(3,'0')] != 0:
                    flux_obs[ind_line] *= ccf_header['CCFWT'+str(173-ind_line).rjust(3,'0')]
        #------------------------------Retrive Windows-------------------------------
        ind_where_obs = np.where(wave_obs[ind_line] == wave_obs[ind_line][(np.abs(wave_obs[ind_line] - wavelength)).argmin()])[0][0]
        ind_where_temp = np.where(temp_wave[ind_line] == wave_obs[ind_line][(np.abs(temp_wave[ind_line] - wavelength)).argmin()])[0][0]

        if ind_line < 55:
            ind_where_obs += 3
            ind_where_temp += 3
        
        #looking at wavelength difference -- none because shifting! it fixed it!
        wave_obs_window = wave_obs[ind_line,ind_where_obs-8:ind_where_obs+8]
        temp_wave_window = temp_wave[ind_line, ind_where_temp-8:ind_where_temp+8]

        flux_obs_window = flux_obs[ind_line,ind_where_obs-8:ind_where_obs+8]
        temp_flux_window = temp_flux[ind_line, ind_where_temp-8:ind_where_temp+8]

        # if wavelength > 460:
        #     plt.figure()
        #     plt.scatter(wave_obs[ind_line,ind_where_obs-20:ind_where_obs+20], flux_obs[ind_line,ind_where_obs-20:ind_where_obs+20]/np.nanmedian(flux_obs[ind_line,ind_where_obs-20:ind_where_obs+20]), label="observation")
        #     plt.scatter(temp_wave[ind_line, ind_where_temp-20:ind_where_temp+20], temp_flux[ind_line, ind_where_temp-20:ind_where_temp+20]/np.nanmedian(temp_flux[ind_line, ind_where_temp-20:ind_where_temp+20]), label="template")
        #     plt.legend()
        #     plt.show()

        #------------------------------Calculate Coefficients: Variance (non-blaze corrected)-------------------------------
        coeffs, var = get_coefficients(flux_obs_window, temp_wave_window, temp_flux_window)
        
        if len(coeffs) > 1:
            A = coeffs[0]
            deltav = coeffs[1]
            perr = np.sqrt(np.diag(var))
            sigmaA = perr[0]
            sigmadeltav = perr[1]
            #print(perr)
            rvi = rv_shift(A, deltav, wavelength)
            errori = rv_error(A,deltav, sigmaA, sigmadeltav, wavelength)
            # if -1 < rvi < 1:
            #     df['rv'][i] = np.nan
            #     df['error'][i] = np.nan
            # else:
            df['rv'][i] = rv_shift(A, deltav, wavelength)
            df['error'][i] = rv_error(A,deltav, sigmaA, sigmadeltav, wavelength)
            #return rv_shift(A, deltav, wavelength), rv_error(A,deltav, sigmaA, sigmadeltav, wavelength)
        else:
            df['rv'][i] = np.nan
            df['error'][i] = np.nan
            #return np.nan, np.nan

    return df

#-------------------------------Data Files------------------------------
# data files
local_path = os.getcwd()
data_path = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/15/"
spec_fits_files = glob.glob(data_path + '/*.fits', recursive=True)

#------------------------------Wavelengths File-------------------------------
wavelengths = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
wavelengths_vac = []
for i in range(len(wavelengths)):
    wavelengths_vac.append(airtovac(wavelengths[i]))
wavelengths_vacDF = pd.DataFrame(wavelengths_vac, columns=['wavelength'])
wavelengths_vac = wavelengths_vacDF['wavelength']

#------------------------------Read-in Masterwave and Template DF-------------------------------
templatePath = os.path.join(local_path, 'test_temp.csv')
template = np.loadtxt(open(templatePath, "rb"), delimiter=",")

masterwavePath = os.path.join(local_path, 'test_wave.csv')
masterwave = np.loadtxt(open(masterwavePath, "rb"), delimiter=",")

#------------------------------Initialize Saving Path-------------------------------
rvs_path = os.path.join(local_path, 'test_rvs.csv')
sindex_path = os.path.join(local_path, 'test_sindex.csv')
#------------------------------Loop-------------------------------
cahk = []
i = 0

for fits_file in spec_fits_files[0:-146]:
    i += 1
    # read in headers
    main_header = fits.getheader(fits_file)
    ccf_header= fits.getheader(fits_file,'CCFS')
    # get jd and bc vel
    jd = fits.getheader(fits_file)['OBSJD']
    bc_vel_file = main_header['SSBRV160'] # km/s # 'Barycentric corr. (km/s) for order 160' why do you take the 160th order? # barrycentric correction for each order WE CAN USE ONE FOR ALL ORDERS BC SUN
    # read in cahk
    activity = fits.getdata(fits_file, 'ACTIVITY')
    cahk.append((jd, activity[0][1], activity[0][2]))

    # calculate rvs
    temp_refline_rvs = calculate_rvs(jd, fits_file, masterwave, template, wavelengths_vac, plot=False)
    temp_refline_rvs['wavelength'] = wavelengths # try to add the og wavelengths instead of the converted ones
    temp_refline_rvs = temp_refline_rvs.dropna()

    # calculate integrated RV
    my_rv = total_rv(np.array(temp_refline_rvs['rv']), np.array(temp_refline_rvs['error'])) * 10**(-3)
    my_rv_err = total_err(np.array(temp_refline_rvs['error'])) * 10**(-3)

    # add dfs
    total_df = pd.DataFrame(data=[['total',my_rv,my_rv_err,'total',jd]], columns=['wavelength', 'rv', 'error', 'order', 'date'])
    neid_df = pd.DataFrame(data=[['neid total',fits.getheader(fits_file,'CCFS')['CCFRVMOD'], fits.getheader(fits_file,'CCFS')['DVRMSMOD'], 'total', jd]], columns=['wavelength', 'rv', 'error', 'order', 'date'])
    temp_refline_rvs = temp_refline_rvs._append(total_df)
    temp_refline_rvs = temp_refline_rvs._append(neid_df)

    # save
    if i == 1:
        temp_refline_rvs.to_csv(rvs_path, mode='a', index=False, header=True)
    else:
        temp_refline_rvs.to_csv(rvs_path, mode='a', index=False, header=False) 