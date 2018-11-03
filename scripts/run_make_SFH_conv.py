# coding: utf-8

#!/usr/bin/env python2

"""
Created on October 14, 2018

@author: sandro.tacchella@cfa.harvard.edu

Comments:
run batch file:
sbatch --array=1-XX submission_script_make_SFH.sh, with XX given by number_of_bins

"""

# import modules

import numpy as np
from DELCgen import *
import os
import h5py
import argparse
import make_SFH_from_PSD
from astropy.table import Table


from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
print cosmo


# define working directory

path_main = os.environ['WDIR_SFH_variability']
path_SFH_cat = path_main + 'catalogs/SFH/'


# read in response function file for convolution

path_response_file = path_main + 'scripts/variability_tools/response_function.dat'
response_function = Table.read(path_response_file, format='ascii')


# read in command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("--idx_key", type=int, help="iteration variable")
parser.add_argument("--filename_SFH", type=str, help="filename of SFH file")
parser.add_argument("--sfh_res", type=float, help="time resolution of SFH in Gyr")
parser.add_argument("--redshift_start", type=float, help="start of SFH in redshift space")
parser.add_argument("--redshift_end", type=float, help="end of SFH in redshift space")
parser.add_argument("--number_galaxies", type=int, help="Number of SFHs/galaxies")
parser.add_argument("--scatter_MS_0", type=float, help="scatter of MS (on shortest timescale)")
parser.add_argument("--list_of_slopes", nargs='+', type=float, help="slope on the short time scales")
parser.add_argument("--list_of_breaks", nargs='+', type=float, help="timescale of break")
parser.add_argument("--aliasTbin", type=float, help="oversampling of time series (resolution)")

args = parser.parse_args()


# get grid of times

time_cosmo = np.arange(cosmo.age(args.redshift_start).value, cosmo.age(args.redshift_end).value, args.sfh_res)
number_steps = len(time_cosmo)
print('number of time steps: ' + str(number_steps))
print('number of galaxies: ' + str(args.number_galaxies))

time_SFH = time_cosmo-time_cosmo[0]
print time_SFH
print len(time_SFH)

### set parameters of light curve model ###
# amplitude => seems not to work
amp = args.scatter_MS_0

# slope on the short time scales
list_of_slopes = np.array(args.list_of_slopes)
print list_of_slopes
list_of_slopes_str = list_of_slopes.astype('str')

#slope on the long time scales
slope_high = 0.0

#frequency of the break
list_of_break_timescales = np.array(args.list_of_breaks)
print list_of_break_timescales
list_of_breaks = 1.0/list_of_break_timescales
list_of_breaks_str = list_of_break_timescales.astype('str')

# intercept/offset: adds constant to the final output
# c = 10.0


# define function for translating key

def translate_key(idx):
    '''
    Translate key to indices.
    '''
    idx_slope = (idx-1) / len(list_of_breaks)
    idx_break = (idx-1) % len(list_of_breaks)
    return(idx_slope, idx_break)


# generate SFH from power-spectrum density
ii_slope, jj_break = translate_key(args.idx_key)
print 'slope = ', list_of_slopes[ii_slope]
print 'break = ', list_of_breaks[jj_break]
array_of_SFH = make_SFH_from_PSD.create_family_SFHs(args.number_galaxies, list_of_slopes[ii_slope], list_of_breaks[jj_break], args.scatter_MS_0, args.aliasTbin, number_steps, amp, slope_high)


# compute SFRs weighted by luminosity evolution

# define high res time bins for convolution

N = 10  # sub-resolution (10 times)
time_hr = np.arange(0.0, np.max(time_SFH)+args.sfh_res/N, args.sfh_res/N)


# make high res grid for response function

response_function_hr_dict = {}

for ii_key in response_function.keys()[1:]:
    response_function_hr_dict[ii_key] = np.interp(time_hr, response_function['time'], response_function[ii_key])
    response_function_hr_dict[ii_key] = response_function_hr_dict[ii_key]/np.sum(response_function_hr_dict[ii_key])


def convolve(SFR_hr, indicator):
    SFR_conv = np.convolve(SFR_hr, response_function_hr_dict[indicator], mode='full')
    SFR_conv = np.convolve(SFR_conv, np.ones((N,))/N, mode='valid')[::N][:len(time_SFH)]
    return(SFR_conv)


# compute array

array_i1500 = []
array_i2800 = []
array_u = []
array_v = []
array_2mass_j = []
array_wise_w3 = []
array_UVIR = []
array_Ha = []


for ii in range(args.number_galaxies):
    print 'progress (%) :  ', np.round(100.0*(ii+1)/args.number_galaxies, 1)
    SFR = 10**array_of_SFH[ii]  # assume MS is 1 Msun/yr
    SFR_hr = np.interp(time_hr, time_SFH, SFR)
    # set SFH parameters
    if (ii == 0):
        array_i1500 = np.log10(convolve(SFR_hr, 'FUV'))
        array_i2800 = np.log10(convolve(SFR_hr, 'NUV'))
        array_u = np.log10(convolve(SFR_hr, 'u'))
        array_v = np.log10(convolve(SFR_hr, 'V'))
        array_2mass_j = np.log10(convolve(SFR_hr, 'J'))
        array_wise_w3 = np.log10(convolve(SFR_hr, 'W3'))
        array_UVIR = np.log10(convolve(SFR_hr, 'UV+IR'))
        array_Ha = np.log10(convolve(SFR_hr, 'Ha'))
    else:
        array_i1500 = np.vstack([array_i1500, np.log10(convolve(SFR_hr, 'FUV'))])
        array_i2800 = np.vstack([array_i2800, np.log10(convolve(SFR_hr, 'NUV'))])
        array_u = np.vstack([array_u, np.log10(convolve(SFR_hr, 'u'))])
        array_v = np.vstack([array_v, np.log10(convolve(SFR_hr, 'V'))])
        array_2mass_j = np.vstack([array_2mass_j, np.log10(convolve(SFR_hr, 'J'))])
        array_wise_w3 = np.vstack([array_wise_w3, np.log10(convolve(SFR_hr, 'W3'))])
        array_UVIR = np.vstack([array_UVIR, np.log10(convolve(SFR_hr, 'UV+IR'))])
        array_Ha = np.vstack([array_Ha, np.log10(convolve(SFR_hr, 'Ha'))])


# save SFH dictionary (contains SFH)
file_name = path_SFH_cat + args.filename_SFH + 'alpha_' + list_of_slopes_str[ii_slope] + '_tau_' + list_of_breaks_str[jj_break] + '.hdf5'

try:
    os.remove(file_name)
except OSError:
    pass


f = h5py.File(file_name, 'w')
# add SFH
grp_SFH = f.create_group("SFH")
grp_SFH.create_dataset('time_cosmo', data=time_cosmo)
grp_SFH.create_dataset('time', data=time_SFH)
grp_SFH.create_dataset('SFR', data=array_of_SFH)
grp_lum = f.create_group("SFH_conv")
grp_lum.create_dataset('SFR_FUV', data=array_i1500)
grp_lum.create_dataset('SFR_NUV', data=array_i2800)
grp_lum.create_dataset('SFR_u', data=array_u)
grp_lum.create_dataset('SFR_v', data=array_v)
grp_lum.create_dataset('SFR_J', data=array_2mass_j)
grp_lum.create_dataset('SFR_W3', data=array_wise_w3)
grp_lum.create_dataset('SFR_UVIR', data=array_UVIR)
grp_lum.create_dataset('SFR_Ha', data=array_Ha)

f.close()
