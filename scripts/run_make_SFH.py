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
import fsps
import compute_luminosities


from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
print cosmo


# define working directory

path_main = os.environ['WDIR_SFH_variability']
path_SFH_cat = path_main + 'catalogs/SFH/'


# read in command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("--idx_key", type=int, help="iteration variable")
parser.add_argument("--filename_SFH", type=str, help="filename of SFH file")
parser.add_argument("--sfh_res", type=float, help="time resolution of SFH in Gyr")
parser.add_argument("--redshift_start", type=float, help="start of SFH in redshift space")
parser.add_argument("--redshift_end", type=float, help="end of SFH in redshift space")
parser.add_argument("--redshift_lum_start", type=float, help="start of luminosity computation in redshift space")
parser.add_argument("--redshift_lum_end", type=float, help="end of luminosity computation in redshift space")
parser.add_argument("--number_galaxies", type=int, help="Number of SFHs/galaxies")
parser.add_argument("--scatter_MS_0", type=float, help="scatter of MS (on shortest timescale)")
parser.add_argument("--list_of_slopes", nargs='+', type=float, help="slope on the short time scales")
parser.add_argument("--list_of_breaks", nargs='+', type=float, help="timescale of break")
parser.add_argument("--aliasTbin", type=float, help="oversampling of time series (resolution)")
parser.add_argument("--logzsol", type=float, help="metallicity: logZsol in FSPS")
parser.add_argument("--dust2", type=float, help="dust2 in FSPS")

args = parser.parse_args()


# get grid of times

time_cosmo = np.arange(cosmo.age(args.redshift_start).value, cosmo.age(args.redshift_end).value, args.sfh_res)
number_steps = len(time_cosmo)
print('number of time steps: ' + str(number_steps))
print('number of galaxies: ' + str(args.number_galaxies))

time_SFH = time_cosmo-time_cosmo[0]
print time_SFH
print len(time_SFH)

# get time of luminosity

t_start = cosmo.age(args.redshift_lum_start).value-time_cosmo[0]
t_end = cosmo.age(args.redshift_lum_end).value-time_cosmo[0]
idx_t = (time_SFH >= t_start) & (time_SFH <= t_end)
time_lum = time_SFH[idx_t]

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
list_of_breaks = 1.0/(list_of_break_timescales*2.0*np.pi)
list_of_breaks_str = list_of_break_timescales.astype('str')

# intercept/offset: adds constant to the final output
# c = 10.0


# define function for translating key

def translate_key(idx):
    '''
    Translate key to indices.
    '''
    idx_slope = idx / len(list_of_breaks)
    idx_break = idx % len(list_of_breaks)
    return(idx_slope, idx_break)


# generate SFH from power-spectrum density
ii_slope, jj_break = translate_key(args.idx_key)
array_of_SFH = make_SFH_from_PSD.create_family_SFHs(args.number_galaxies, list_of_slopes[ii_slope], list_of_breaks[jj_break], args.scatter_MS_0, args.aliasTbin, number_steps, amp, slope_high)


# compute luminosities

# set up fsps
sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1, imf_type=1, add_neb_emission=True, sfh=3, logzsol=args.logzsol, dust_type=2, dust2=args.dust2)


# define filters
filters = ['i1500', 'i2800', 'u', 'v', '2mass_j', 'wise_w3']
idx_Halpha = (np.abs(sp.emline_wavelengths-6564.61)).argmin()


array_of_luminosities = []

for ii in range(number_galaxies):
    SFR = 10**array_of_SFH[ii]  # assume MS is 1 Msun/yr
    if (ii == 0):
        array_of_luminosities = compute_luminosities.get_magnitude_SFH(sp, time_SFH, SFR, time_lum, filters, idx_Halpha)
    else:
        array_of_luminosities = np.vstack([array_of_luminosities, compute_luminosities.get_magnitude_SFH(sp, time_SFH, SFR, time_lum, filters, idx_Halpha)])


# save SFH dictionary (contains SFH)
file_name = path_SFH_cat + args.filename_SFH + 'alpha_' + list_of_slopes_str[ii_slope] + '_tau_' + list_of_breaks_str[jj_break] + '.hdf5'

try:
    os.remove(file_name)
except OSError:
    pass

f = h5py.File(file_name, 'w')
# add SFH
grp_SFH = f.create_group("SFH")
grp_SFH.create_dataset('SFH_time_cosmo', data=time_cosmo)
grp_SFH.create_dataset('SFH_time', data=time_SFH)
grp_SFH.create_dataset('SFH_SFR', data=array_of_SFH)
grp_lum = f.create_group("lum")
grp_lum.create_dataset('lum_time', data=time_lum)
grp_lum.create_dataset('lum_luminosities', data=array_of_luminosities)

f.close()

