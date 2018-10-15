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
import pickle
import argparse
import make_SFH_from_PSD


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
parser.add_argument("--number_of_galaxies", type=int, help="Number of SFHs/galaxies")
parser.add_argument("--scatter_MS_0", type=float, help="scatter of MS (on shortest timescale)")
parser.add_argument("--list_of_slopes", type=list, help="slope on the short time scales")
parser.add_argument("--list_of_breaks", type=list, help="timescale of break")
parser.add_argument("--aliasTbin", type=float, help="oversampling of time series (resolution)")
args = parser.parse_args()


# get grid of times

time_SFH = np.arange(cosmo.age(args.redshift_start).value, cosmo.age(args.redshift_end).value, args.sfh_res)
number_steps = len(time_SFH)  # each step coresponds to 1.0 Myr, so total is 10 Gyr
print('number of time steps: ' + str(number_steps))
print('number of galaxies: ' + str(args.number_galaxies))


### set parameters of light curve model ###
# amplitude => seems not to work
amp = args.scatter_MS_0

# slope on the short time scales
list_of_slopes = np.array(args.list_of_slopes)
list_of_slopes_str = list_of_slopes.astype('str')

#slope on the long time scales
slope_high = 0.0

#frequency of the break
list_of_break_timescales = np.array(args.list_of_breaks)
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
array_of_SFH = make_SFH_from_PSD.create_family_SFHs(number_galaxies, list_of_slopes[ii_slope], list_of_breaks[jj_break], args.scatter_MS_0, args.aliasTbin, number_steps, amp, slope_high)

# make dictionary
SFH_dict = {}
SFH_dict['time'] = time_SFH
SFH_dict['SFR_array'] = array_of_SFH

# save SFH dictionary (contains SFH)
file_name = path_SFH_cat + args.filename_SFH + 'alpha_' + list_of_slopes_str[idx_slope] + '_tau_' + list_of_breaks_str[idx_break] + '.hdf5'
with open(file_name, 'wb') as f:
    pickle.dump(SFH_dict, f, pickle.HIGHEST_PROTOCOL)




