# coding: utf-8

#!/usr/bin/env python2

"""
Created on October 27, 2018

@author: sandro.tacchella@cfa.harvard.edu

"""

import numpy as np
import h5py
import glob


def read_in_files(path_files, mode='r'):
    '''
    This function reads all hdf5 files for
    one SFH run. It returns a dictionary with
    all hdf5 files.
    '''
    # get list of all files
    list_all_files = np.array(glob.glob(path_files))
    # create list for alpha and tau
    list_alpha = []
    list_tau = []
    for ii_file in list_all_files:
        list_alpha = np.append(list_alpha, float(ii_file.split('/')[-1].split('_')[-3]))
        list_tau = np.append(list_tau, float(ii_file.split('/')[-1].split('_')[-1][:-5]))
    # create dictionary
    dict_data = {}
    for ii_alpha in np.unique(list_alpha):
        dict_sub = {}
        for jj_tau in np.unique(list_tau):
            dict_sub[jj_tau] = h5py.File(list_all_files[(list_alpha == ii_alpha) & (list_tau == jj_tau)][0], mode)
        dict_data[ii_alpha] = dict_sub
    return(dict_data, list_alpha, list_tau)


#### add SFRs from different indicators ###

# use Driver+15 conversions
# log10[SFR(M_sun yr^−1)] = m(log10[L(unit)] − Lf) + C.
conversion_type = np.array(['Ha', 'FUV', 'NUV', 'u', 'W3', 'UV+TIR'])
m_list = np.array([0.59, 0.75, 0.62, 0.92, 0.66, 0.83])
Lf_list = np.array([34.0, 21.5, 21.5, 21.25, 22.25, 37.0])
C_list = np.array([0.031, 0.17, 0.014, -0.079, 0.16, 0.20])

# SFR indicators considered
name_indicator = ['Ha', 'FUV', 'NUV', 'u', 'W3', 'UV+TIR', 'UV+TIR']
name_lum = ['lum/lum_Ha', 'lum/lum_i1500', 'lum/lum_i2800', 'lum/lum_u', 'lum/lum_wise_w3', 'lum/lum_UVIR', 'lum/lum_UVIR2']
name_SFR = ['SFR_Ha', 'SFR_FUV', 'SFR_NUV', 'SFR_u', 'SFR_W3', 'SFR_UVIR', 'SFR_UVIR2']


def get_SFR(SF_type, lum):
    idx_type = (conversion_type == SF_type)
    return(m_list[idx_type]*(np.log10(lum) - Lf_list[idx_type]) + C_list[idx_type])


def append_SFR(dict_data, list_alpha, list_tau):
    for ii_slope in np.unique(list_alpha):
        print 'working on slope = ', ii_slope
        for jj_break in np.unique(list_tau):
            # save SFRs
            f = dict_data[ii_slope][jj_break]
            try:
                del f["SFR"]
            except KeyError:
                pass
            grp_SFR = f.create_group("SFR")
            grp_SFR.create_dataset('SFR_time', data=dict_data[ii_slope][jj_break]['lum/lum_time'][:])
            # compute SFRs
            for ii in range(len(name_indicator)):
                SFR_data = 10**get_SFR(name_indicator[ii], dict_data[ii_slope][jj_break][name_lum[ii]][:])
                SFR_data = np.log10(SFR_data/np.nanmedian(SFR_data))
                grp_SFR.create_dataset(name_SFR[ii], data=SFR_data)
            f.close()

