# coding: utf-8

#!/usr/bin/env python2

"""
Created on October 27, 2018

@author: sandro.tacchella@cfa.harvard.edu

"""

import numpy as np
import h5py
import glob


def read_in_files(path_files):
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
        list_alpha = np.append(list_alpha, float(ii_file.split('/')[-1].split('_')[2]))
        list_tau = np.append(list_tau, float(ii_file.split('/')[-1].split('_')[-1][:-5]))
    # create dictionary
    dict_data = {}
    for ii_alpha in np.unique(list_alpha):
        dict_sub = {}
        for jj_tau in np.unique(list_tau):
            dict_sub[jj_tau] = h5py.File(list_all_files[(list_alpha == ii_alpha) & (list_tau == jj_tau)][0], 'r')
        dict_data[ii_alpha] = dict_sub
    return(dict_data, list_alpha, list_tau)


