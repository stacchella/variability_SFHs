# coding: utf-8

#!/usr/bin/env python2

"""
Created on October 14, 2018

@author: sandro.tacchella@cfa.harvard.edu

"""

import numpy as np
from DELCgen import *


def generate_SFH(a_low, v_bend, aliasTbin, number_steps, amp, a_high):
    '''
    generate random walk SFH with slope of the
    power-spectrum density (a_low) and length
    '''
    return Simulate_TK_Lightcurve(BendingPL, (amp, v_bend, a_low, a_high, 0), RedNoiseL=10, aliasTbin=aliasTbin, tbin=1, length=number_steps)


def get_stats_DEL(SFH_array, number_galaxies, only_last=False):
    flux_all = []
    std_individual_list = []
    for ii in range(number_galaxies):
        if only_last:
            std_individual_list.append(np.std(SFH_array[ii].flux[-1]))
            flux_all.append(SFH_array[ii].flux[-1])
        else:
            std_individual_list.append(np.std(SFH_array[ii].flux))
            flux_all.append(SFH_array[ii].flux)
    flux_all = np.array(flux_all).flatten()
    return(np.mean(flux_all), np.std(flux_all), np.mean(std_individual_list), np.std(std_individual_list))


def create_family_SFHs(number_galaxies, slope, v_bend, scatter_MS, aliasTbin, number_steps, amp, a_high):
    print 'making SFHs...'
    array_of_DELs = np.array([generate_SFH(slope, v_bend, aliasTbin, number_steps, amp, a_high) for ii in range(number_galaxies)])
    mean, std, std_mean_individual, std_std_individual = get_stats_DEL(array_of_DELs, number_galaxies, only_last=True)
    # get normalization so that scatter is consistent with MS
    normalization = scatter_MS/std
    print '#####################'
    print 'slope = ', slope
    print 'slope = ', v_bend
    print 'mean(global) = ', np.round(mean, 3)
    print 'std(global) = ', np.round(std, 3)
    print 'std(individual) = ', np.round(std_mean_individual, 3), ' +- ', np.round(std_std_individual, 3)
    for ii in range(number_galaxies):
        if (ii == 0):
            array_of_SFH = normalization*(array_of_DELs[ii].flux-mean)
        else:
            array_of_SFH = np.vstack([array_of_SFH, normalization*(array_of_DELs[ii].flux-mean)])
    return(array_of_SFH)




