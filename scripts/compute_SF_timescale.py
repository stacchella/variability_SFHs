# read in modules

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy import signal
from astropy.table import Table
from DELCgen import *
import nolds
import pickle
from scipy.ndimage.filters import gaussian_filter


import toolbox_variability
import MS_Variability


# define paths

path_data = '/Volumes/Tacchella/Work/Postdoc/variability/catalogs/SFH/'
path_figures = '/Users/sandrotacchella/ASTRO/SFH_Variability/variability_SFHs/figures/'
path_analytical = '/Users/sandrotacchella/ASTRO/SFH_Variability/variability_SFHs/data/'


# read in files

dict_data, list_alpha, list_tau = toolbox_variability.read_in_files(path_data + 'SFH_04_*.hdf5', mode='r')


# create unique list

list_of_slopes = np.unique(list_alpha)
list_of_tau = np.unique(list_tau)

print list_of_slopes
print list_of_tau


# average SFR over time
# which time does avg SFR agree with indicator
# same for scatter: which time does scatter of averaged SFR agree with scatter of indicator


def get_delta_indicator(slope, tau, ave_time_list, indicator, num_rand_samples=1):
    random_times = np.random.randint(0, len(dict_data[slope][tau]['SFH/time'][:]), num_rand_samples)
    delta_SFR_P16 = []
    delta_SFR_P50 = []
    delta_SFR_P84 = []
    for dt in ave_time_list:
        delta_SFR_sublist = []
        for time_0 in random_times[random_times > dt]:
            # compute SFR difference
            SFR = np.mean(10**dict_data[slope][tau]['SFH/SFR'][:][:, time_0-dt:time_0], axis=1)
            SFR_indicator = 10**dict_data[slope][tau][indicator][:][:, time_0]
            delta_SFR_sublist.append(np.abs(np.median(np.log10(SFR)-np.log10(SFR_indicator))))
        delta_SFR_P16.append(np.percentile(delta_SFR_sublist, 16))
        delta_SFR_P50.append(np.percentile(delta_SFR_sublist, 50))
        delta_SFR_P84.append(np.percentile(delta_SFR_sublist, 84))
    return([delta_SFR_P16, delta_SFR_P50, delta_SFR_P84])


def get_avg_timescale(slope, tau, ave_time_list, indicator, factor_range, num_rand_samples):
    delta_SFR_list = get_delta_indicator(slope, tau, ave_time_list, indicator, num_rand_samples)
    timescale_SFR = ave_time_list[np.argmin(delta_SFR_list[1])]
    timescale_SFR_low = np.min(ave_time_list[factor_range*np.min(delta_SFR_list[1]) > delta_SFR_list[1]])
    timescale_SFR_high = np.max(ave_time_list[factor_range*np.min(delta_SFR_list[1]) > delta_SFR_list[1]])
    return([timescale_SFR, timescale_SFR_low, timescale_SFR_high])


def get_avg_timescale_profiles(indicator, text_indicator, ave_time_list, num_rand_samples, xlim_in, factor_range=1.2):
    timescale_SFR = []
    timescale_SFR_range =[]
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10.0, 5.0))
    for ii_slope in [1.0, 2.0, 3.0]:
        for jj_tau in [10.0, 100.0, 1000.0]:
            delta_SFR_list = get_delta_indicator(ii_slope, jj_tau, ave_time_list, indicator, num_rand_samples=num_rand_samples)
            ax1.plot(ave_time_list, delta_SFR_list[1], '-', lw=2, label=r'$\alpha=%.1f,\ \tau=%.0f$' %(ii_slope, jj_tau))
            ax1.fill_between(ave_time_list, delta_SFR_list[0], delta_SFR_list[2], alpha=0.4)
            timescale_SFR = np.append(timescale_SFR, ave_time_list[np.argmin(delta_SFR_list[1])])
            timescale_SFR_range = np.append(timescale_SFR_range, [np.min(ave_time_list[factor_range*np.min(delta_SFR_list[1]) > delta_SFR_list[1]]), np.max(ave_time_list[factor_range*np.min(delta_SFR_list[1]) > delta_SFR_list[1]])])
    #ax1.text(30, 0.7, text_indicator, fontsize=16)
    ax1.set_xlim(xlim_in)
    #ax1.set_ylim([-0.02, 0.8])
    ax1.set_xlabel('avg time', fontsize=18)
    #ax1.xaxis.set_ticklabels([])
    ax1.set_ylabel(r'$\Delta_{\rm SFR}$', fontsize=18)
    ax1.legend(frameon=False, ncol=2)
    plt.subplots_adjust(hspace=0.05)
    #plt.savefig(path_Figs + 'SFR_about_MS_b ' + list_sample_SFHs_str_simple[ii_SFH] + '.pdf', bbox_inches='tight', dpi=400)
    plt.show()
    return(timescale_SFR, timescale_SFR_range)


def get_avg_timescale_matrix(ave_time_list, indicator, factor_range, num_rand_samples):
    list_of_slopes_plot = list_of_slopes[2:][::2]
    list_of_tau_plot = list_of_tau[::2]
    t_mat_SFR = np.zeros((len(list_of_slopes_plot), len(list_of_tau_plot)))
    t_mat_SFR_l = np.zeros((len(list_of_slopes_plot), len(list_of_tau_plot)))
    t_mat_SFR_h = np.zeros((len(list_of_slopes_plot), len(list_of_tau_plot)))
    for ii_slope in range(len(list_of_slopes_plot)):
        print 'progress (%) = ', np.round(100.0*ii_slope/len(list_of_slopes_plot), 0)
        for jj_tau in range(len(list_of_tau_plot)):
            t_SFR = get_avg_timescale(list_of_slopes_plot[ii_slope], list_of_tau_plot[jj_tau], ave_time_list, indicator, factor_range=factor_range, num_rand_samples=num_rand_samples)
            t_mat_SFR[ii_slope][jj_tau] = t_SFR[0]
            t_mat_SFR_l[ii_slope][jj_tau] = t_SFR[1]
            t_mat_SFR_h[ii_slope][jj_tau] = t_SFR[2]
    return(t_mat_SFR, t_mat_SFR_l, t_mat_SFR_h)


ave_time_list = np.unique(np.logspace(0.0, 3.0, num=50).astype('int'))
print ave_time_list


factor_range = 1.2
num_rand_samples = 100

t_mat_SFR_Ha = get_avg_timescale_matrix(ave_time_list, 'SFH_conv/SFR_Ha', factor_range=factor_range, num_rand_samples=num_rand_samples)
#t_mat_SFR_UVIR = get_avg_timescale_matrix(ave_time_list, 'SFH_conv/SFR_UVIR', factor_range=factor_range, num_rand_samples=num_rand_samples)
#t_mat_SFR_u = get_avg_timescale_matrix(ave_time_list, 'SFH_conv/SFR_u', factor_range=factor_range, num_rand_samples=num_rand_samples)
#t_mat_SFR_FUV = get_avg_timescale_matrix(ave_time_list, 'SFH_conv/SFR_FUV', factor_range=factor_range, num_rand_samples=num_rand_samples)
#t_mat_SFR_NUV = get_avg_timescale_matrix(ave_time_list, 'SFH_conv/SFR_NUV', factor_range=factor_range, num_rand_samples=num_rand_samples)
#t_mat_SFR_W3 = get_avg_timescale_matrix(ave_time_list, 'SFH_conv/SFR_W3', factor_range=factor_range, num_rand_samples=num_rand_samples)


# save matrix

np.save(path_data + 't_mat_SFR_Ha', t_mat_SFR_Ha)
#np.save(path_data + 't_mat_SFR_UVIR', t_mat_SFR_UVIR)
#np.save(path_data + 't_mat_SFR_u', t_mat_SFR_u)
#np.save(path_data + 't_mat_SFR_FUV', t_mat_SFR_FUV)
#np.save(path_data + 't_mat_SFR_NUV', t_mat_SFR_NUV)
#np.save(path_data + 't_mat_SFR_W3', t_mat_SFR_W3)
