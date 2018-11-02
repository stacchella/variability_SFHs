# coding: utf-8

#!/usr/bin/env python2

"""
Created on October 17, 2018

@author: sandro.tacchella@cfa.harvard.edu

"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import simps

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
print cosmo

def integrate_mag(spec_lam, spectra, filter, z=None):
    '''
    borrowed from calc_ml
    given a filter name and spectrum, calculate magnitude/luminosity in filter (see alt_file for filter names)
    INPUT:
        SPEC_LAM: must be in angstroms. this will NOT BE corrected for reddening even if redshift is specified. this
        allows you to calculate magnitudes in rest or observed frame.
        SPECTRA: must be in Lsun/Hz (FSPS standard). if redshift is specified, the normalization will be taken care of.
    OUTPUT:
        MAG: comes out as absolute magnitude
        LUMINOSITY: comes out in erg/s
            NOTE: if redshift is specified, INSTEAD RETURN apparent magnitude and flux [erg/s/cm^2]
    '''
    resp_lam = filter[0][0]
    res = filter[1][0]
    # physical units, in CGS
    pc2cm = 3.08568E18
    lsun = 3.839E33
    c = 2.99E10
    # interpolate filter response onto spectral lambda array
    # when interpolating, set response outside wavelength range to be zero.
    response_interp_function = interp1d(resp_lam, res, bounds_error=False, fill_value=0)
    resp_interp = response_interp_function(spec_lam)
    # integrate spectrum over filter response
    # first calculate luminosity: convert to flambda (factor of c/lam^2, with c converted to AA/s)
    # then integrate over flambda [Lsun/AA] to get Lsun
    spec_flam = spectra*(c*1e8/(spec_lam**2))
    luminosity = simps(spec_flam*resp_interp, spec_lam)
    # now calculate luminosity density [erg/s/Hz] in filter
    # this involves normalizing the filter response by integrating over wavelength
    norm = simps(resp_interp/spec_lam, spec_lam)
    luminosity_density = simps(spectra*(resp_interp/norm)/spec_lam, spec_lam)
    # if redshift is specified, convert to flux and apparent magnitude
    if z is not None:
        dfactor = (cosmo.luminosity_distance(z).value*1e5)**(-2)*(1+z)
        luminosity = luminosity*dfactor
        luminosity_density = luminosity_density*dfactor
    # convert luminosity density to flux density
    # the units of the spectra are Lsun/Hz; convert to
    # erg/s/cm^2/Hz, at 10pc for absolute mags
    flux_density = luminosity_density*lsun/(4.0*np.pi*(pc2cm*10)**2)
    luminosity = luminosity*lsun
    # convert flux density to magnitudes in AB system
    mag = -2.5*np.log10(flux_density)-48.60
    #print 'maggies: {0}'.format(10**(-0.4*mag)*1e10)
    return mag, luminosity


def return_lir(lam, spec, z=None):
    """ returns IR luminosity (8-1000 microns) in erg/s
    input spectrum must be Lsun/Hz, wavelength in \AA
    """
    botlam = np.atleast_1d(8e4-1)
    toplam = np.atleast_1d(1000e4+1)
    edgetrans = np.atleast_1d(0)
    lir_filter = [[np.concatenate((botlam, np.linspace(8e4, 1000e4, num=100), toplam))],
                  [np.concatenate((edgetrans, np.ones(100), edgetrans))]]
    # calculate integral
    _, lir = integrate_mag(lam, spec, lir_filter, z=z)  # comes out in ergs/s
    return lir


def return_luv(lam, spec, z=None):
    """ returns UV luminosity (1216-3000 \AA) in erg/s
    input spectrum must be Lsun/Hz, wavelength in \AA
    """
    botlam = np.atleast_1d(1216)
    toplam = np.atleast_1d(3000)
    edgetrans = np.atleast_1d(0)
    luv_filter = [[np.concatenate((botlam-1, np.linspace(botlam, toplam, num=100), toplam+1))],
                  [np.concatenate((edgetrans, np.ones(100), edgetrans))]]
    # calculate integral
    _, luv = integrate_mag(lam, spec, luv_filter, z=z)  # comes out in ergs/s
    return luv


def get_magnitude_SFH(sp_in, time_SFH, SFR_SFH, tage_list, filters, idx_Halpha):
    for tage in tage_list:
        if (tage == tage_list[0]):
            L_mat = get_magnitude_SFH_single(sp_in, time_SFH, SFR_SFH, tage, filters, idx_Halpha)
        else:
            L_mat = np.vstack([L_mat, get_magnitude_SFH_single(sp_in, time_SFH, SFR_SFH, tage, filters, idx_Halpha)])
    return(L_mat.T)


def get_magnitude_SFH_single(sp_in, time_SFH, SFR_SFH, tage_in, filters, idx_Halpha):
    sp_in.set_tabular_sfh(time_SFH, SFR_SFH)
    # get luminosities
    mag_list = sp_in.get_mags(tage=tage_in, bands=filters)
    L_list = 4*np.pi*(3.086e+19)**2*np.power(10, -0.4*(mag_list+48.6))
    wave, spec = sp_in.get_spectrum(tage=tage_in)
    LIR = return_lir(wave, spec, z=None)
    L_tot = LIR + 2.2*return_luv(wave, spec, z=None)  # from total UV
    L_tot2 = LIR + 2.2*1.5*L_list[1]*3e8/(2800*10**-10)  # from 2800
    L_list = np.append(np.append(np.append(L_list, L_tot), L_tot2), 3.839*10**33*sp_in.emline_luminosity[idx_Halpha])
    return(L_list)

