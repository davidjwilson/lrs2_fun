import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const

"""
Quick script for coadding LRS2 data and saving each band as ecsv files. Asumes that all of your spectra has been dumped into one place. Coadds all data for each star. No date dependence yet. 

Usage - last line, adjust the lrs2_add() function. Defaults are set to look for the fits files in a directory called "data", and save them to a directory called "spectra"

lrs2_add(path='data/', plot=True, save=True, bands = ['uv', 'orange', 'red', 'farred'], savepath='spectra/')

"""

def coadd_flux(f_array, e_array):
    """
    Returns a variance-weighted coadd with standard error of the weighted mean (variance weights, scale corrected).
    f_array and e_arrays are collections of flux and error arrays, which should have the same lenth and wavelength scale
    """
    weights = 1 / (e_array**2)
    flux = np.average(f_array, axis =0, weights = weights)
    var = 1 / np.sum(weights, axis=0)
    rcs = np.sum((((flux - f_array)**2) * weights), axis=0) / (len(f_array)-1) #reduced chi-squared
    error = (var * rcs)**0.5
    return flux, error


def save_ecsv(w, f, e, star, band, savepath):
    """
    saves coadded spectrum to ecsv
    """
    save_dat = Table([w*u.AA,f*u.erg/u.s/u.cm**2/u.AA, e*u.erg/u.s/u.cm**2/u.AA], names = ['WAVELENGTH', 'FLUX', 'ERROR'])
    name = savepath + star+'_LRS2_'+band+'.ecsv'
    ascii.write(save_dat, name, format='ecsv', overwrite=True)
    
def coadd_by_star(path, star, band):
    """
    combines all spectra for a given star and band
    """
    fs = []
    es = []
    spectra = glob.glob(path+'*'+band+'.fits')
    for sp in spectra:
        if fits.getheader(sp,0)['OBJECT'] == star:
            data = fits.getdata(sp,0)
            fs.append(data[1])
            es.append(data[3]) #2 is sky
            w = data[0] #fixed w scale, no need to interpolate
    f, e = coadd_flux(np.array(fs),np.array(es))
    return w, f, e

def plot_spectrum(ws, fs, star, figsave, figsavepath):
    """
    plots all bands for a star. ws and fs are collections of w and f arrays
    """
    plt.figure(star)
    for w, f in zip(ws, fs):
        plt.plot(w,f)
    plt.xlabel('Wavelength (\AA)', size=20)
    plt.ylabel('Flux (erg s$^{-1}$cm$^{-2}$\AA$^{-1}$)', size=20)
    plt.tight_layout()
    if figsave:
        plt.savefig(figsavepath+star+'_lrs2.png', dpi=150)
        plt.savefig(figsavepath+star+'_lrs2.pdf', dpi=150)
    plt.show()
    plt.close()

def stars_finder(path):
    """
    identifies all targets in a file 
    """
    stars  = []
    spectra = os.listdir(path)
    for spectrum in spectra:
        stars.append(fits.getheader(path+spectrum,0)['OBJECT'])
    return np.unique(stars)
    
def lrs2_add(path='data/', plot=True, save=True, bands = ['uv', 'orange', 'red', 'farred'], savepath='spectra/', figsave=True, figsavepath='plots/'):
    """
    Main function. Does not use dates yet.
    """
    stars = stars_finder(path)
    for star in stars:
        print(star)
        ws, fs, es = [], [], []
        for band in bands:
            print(band)
            w, f, e = coadd_by_star(path, star, band)
            if save:
                save_ecsv(w, f, e, star, band, savepath)
            ws.append(w)
            fs.append(f)
            es.append(e)
        if plot:
            plot_spectrum(ws, fs, star, figsave, figsavepath)

                


bands = ['uv', 'orange'] #we only used the blue arm so far
lrs2_add(bands=bands)