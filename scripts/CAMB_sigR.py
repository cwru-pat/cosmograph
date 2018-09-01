# Referencing: http://home.strw.leidenuniv.nl/~jarle/Teaching/GalaxyFormation/Lectures/lecture04.pdf
# And https://camb.readthedocs.io/en/latest/CAMBdemo.html

import numpy as np
import camb
from scipy import integrate
from camb import model

#Now get matter power spectra and sigma8 at redshift 0 and 0.8
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.set_dark_energy() #re-set defaults
pars.InitPower.set_params(ns=0.965)
#Not non-linear corrections couples to smaller scales than you want
pars.NonLinear = model.NonLinear_both

Om0 = 0.316049
zI = 50.0
aI = 1.0
a0 = (1.0 + zI)*aI
H0 = 4440 # Mpc^-1
HI = H0*np.sqrt( Om0 * (aI/a0)**3 )
OL0 = 1.0 - Om0
OLI = OL0/( OL0 + Om0*(a0/aI)**3 )
print("OLI = ", OLI)

Ls = np.logspace(np.log10(0.1), np.log10(12), 20)
LphysIs = HI*Ls
Lphys0s = LphysIs*a0/aI


pars.set_matter_power(redshifts=[zI], kmax=2.0)
results = camb.get_results(pars)
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
ks = kh_nonlin/0.675
Dks = pk_nonlin[0,:]*1.0/2.0/np.pi*ks**3

def W(kR) :
    return np.exp(-(kR)**2/2)

def getSig(R):
    sig2 = integrate.trapz(Dks*W(ks*R)**2/ks, ks)
    return np.sqrt(sig2)

Sigs = np.array([getSig(R) for R in Lphys0s])
print(Ls)
print(Sigs)
print(Lphys0s)

