import numpy as np
import healpy as hp
import h5py
import re
from scipy.interpolate import interp1d

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

try :
  lin_data.shape
except :
  lindir = '/home/jmertens/cosmograph/build/dust_lensing_R-128_F-800_kcut-17_Nside-16_S-7/job-11582/dust_lensing/'
  f_lin = h5py.File(lindir+'raytracedata.values.h5.gz', 'r')
  lin_data = []
  for k in natural_sort(list(f_lin.keys())) :
    try :
      lin_data.append(np.array(f_lin[k]))
    except :
      print("Unable to use data with key:", k)
lin_data = np.array(lin_data)

try :
  nl_data.shape
except :
  nldir = '/home/jmertens/cosmograph/build/dust_lensing_R-128_F-800_kcut-17_Nside-16_S-7/job-11560/dust_lensing/'
  f_nl = h5py.File(nldir+'raytracedata.values.h5.gz', 'r')
  nl_data = []
  for k in natural_sort(list(f_nl.keys())) :
    try :
      nl_data.append(np.array(f_nl[k]))
    except :
      print("Unable to use data with key:", k)
  nl_data = np.array(nl_data)

def zofE(E) :
  return E - 1

rhoz_lin_data = [interp1d(zofE(lin_data[:-2,n,0]), lin_data[:-2,n,8]) for n in range(lin_data.shape[1]-1)]
rhoz_nl_data = [interp1d(zofE(nl_data[:-2,n,0]), nl_data[:-2,n,8]) for n in range(nl_data.shape[1]-1)]

def rho(data, z) :
  return np.array([data[n](z) for n in range(len(data))])

rhos_lin_z1 = rho(rhoz_lin_data, .1)
rhos_nl_z1 = rho(rhoz_nl_data, .1)

cls_lin = hp.anafast(rhos_lin_z1)
cls_nl = hp.anafast(rhos_nl_z1)

