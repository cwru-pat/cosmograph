import numpy as np
import healpy as hp
import re
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import colors as mcolors

NSIDE=32
NPIX=12*NSIDE**2

try :
  len(DATA)
except :
  DATA={}

def loadData(dir, scale=1.0) :
  global DATA
  try :
    DATA[dir].shape
  except :
    try :
      print("Loading npz data from", dir)
      data = np.load(dir+'raysheet.dat.npz')
      print("Done.")
    except :
      print("Loading raw data from", dir)
      data = np.loadtxt(dir+"raysheet.dat.gz")
      print("Done.")
      data = np.reshape(data, (data.shape[0],NPIX,3))
      if scale != 1.0 :
        # Scale up linear/small-amplitude data
        for i in range(3) :
          data[:,:,i] = np.array([
            scale*(data[step,:,i]-np.mean(data[step,:,i]))+np.mean(data[step,:,i])
            for step in range(data.shape[0]) ])
      np.savez(dir+'raysheet.dat.npz', E=data[:,:,0], rho=data[:,:,1], DA=data[:,:,2])
      data = np.load(dir+'raysheet.dat.npz')
    DATA[dir] = data
  return DATA[dir]

def rescaleData(dir, scale=1.0/1000.0) :
  data = np.load(dir+'raysheet.dat.npz')
  print("Working in", dir)
  rho_data = data['rho']
  E_data = data['E']
  rho_data = np.array([scale*(rho_data[step,:]-np.mean(rho_data[step,:]))+np.mean(rho_data[step,:])
    for step in range(rho_data.shape[0]) ])
  E_data = np.array([scale*(E_data[step,:]-np.mean(E_data[step,:]))+np.mean(E_data[step,:])
    for step in range(E_data.shape[0]) ])
  np.savez(dir+'raysheet.dat.npz', rho=rho_data, E=E_data)


def zofE(E) :
  return E/E[0] - 1

def rho_interp(data) :
  print("Loading E data...")
  E = data['E']
  print("Loading rho data...")
  rho = data['rho']
  print("Interpolating...")
  interps = [interp1d(zofE(E[:,n]), rho[:,n]) for n in range(rho.shape[1])]
  print("Done.")
  return interps

def rho(rho_interp_fn, z) :
  return np.array([rho_interp_fn[n](z) for n in range(len(rho_interp_fn))])




if False :
  lindir_p_96 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r96_P0.00000013/dust_lensing/'
  rho_data_lin_96 = loadData(lindir_p_96, 1000.0)
  rho_interp_lin_96 = rho_interp(rho_data_lin_96)
  rhos_lin_96 = rho(rho_interp_lin_96, 0.25)
  cls_lin_96 = hp.anafast(rhos_lin_96)

  nldir_p_96 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r96_P0.13/dust_lensing/'
  rho_data_nl_96 = loadData(nldir_p_96, 1.0)
  rho_interp_nl_96 = rho_interp(rho_data_nl_96)
  rhos_nl_96 = rho(rho_interp_nl_96, 0.25)
  cls_nl_96 = hp.anafast(rhos_nl_96)


  lindir_p_128 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r128_P0.00000013/dust_lensing/'
  rho_data_lin_128 = loadData(lindir_p_128, 1000.0)
  rho_interp_lin_128 = rho_interp(rho_data_lin_128)
  rhos_lin_128 = rho(rho_interp_lin_128, 0.25)
  cls_lin_128 = hp.anafast(rhos_lin_128)

  nldir_p_128 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r128_P0.13/dust_lensing/'
  rho_data_nl_128 = loadData(nldir_p_128, 1.0)
  rho_interp_nl_128 = rho_interp(rho_data_nl_128)
  rhos_nl_128 = rho(rho_interp_nl_128, 0.25)
  cls_nl_128 = hp.anafast(rhos_nl_128)


  lindir_p_144 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r144_P0.00000013/dust_lensing/'
  rho_data_lin_144 = loadData(lindir_p_144, 1000.0)
  rho_interp_lin_144 = rho_interp(rho_data_lin_144)
  rhos_lin_144 = rho(rho_interp_lin_144, 0.25)
  cls_lin_144 = hp.anafast(rhos_lin_144)

  nldir_p_144 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r144_P0.13/dust_lensing/'
  rho_data_nl_144 = loadData(nldir_p_144, 1.0)
  rho_interp_nl_144 = rho_interp(rho_data_nl_144)
  rhos_nl_144 = rho(rho_interp_nl_144, 0.25)
  cls_nl_144 = hp.anafast(rhos_nl_144)

  diff96 = ((cls_nl_96-cls_lin_96)/cls_lin_96)[:17]
  diff128 = ((cls_nl_128-cls_lin_128)/cls_lin_128)[:17]
  diff144 = ((cls_nl_144-cls_lin_144)/cls_lin_144)[:17]

def getCls(dir, mlt, smooth=0.0) :
  print("Working on data in", dir)
  rho_data = loadData(dir, mlt)
  rho_interps = rho_interp(rho_data)
  rho_map = rho(rho_interps, 0.25)
  if smooth:
    rho_map = hp.sphtfunc.smoothing(rho_map, smooth)
  Cls = hp.anafast(rho_map)
  return Cls

def getalms(dir, mlt, smooth=0.0) :
  print("Working on data in", dir)
  rho_data = loadData(dir, mlt)
  rho_interps = rho_interp(rho_data)
  rho_map = rho(rho_interps, 0.25)
  if smooth:
    rho_map = hp.sphtfunc.smoothing(rho_map, smooth)
  alms = hp.map2alm(rho_map)
  return alms

lin_dirs = [ "r96_P0.00000013_M15", "r96_P-0.00000013_M15", "r96_P0.00000013_M14", "r96_P-0.00000013_M14", "r96_P0.00000013_M13", "r96_P-0.00000013_M13", "r96_P0.00000013_M12", "r96_P-0.00000013_M12", "r96_P0.00000013_M11", "r96_P-0.00000013_M11", "r96_P0.00000013_M10", "r96_P-0.00000013_M10", "r96_P0.00000013_M9", "r96_P-0.00000013_M9", "r96_P0.00000013_M8", "r96_P-0.00000013_M8", "r96_P0.00000013", "r96_P-0.00000013" ]
nl_dirs = [ "r96_P0.13_M15", "r96_P-0.13_M15", "r96_P0.13_M14", "r96_P-0.13_M14", "r96_P0.13_M13", "r96_P-0.13_M13", "r96_P0.13_M12", "r96_P-0.13_M12", "r96_P0.13_M11", "r96_P-0.13_M11", "r96_P0.13_M10", "r96_P-0.13_M10", "r96_P0.13_M9", "r96_P-0.13_M9", "r96_P0.13_M8", "r96_P-0.13_M8", "r96_P0.13", "r96_P-0.13" ]
if False :
  smooth = 0.157
  all_lin_cls = np.array([
    getCls("/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"+dir+"/dust_lensing/", 1000.0, smooth)
    for dir in lin_dirs])
  all_nl_cls = np.array([
    getCls("/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"+dir+"/dust_lensing/", 1.0, smooth)
    for dir in nl_dirs])
if False :
  smooth = 0.157
  all_lin_alms = np.array([
    getalms("/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"+dir+"/dust_lensing/", 1000.0, smooth)
    for dir in lin_dirs])
  all_nl_alms = np.array([
    getalms("/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"+dir+"/dust_lensing/", 1.0, smooth)
    for dir in nl_dirs])

if False :
  smooth = 0.157
  for M in range(1,10) :
    base_dir = "/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"
    lin_dir = base_dir + "r96_P0.00000013_M"+str(M)+"/dust_lensing/"
    nl_dir = base_dir + "r96_P0.13_M"+str(M)+"/dust_lensing/"
    loadData(lin_dir, 1000.0)
    loadData(nl_dir, 1.0)


lindir_p_96 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r96_P0.00000013_M1/dust_lensing.1/'
rho_data_lin_96 = loadData(lindir_p_96, 1000.0)
rho_interp_lin_96 = rho_interp(rho_data_lin_96)
rhos_lin_96 = rho(rho_interp_lin_96, 0.5)
cls_lin_96 = hp.anafast(rhos_lin_96)

nldir_p_96 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r96_P0.13_M1/dust_lensing.1/'
rho_data_nl_96 = loadData(nldir_p_96, 1.0)
rho_interp_nl_96 = rho_interp(rho_data_nl_96)
rhos_nl_96 = rho(rho_interp_nl_96, 0.5)
cls_nl_96 = hp.anafast(rhos_nl_96)

hpmin = np.min(rhos_lin_96)
hpmax = np.max(rhos_lin_96)
hp.mollview(rhos_lin_96, min=hpmin, max=hpmax)
plt.savefig('rhos_lin_96.png')
plt.close()
hp.mollview(rhos_nl_96, min=hpmin, max=hpmax)
plt.savefig('rhos_nl_96.png')
plt.close()

# print( ( (rhos_nl_96-np.mean(rhos_nl_96)) - (rhos_lin_96-np.mean(rhos_lin_96)) ) / np.std(rhos_lin_96) )
delta_lin_96 = (rhos_lin_96-np.mean(rhos_lin_96))/np.mean(rhos_lin_96)
delta_nl_96 = (rhos_nl_96-np.mean(rhos_nl_96))/np.mean(rhos_nl_96)


lindir_p_128 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r128_P0.00000013_M1/dust_lensing.1/'
rho_data_lin_128 = loadData(lindir_p_128, 1000.0)
rho_interp_lin_128 = rho_interp(rho_data_lin_128)
rhos_lin_128 = rho(rho_interp_lin_128, 0.5)
cls_lin_128 = hp.anafast(rhos_lin_128)

nldir_p_128 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r128_P0.13_M1/dust_lensing.1/'
rho_data_nl_128 = loadData(nldir_p_128, 1.0)
rho_interp_nl_128 = rho_interp(rho_data_nl_128)
rhos_nl_128 = rho(rho_interp_nl_128, 0.5)
cls_nl_128 = hp.anafast(rhos_nl_128)

hpmin = np.min(rhos_lin_128)
hpmax = np.max(rhos_lin_128)
hp.mollview(rhos_lin_128, min=hpmin, max=hpmax)
plt.savefig('rhos_lin_128.png')
plt.close()
hp.mollview(rhos_nl_128, min=hpmin, max=hpmax)
plt.savefig('rhos_nl_128.png')
plt.close()

delta_lin_128 = (rhos_lin_128-np.mean(rhos_lin_128))/np.mean(rhos_lin_128)
delta_nl_128 = (rhos_nl_128-np.mean(rhos_nl_128))/np.mean(rhos_nl_128)

# print( ( (rhos_nl_128-np.mean(rhos_nl_128)) - (rhos_lin_128-np.mean(rhos_lin_128)) ) / np.std(rhos_lin_128) )


print(np.std(delta_lin_96-delta_nl_96)/np.std(delta_lin_96))
print(np.std(delta_lin_128-delta_nl_128)/np.std(delta_lin_128))
