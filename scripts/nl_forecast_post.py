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

def zofE(E) :
  return E/E[0] - 1

def map_interp(data, field='rho') :
  print("Loading E data...")
  E = data['E']
  print("Loading `"+field+"` data...")
  field_map = data[field]
  print("Interpolating...")
  interps = [interp1d(zofE(E[:,n]), field_map[:,n]) for n in range(field_map.shape[1])]
  print("Done.")
  return interps

def eval_interps(interp_fn, z) :
  return np.array([interp_fn[n](z) for n in range(len(interp_fn))])

nsz_data=np.load('ns_of_z.npz')
nsz = interp1d(nsz_data['zs'], nsz_data['ns'])
def bG(z) :
  return 0.95 + 0.67*z


def getmap(dir, mlt=1.0, smooth=0.0, z=0.25, field='rho') :
  if dir[0] != "/" :
    dir = "/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"+dir+"/dust_lensing/"
  print("Working on data in", dir)
  data = loadData(dir, mlt)
  interps = map_interp(data, field)
  field_map = eval_interps(interps, z)
  return field_map

def laplacian(mymap) :
  nside = hp.npix2nside(len(mymap))
  alms = hp.map2alm(mymap)
  lmax = hp.Alm.getlmax(len(alms))
  (l,m) = hp.Alm.getlm(lmax)
  alms = -alms/l/(l+1.0)
  alms[0] = 0
  return hp.alm2map(alms, nside)

def gradmapsquared(mymap) :
  nside = hp.npix2nside(len(mymap))
  lmax = hp.Alm.getlmax(len(mymap))
  (l,m) = hp.Alm.getlm(lmax)
  alms = hp.map2alm(mymap)
  mymap, dphi_map, dtheta_map = hp.alm2map_der1(alms, nside)
  return dphi_map*dphi_map + dtheta_map*dtheta_map


def getalms(dir, mlt=1.0, smooth=0.0, z=0.25) :
  rho_map = getmap(dir, mlt, smooth, z, 'rho')
  N_map = (rho_map/np.average(rho_map)-1.0)*nsz(z)*bG(z)
  if smooth:
    N_map = hp.sphtfunc.smoothing(N_map, smooth)
  alms = hp.map2alm(N_map)
  return alms

def getCls(dir, mlt=1.0, smooth=0.0, z=0.25) :
  alms = getalms(dir, mlt, smooth)
  return hp.alm2cl(alms)[:30]

alin_dirs_harm = [ "v1_r96_P0.00000013_M"+str(i)+"_resc1.0_Harmonic_nlsrc" for i in range(1,20) ] # All linear
rlin_dirs_harm = [ "v1_r96_P0.13_M"+str(i)+"_resc1.0e-3_Harmonic_rsc1.0e-3" for i in range(1,20) ] # metric & raysheet linear
mlin_dirs_harm = [ "v1_r96_P0.13_M"+str(i)+"_resc1.0e-3_Harmonic_nlsrc" for i in range(1,20) ] # metric (only) linear
nl_dirs_harm = [ "v1_r96_P0.13_M"+str(i)+"_resc1.0_Harmonic_nlsrc" for i in range(1,20) ] # full nonlinear

lin_dirs_cflrw = [ "r96_P0.13_M"+str(i)+"_resc1.0e-3_ConformalFLRW" for i in range(1,20) ]
nl_dirs_cflrw = [ "r96_P0.13_M"+str(i)+"_resc1.0_ConformalFLRW" for i in range(1,20) ]

if True :
  smooth = 0.0 #0.157

  # all_alin_rho = np.array([ getmap(dir, mlt=1000.0, field='rho') for dir in alin_dirs_harm])
  # all_mlin_rho = np.array([ getmap(dir, field='rho') for dir in mlin_dirs_harm])
  # all_rlin_rho = np.array([ getmap(dir, field='rho') for dir in rlin_dirs_harm])
  # all_nl_rho = np.array([ getmap(dir, field='rho') for dir in nl_dirs_harm])

  # all_alin_DA = np.array([ getmap(dir, mlt=1000.0, field='DA') for dir in alin_dirs_harm])
  # all_mlin_DA = np.array([ getmap(dir, field='DA') for dir in mlin_dirs_harm])
  # all_nl_DA = np.array([ getmap(dir, field='DA') for dir in nl_dirs_harm])

  np.savetxt('all_alin_rho.txt', all_alin_rho)
  np.savetxt('all_rlin_rho.txt', all_rlin_rho)
  np.savetxt('all_mlin_rho.txt', all_mlin_rho)
  np.savetxt('all_nl_rho.txt', all_nl_rho)

  # np.savetxt('all_alin_d2rho.txt', np.array([ laplacian(alin) for alin in all_alin_rho ]) )
  # np.savetxt('all_mlin_d2rho.txt', np.array([ laplacian(alin) for alin in all_mlin_rho ]) )
  # np.savetxt('all_nl_d2rho.txt', np.array([ laplacian(alin) for alin in all_nl_rho ]) )

  # np.savetxt('all_alin_drho2.txt', np.array([ laplacian(alin) for alin in all_alin_rho ]) )
  # np.savetxt('all_mlin_drho2.txt', np.array([ gradmapsquared(alin) for alin in all_mlin_rho ]) )
  # np.savetxt('all_nl_drho2.txt', np.array([ gradmapsquared(alin) for alin in all_nl_rho ]) )

  # np.savetxt('all_alin_DA.txt', all_alin_DA)
  # np.savetxt('all_mlin_DA.txt', all_mlin_DA)
  # np.savetxt('all_nl_DA.txt', all_nl_DA)
  # np.savetxt('all_alin_d2DA.txt', np.array([ laplacian(alin) for alin in all_alin_DA ]) )
  # np.savetxt('all_mlin_d2DA.txt', np.array([ laplacian(alin) for alin in all_mlin_DA ]) )
  # np.savetxt('all_nl_d2DA.txt', np.array([ laplacian(alin) for alin in all_nl_DA ]) )
  # np.savetxt('all_alin_d2DA.txt', np.array([ laplacian(alin) for alin in all_alin_DA ]) )
  # np.savetxt('all_mlin_dDA2.txt', np.array([ gradmapsquared(alin) for alin in all_mlin_DA ]) )
  # np.savetxt('all_nl_dDA2.txt', np.array([ gradmapsquared(alin) for alin in all_nl_DA ]) )

  # all_flin_harm = np.array([ getmap(dir) for dir in flin_dirs_harm])
  # np.savetxt('all_flin_harm.txt', all_flin_harm)
  # all_lin_cls_cflrw = np.array([ getCls(dir) for dir in lin_dirs_cflrw])
  # np.savetxt('all_lin_cls_cflrw.txt', all_lin_cls_cflrw)
  # all_nl_cls_cflrw = np.array([ getCls(dir) for dir in nl_dirs_cflrw])
  # np.savetxt('all_nl_cls_cflrw.txt', all_nl_cls_cflrw)

if False :
  alin_map = getmap("v1_"+alin_dirs_harm[0], z=0.250, mlt=1000.0) # All linearized
  # flin_map = getmap("v1_"+flin_dirs_harm[0], z=0.50) # Fluid + metric linearized
  mlin_map = getmap("v1_"+mlin_dirs_harm[0], z=0.250) # linearized metric only
  nl_map = getmap("v1_"+nl_dirs_harm[0], z=0.250)
  print(np.mean(np.abs((nl_map - alin_map))/np.std(nl_map)))
  # print(np.mean(np.abs((nl_map - flin_map))/np.std(nl_map)))
  print(np.mean(np.abs((nl_map - mlin_map))/np.std(nl_map)))

  # hpmin = np.min((nl_map - flin_map)/np.std(nl_map))
  # hpmax = np.max((nl_map - flin_map)/np.std(nl_map))
  # hp.mollview(np.abs((nl_map - alin_map))/np.std(alin_map), min=hpmin, max=hpmax)
  # plt.savefig('alin_map.png')
  # plt.close()
  # hp.mollview(np.abs((nl_map - flin_map))/np.std(flin_map), min=hpmin, max=hpmax)
  # plt.savefig('flin_map.png')
  # plt.close()
  # hp.mollview(np.abs((nl_map - mlin_map))/np.std(mlin_map), min=hpmin, max=hpmax)
  # plt.savefig('mlin_map.png')
  # plt.close()

if False :
  map128 = getmap("r128_P0.13_M1_resc1.0_Harmonic")
  map96 = getmap("r96_P0.13_M1_resc1.0_Harmonic")
  print(np.mean(np.abs((map128 - map96))/np.std(map128)))

  map128 = getmap("r128_P0.13_M1_resc1.0_Harmonic_nlsrc")
  map96 = getmap("r96_P0.13_M1_resc1.0_Harmonic_nlsrc")
  print(np.mean(np.abs((map128 - map96))/np.std(map128)))

if False :
  smooth = 0.157
  for M in range(1,10) :
    base_dir = "/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/"
    lin_dir = base_dir + "r96_P0.00000013_M"+str(M)+"/dust_lensing/"
    nl_dir = base_dir + "r96_P0.13_M"+str(M)+"/dust_lensing/"
    loadData(lin_dir, 1000.0)
    loadData(nl_dir, 1.0)


if False:
  lindir_p_96 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r96_P0.13_M1_resc1.0e-3/dust_lensing.1/'
  rho_data_lin_96 = loadData(lindir_p_96, 1.0)
  rho_interp_lin_96 = rho_interp(rho_data_lin_96)
  rhos_lin_96 = rho(rho_interp_lin_96, 0.5)
  cls_lin_96 = hp.anafast(rhos_lin_96)

  nldir_p_96 = '/home/jbm120/nonlinear_forecast/cosmograph/build/nl_runs/r96_P0.13_M1_resc1.0/dust_lensing/'
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


