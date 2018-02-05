import numpy as np
import healpy as hp

nside = 16

pix_nums = np.arange(nside*nside*12)
vecs = np.transpose(np.array(hp.pix2vec(nside, pix_nums)))

np.savetxt("nside_" + str(nside) + ".vecs", vecs)
