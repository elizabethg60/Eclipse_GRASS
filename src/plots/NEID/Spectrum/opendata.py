import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

specdata = fits.open('Data/neidL2_20231014T223048.fits')
# specdata.info() # shows file extensions/info
# print(specdata[0].header)
# print(specdata[0].header["EXTSNR"])

sciflux=specdata[1].data#[]
scilambda=specdata[7].data#[]
# print(sci.shape)

# plt.plot(np.concatenate(scilambda),np.concatenate(sciflux),label="original")
# plt.legend()
# plt.show()

plt.figure()
plt.imshow(scilambda, aspect='auto', interpolation='nearest')
plt.colorbar()
plt.savefig('wavelength.png')

plt.figure()
plt.imshow(sciflux, aspect='auto', vmin=-10, vmax=10)
plt.colorbar()
plt.savefig('flux.png')