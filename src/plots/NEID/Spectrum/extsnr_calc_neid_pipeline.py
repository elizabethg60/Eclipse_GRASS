import numpy  as np
from astropy.io import fits

def calc_snr_at_wav(inputSpectrum,wavl,Fiber,calculate_at_center=True):
    nd = inputSpectrum[0].header
    
    if wavl is not None:
        Wavl = wavl
    else:
        if isinstance(nd['QSNRWL'],(float,int)) and (3500 < nd['QSNRWL'] < 12000):
            Wavl = nd['QSNRWL']
        else:
            print('No sensible wavelength value in QSNRWL={0}. Defaulting to 5529.7 A'.format(nd.pri0['QSNRWL']))
            Wavl = 5529.7
    
    fluxarray = inputSpectrum[1].data
    wavearray = inputSpectrum[7].data
    vararray = inputSpectrum[4].data
    
    # Find the index of the wavl array which contains the args.Wavl closest to the center of the array
    try:
        disp_pixels = np.argmin(np.abs(wavearray - Wavl),axis=1)
    except TypeError:
        print('{} No wavelength solution. Skipping SNR calculation.'.format(inputSpectrum))
        return
    xd_index = np.argmin(np.abs(disp_pixels - wavearray.shape[1]/2))
    d_pix = int(wavearray.shape[1]/2) if calculate_at_center else disp_pixels[xd_index]
    
    print('Pixel coordinate for calculating SNR of {0} is at ({1},{2})'.format(Wavl,xd_index,d_pix))
    # Calculate the median S/N inside a window of 50 pixels around the wavelength in spectrum
    window = 25  #pixels
    median_SbyN = np.nanmedian(fluxarray[xd_index,d_pix-window:d_pix+window+1]/np.sqrt(vararray[xd_index,d_pix-window:d_pix+window+1]))
    print('S/N of {0} at {1} A = {2}'.format(inputSpectrum,Wavl,median_SbyN))
    return median_SbyN
    
inputSpectrum = fits.open('Data/Level1/neidL1_20231014T151516.fits')
wavl = 5529.7
Fiber = 'SCI'

print(calc_snr_at_wav(inputSpectrum,wavl,Fiber,calculate_at_center=True))