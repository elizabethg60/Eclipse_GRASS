import os
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

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
    
    #print('Pixel coordinate for calculating SNR of {0} is at ({1},{2})'.format(Wavl,xd_index,d_pix))
    # Calculate the median S/N inside a window of 50 pixels around the wavelength in spectrum
    window = 25  #pixels
    median_SbyN = np.nanmedian(fluxarray[xd_index,d_pix-window:d_pix+window+1]/np.sqrt(vararray[xd_index,d_pix-window:d_pix+window+1]))
    #print('S/N of {0} at {1} A = {2}'.format(inputSpectrum,Wavl,median_SbyN))
    return median_SbyN

Fiber = 'SCI'
directory = os.fsencode("Data")
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
df = pd.DataFrame(columns= ['datetime','airmass']+ wavelength)

ind = 0
airmass = [2.572425367171159, 2.5448735167307444, 2.5180129852674447, 2.4918196512704975, 2.4662705139035777, 2.4413436288441934, 2.4170177578145795, 2.3932734823889623, 2.37009138693885, 2.3474531938592382, 2.3253414202609437, 2.3037393354077356, 2.2826306684738125, 2.2620005864480275, 2.2418338863836254, 2.222117195250764, 2.202836269389674, 2.1839783769544594, 2.165530844281854, 2.147482158803499, 2.1298206065409873, 2.1125351465601323, 2.0956149615094435, 2.079050251585693, 2.062830987750628, 2.046947716592267, 2.0313913365785403, 2.016153081856585, 2.0012245071345185, 1.986597298269381, 1.972263962461992, 1.958216754226753, 1.9444483742854706, 1.9309517792954787, 1.9177200124152, 1.9047468298365413, 1.8920257307863801, 1.8795505914120514, 1.8673154942382788, 1.8553147196941986, 1.8435427380566671, 1.831994063434159, 1.820663666748346, 1.809546809503011, 1.798638241102934, 1.7879334089554906, 1.777427647870953, 1.7671168159344612, 1.7569965267146577, 1.7470626546708954, 1.7373112019147001, 1.727738293459651, 1.718340060053852, 1.7091130864141855, 1.7000537251133014, 1.691158549369504, 1.6824242345351446, 1.6738474516746664, 1.665425277118639, 1.6571545668250862, 1.6490323699354348, 1.6410558210923736, 1.6332221375080378, 1.6255285239477424, 1.6179725405672654, 1.6105515420578271, 1.6032630490559034, 1.5961046518781352, 1.5890740082488395, 1.5821688411184676, 1.5753868552864425, 1.5687260619706278, 1.5621842847995455, 1.5557594874908791, 1.5494496136755136, 1.5432528890610906, 1.5371673619055677, 1.5311912080589878, 1.5253226524609174, 1.5195599676196931, 1.5139014722332405, 1.5083454632407058, 1.5028904820389442, 1.49753491009505, 1.4922772377184867, 1.487115995081263, 1.4820497510829715, 1.4770771122931756, 1.4721966633832764, 1.4674072011461599, 1.4627073793625198, 1.4580959450772295, 1.4535716239098995, 1.449133337013963, 1.444779871923669, 1.4405100505253015, 1.436322928737646, 1.432217235163834, 1.428192079430581, 1.4242463978391648, 1.4203792034339793, 1.4165895351965723, 1.412876457412079, 1.4092390590575627, 1.4056764105696757, 1.402187734730558, 1.398772147589441, 1.3954288311855378, 1.3921569503535707, 1.388955847927697, 1.3858246154153044, 1.3827626345546815, 1.3797691539915007, 1.376843480268806, 1.373984904963574, 1.3711928417826138, 1.368466618713325, 1.3658056167673367, 1.3632092346354667, 1.3606768883388554, 1.3582080108923182, 1.3558020519795388, 1.3534584496721638, 1.3511767427430967, 1.348956373841781, 1.3467969357634277, 1.3446978525698496, 1.342658745016126, 1.3406791455438265, 1.3387586268041718, 1.3368967756697012, 1.3350931930210945, 1.3333474727522334, 1.331659285322119, 1.3300282508874752, 1.32845402425539, 1.3269362732391417, 1.325474678483048, 1.324068916628486, 1.3227187274887562, 1.321423811884163, 1.320183900178243, 1.3189987347982757, 1.3178680701008785, 1.3167916722521957, 1.3157693070586747, 1.3148007886685817, 1.313885905273252, 1.3130244691189576, 1.312216303719761, 1.311461234956191, 1.3107591268798375, 1.3101098268427607, 1.309513202596138, 1.3089691327704955, 1.3084775011853784, 1.3080382199450036, 1.3076511937122919]
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    row = []
    if filename.endswith(".fits"): 
        inputSpectrum = fits.open('Data/{}'.format(filename))
        # inputSpectrum.info() # shows file extensions/info
        # print(inputSpectrum[0].header)
        # sciflux=inputSpectrum[1].data#[]
        # scilambda=inputSpectrum[7].data#[]
        # print(sci.shape)
        # plt.plot(np.concatenate(scilambda),np.concatenate(sciflux),label="original")
        # plt.legend()
        # plt.show()
        # plt.figure()
        # plt.imshow(scilambda, aspect='auto', interpolation='nearest')
        # plt.colorbar()
        # plt.savefig('wavelength.png')

        # plt.figure()
        # plt.imshow(sciflux, aspect='auto', vmin=-10, vmax=10)
        # plt.colorbar()
        # plt.savefig('flux.png')

        row.append(inputSpectrum[0].header["DATE-OBS"])
        row.append(airmass[ind])
        for wavl in wavelength:
            row.append(calc_snr_at_wav(inputSpectrum,wavl,Fiber,calculate_at_center=True))
        df.loc[len(df.index)] = row
        ind += 1

fig = plt.figure()
ax1 = fig.add_subplot()
for i in wavelength:
    ax1.plot(df['airmass'], df[i], label = i)
ax1.set_xlabel("airmass")
ax1.set_ylabel("SNR") 
plt.legend()
plt.savefig("ext_coeff_eclipse.png")