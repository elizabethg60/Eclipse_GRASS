import os
import h5py
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
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
    return median_SbyN, np.nanmedian(fluxarray[xd_index,d_pix-window:d_pix+window+1])

Fiber = 'SCI'
directory = os.fsencode("Data")
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
airmass = [1.397197412189941, 1.3908964898437286, 1.3846972732593263, 1.3785979512128343, 1.3725967617534962, 1.3666922739847287, 1.3608825971962286, 1.3551663794081479, 1.3495420247000356, 1.3440080462174968, 1.338562994232855, 1.3332054550688581, 1.3279339868627922, 1.3227473724868755, 1.3176442361290022, 1.3126232984159738, 1.3076833112061232, 1.3028229984543085, 1.2980412321549541, 1.2933369657133118, 1.2887088955749932, 1.2841559202539607, 1.2796770184135349, 1.2752712437835534, 1.2709374083602234, 1.2666745635205159, 1.262481730743364, 1.2583579547261283, 1.2543022541463438, 1.2503138164793792, 1.2463917031360474, 1.242535045953556, 1.2387429972122834, 1.2350146844216574, 1.231349389438315, 1.227746234241975, 1.2242045776341632, 1.2207234955340522, 1.2173023802873375, 1.2139404716092423, 1.2106370685308412, 1.2073914865417654, 1.2042030189971746, 1.2010710903022697, 1.1979950242101485, 1.1949741984029105, 1.192008005256138, 1.1890958166194054, 1.1862371237435925, 1.1834313256509845, 1.180677870377596, 1.1779762191572636, 1.175325814378598, 1.1727262069264732, 1.1701768636113432, 1.1676772957447967, 1.1652270265355655, 1.1628255908551044, 1.1604725068404105, 1.1581673889309898, 1.1559097768988533, 1.153699223609928, 1.1515353983400873, 1.1494177967090575, 1.147346106919971, 1.1453199236398466, 1.1433388772813957, 1.141402607940455, 1.1395107200031651, 1.1376629860384224, 1.13585898323217, 1.134098390388702, 1.1323809567393825, 1.1307062940852546, 1.1290741674951181, 1.1274842684823778, 1.1259362989774035, 1.1244300430314391, 1.1229651645510026, 1.1215414452533292, 1.1201586707862459, 1.118816550728384, 1.1175148702417272, 1.1162533900963438, 1.115031924092679, 1.11385024742005, 1.1127081580911469, 1.1116054613557487, 1.1105419696107355, 1.1095174901736817, 1.1085318742923187, 1.1075849427372155, 1.10667652500582, 1.1058064901143831, 1.104974680613193, 1.1041809569397172, 1.1034251861412232, 1.1027072333648749, 1.1020270040832578, 1.1013843519450868, 1.1007791939203093, 1.1002114147029067, 1.099680939874652, 1.0991876614583995, 1.0987315042909727, 1.0983123879404313, 1.0979302531453083, 1.0975850313724296, 1.0972766620084753, 1.0970051015170965, 1.0967703010408376, 1.0965722217471665, 1.096410829130289, 1.096286101648452, 1.0961980125685193, 1.0961465510767265, 1.0961317072306893, 1.0961534783901423, 1.096211867929, 1.0963068866074763, 1.0964385475358078, 1.0966068730756326, 1.0968118960978699, 1.0970536339056294]
df = pd.DataFrame(columns= ['datetime','airmass'] + wavelength)

# file_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_April/model_data.jld2", "r")
# intensity_list = file_data["intensity_list"][()]

# intensity_array = [[]]*len(wavelength)
# for i in range(0,len(wavelength)):
#     intensity = (file_data[intensity_list[i]][()][24:,])
#     intensity_array[i].append(intensity/max(intensity))

ind = 0
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    row = []
    if filename.endswith(".fits"): 
        inputSpectrum = fits.open('Data/{}'.format(filename))
        row.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f")) 
        row.append(airmass[ind])
        for wavl_ind in range(0, len(wavelength)):
            row.append(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[1])
        df.loc[len(df.index)] = row
        ind += 1

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for i in wavelength:
#     ax1.plot(df['airmass'], df[i])
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("SNR") 
# plt.savefig("ext_coeff_eclipse.png")

fig = plt.figure()
ax1 = fig.add_subplot()
for i in range(0, len(wavelength)):
    ax1.plot(df['datetime'], df[wavelength[i]]/max(df[wavelength[i]]))
#ax1.set_xlabel("airmass")
ax1.set_xlabel("hour on 4/8")
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_ylabel("relative flux") 
plt.savefig("flux_coeff_eclipse.png")
plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for i in range(0, len(wavelength)):
#     ax1.plot(df['airmass'], intensity_array[0][i])
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative intensity") 
# plt.savefig("intensity_coeff_eclipse.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for i in range(0, len(wavelength)):
#     ax1.plot(df['airmass'], (df[wavelength[i]]/max(df[wavelength[i]]))/(intensity_array[0][i]))
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux/intensity") 
# plt.savefig("coeff_eclipse.png")
# plt.clf()

df_pyrrheliometer = pd.read_csv("Data/neid_ljpyrohelio_chv0_20240408.tel.txt", sep=" ", names=["Date", "Photocell", "Irradiance"])
df_pyrrheliometer["DateTime"] = [datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f") for date in df_pyrrheliometer['Date']]
df_pyrrheliometer["Day"] = [str(datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f").date()) for date in df_pyrrheliometer['Date']]
df_pyrrheliometer["Hour"] = pd.to_datetime(df_pyrrheliometer['DateTime']).dt.hour 
df= df_pyrrheliometer.loc[(df_pyrrheliometer['Day'] == "2024-04-08") & (df_pyrrheliometer['Hour'] >= 16)]

fig = plt.figure()
ax1 = fig.add_subplot()
ax1.plot(np.array(df["DateTime"]), np.array(df['Irradiance']))
ax1.set_ylabel("Irradiance (W/m^2)") 
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
ax1.set_xlabel("hour on 4/8")
plt.savefig("pyrrheliometer.png")
plt.clf()