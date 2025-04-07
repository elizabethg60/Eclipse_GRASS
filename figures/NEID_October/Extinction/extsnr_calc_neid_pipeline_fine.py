import os
import h5py
import csv
import numpy as np
import pandas as pd
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def float_range(start, end, num_steps):
    step_size = (end - start) / (num_steps - 1)
    current = start
    while current <= end:
        yield current
        current += step_size

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
    # Calculate the median S/N inside a window of 120 pixels around the wavelength in spectrum
    window = 60  #pixels
    median_SbyN = np.nanmedian(fluxarray[xd_index,d_pix-window:d_pix+window+1]/np.sqrt(vararray[xd_index,d_pix-window:d_pix+window+1]))
    #print('S/N of {0} at {1} A = {2}'.format(inputSpectrum,Wavl,median_SbyN))
    return median_SbyN, np.nanmedian(fluxarray[xd_index,d_pix-window:d_pix+window+1])

def NEID_snr(Fiber, directory, airmass, wavelength, file):
    df = pd.DataFrame(columns= ['datetime','airmass'] + wavelength)
    ind = 0
    for j in range(0,len(timestamps_october)):
        row = []
        inputSpectrum = fits.open('{}/{}'.format(directory, np.array(file)[j]))
        row.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f")) 
        row.append(airmass[ind])
        for wavl_ind in range(0, len(wavelength)):
            row.append(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[0]) 
        df.loc[len(df.index)] = row
        ind += 1
    return df

def NEID_flux(Fiber, directory, airmass, wavelength, file):
    df = pd.DataFrame(columns= ['datetime','airmass'] + wavelength)
    ind = 0
    for j in range(0,len(timestamps_october)):
        row = []
        inputSpectrum = fits.open('{}/{}'.format(directory, np.array(file)[j]))
        row.append(datetime.strptime(inputSpectrum[0].header["DATE-OBS"], "%Y-%m-%dT%H:%M:%S.%f")) 
        row.append(airmass[ind])
        for wavl_ind in range(0, len(wavelength)):
            row.append(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[1])
        df.loc[len(df.index)] = row
        ind += 1
    return df

def model_flux(file_data):
    intensity_list = file_data["intensity_list"][()]
    intensity_array = [[]]*len(wavelength)
    for i in range(0,len(wavelength)):
        intensity = (file_data[intensity_list[i]][()])
        intensity_array[i].append(intensity/max(intensity))
    return intensity_list, intensity_array

Fiber = 'SCI'
path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")["filename"]
timestamps_october = timestamps_full_october[15:-150]
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
airmass = [2.5632073517047638, 2.5361097748885557, 2.5093661650400723, 2.4832872029346307, 2.458152719768879, 2.4333280353218156, 2.4091026139254765, 2.3857380385003344, 2.362645575113234, 2.3400953849743926, 2.3183323449274487, 2.2968091303568645, 2.27577829301206, 2.2554693819371017, 2.2353725338593415, 2.2157241944955177, 2.1967398191699803, 2.1779435542209677, 2.159557077563393, 2.1415684495759097, 2.1241760200750517, 2.106944668512874, 2.0900782010965315, 2.0737632866181173, 2.05759228654005, 2.0417568347767414, 2.0264327849819934, 2.011237707679954, 1.9963519894425743, 1.9819414507550772, 1.9676468497988242, 1.9536381582317257, 1.9400718707467606, 1.926610129678612, 1.9134132425208898, 1.9006290194999498, 1.8879393729586693, 1.8754956248574794, 1.8634374833021545, 1.8514651818936574, 1.8397216717571738, 1.8282016018286706, 1.817034683533013, 1.8059435946574782, 1.7950609742097614, 1.7845095958663924, 1.7740276105655446, 1.763740511098882, 1.753764560828184, 1.743852354616247, 1.7341226871090905, 1.7246857022512947, 1.715307505275213, 1.7061006017579452, 1.697169374222985, 1.688292589769291, 1.6795768389575176, 1.6711210712224671, 1.6627159428484009, 1.6544624719613135, 1.646454476321036, 1.638493796560436, 1.630676192498115, 1.6230906273882857, 1.6155494821322123, 1.6081435466952902, 1.6008703361742513, 1.5938127284826231, 1.5867962627606653, 1.5799055052643187, 1.573219043936084, 1.5665716579824998, 1.5600435303724371, 1.5537091700100665, 1.5474121144362678, 1.5412283821572987, 1.5352285977595637, 1.5292646143956121, 1.523408487955313, 1.5177271385876636, 1.5120803277757258, 1.5065363351237604, 1.5011585426132867, 1.4958142450515572, 1.4905681178810164, 1.4854801586717268, 1.4804248500221022, 1.4754634222942877, 1.4705945147488082, 1.4658738276317298, 1.4611849448608079, 1.456584725994602, 1.4521258025734258, 1.447698240609205, 1.4433557809843884, 1.43914810522267, 1.43497149743696, 1.4308767024770255, 1.4269105616726943, 1.4229753336016377, 1.419118880539379, 1.4153853052103214, 1.4116826158282119, 1.408055897033766, 1.4045466068708823, 1.4010682957774412, 1.3976633677805828, 1.3943707230207358, 1.3911092629152901, 1.387918800523127, 1.3847985802902372, 1.381784212071226, 1.3788014670576347, 1.3758868280098038, 1.3730735241834369, 1.3702922880917996, 1.3675771947689528, 1.364959158943391, 1.3623737296001717, 1.3598526423845454, 1.3574245581270798, 1.3550297095145418, 1.3526975555563139, 1.350454558686985, 1.3482455138644898, 1.346097660961062, 1.3440353146548754, 1.34200772145745, 1.3400399546266697, 1.3381542271029911, 1.3363041363785326, 1.334512636897446, 1.3327998814765274, 1.331123727652466, 1.329505054361501, 1.3279435167148383, 1.3264565754953044, 1.3250076481928261, 1.323614903675427, 1.3222938236556199, 1.3210119084301528, 1.319785336067303, 1.3186276374092976, 1.3175103337752088, 1.3164476420975404, 1.3154511668285525, 1.3144963965600829, 1.3135956129740973, 1.3127585150081922, 1.3119645125498263, 1.3112239742717597, 1.310544710948751, 1.3099100153508805, 1.3093283618247469, 1.3088056868926414, 1.3083291351789732, 1.307905301914301, 1.3075341005197691]
# Normalize the values to the range [0, 1]
norm = mcolors.Normalize(vmin=np.min(wavelength), vmax=np.max(wavelength))
# Create a colormap from blue to red
cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

def NEID_no_ext_SSD():
    df = NEID_flux(Fiber, path_october, airmass, wavelength, timestamps_october)[0:-28] 

    grass_data_no_cb_1 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_1.jld2", "r")
    brightness_1 = grass_data_no_cb_1["brightness"][()]
    grass_data_no_cb_2 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_2.jld2", "r")
    brightness_2 = grass_data_no_cb_2["brightness"][()]

    fig = plt.figure()
    ax1 = fig.add_subplot()
    for i, value in enumerate(wavelength[0:12]):
        original = grass_data_no_cb_1[brightness_1[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        airmass_arr = df['airmass']
        comp_arr = ((df[wavelength[i]]/max(df[wavelength[i]]))/(array/max(array)))
        ax1.plot(airmass_arr, comp_arr, color=color)

    for i, value in enumerate(wavelength[12:23]):
        original = grass_data_no_cb_2[brightness_2[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        airmass_arr = df['airmass']
        comp_arr = ((df[wavelength[i]]/max(df[wavelength[i]]))/(array/max(array)))
        ax1.plot(airmass_arr, comp_arr, color=color)
    ax1.set_xlabel("airmass")
    ax1.set_ylabel("relative flux/intensity") 
    ax1.axvline(x = df['airmass'][80])
    plt.savefig("figures/no_ext/fine_time/coeff_eclipse_SSD_fine.png")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot()  
    for i, value in enumerate(wavelength[0:12]):
        original = grass_data_no_cb_1[brightness_1[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], array/max(array), label = wavelength[i], color=color, s = 1)
        ax1.plot(df['airmass'], df[wavelength[i]]/max(df[wavelength[i]]), color=color)

    for i, value in enumerate(wavelength[12:23]):
        original = grass_data_no_cb_2[brightness_2[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], array/max(array), label = wavelength[i], color=color, s = 1)
        ax1.plot(df['airmass'], df[wavelength[i]]/max(df[wavelength[i]]), color=color)
    rms = round(np.sqrt((np.nansum(((df[wavelength[8]]/max(df[wavelength[8]])) - array/max(array))**2))/len((df[wavelength[8]]/max(df[wavelength[8]])) - array/max(array))),4)
    ax1.text(np.array(df['airmass'])[-10], 0.5, "{} RMS {}".format(wavelength[8], rms))    
    ax1.set_xlabel("airmass")
    ax1.set_ylabel("relative flux") 
    plt.savefig("figures/no_ext/fine_time/flux_comp_SSD_fine.png")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot()  
    for i, value in enumerate(wavelength[0:12]):
        original = grass_data_no_cb_1[brightness_1[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], (df[wavelength[i]]/max(df[wavelength[i]])) - array/max(array), label = wavelength[i], color=color, s = 1)
    for i, value in enumerate(wavelength[12:23]):
        original = grass_data_no_cb_2[brightness_2[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], (df[wavelength[i]]/max(df[wavelength[i]])) - array/max(array), label = wavelength[i], color=color, s = 1)

    ax1.set_xlabel("airmass")
    ax1.set_ylabel("relative flux") 
    ax1.axvline(x = df['airmass'][80])
    plt.savefig("figures/no_ext/fine_time/flux_comp_diff_SSD_fine.png")
    plt.clf()

def determine_ext_SSD():
    df = NEID_flux(Fiber, path_october, airmass, wavelength, timestamps_october)[0:-28] 

    grass_data_no_cb_1 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_1.jld2", "r")
    brightness_1 = grass_data_no_cb_1["brightness"][()]
    grass_data_no_cb_2 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_2.jld2", "r")
    brightness_2 = grass_data_no_cb_2["brightness"][()]

    min_ext_arr = []

    for i, value in enumerate(wavelength[0:12]):
        color = cmap(norm(value))
        data_no_clouds = df[wavelength[i]]
        data_arr = data_no_clouds / np.max(data_no_clouds) 

        original = grass_data_no_cb_1[brightness_1[0]][()]
        model_intensity = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]

        min_chi2_ext = 10
        min_ext = 0
        for ext in float_range(0.09, 0.55, 50):
                I_ext = model_intensity * np.exp(-np.array(df["airmass"])*ext)
                ext_model = I_ext / np.max(I_ext)

                chi2_iterate = np.sum((data_arr - ext_model) ** 2 / ext_model)
                if chi2_iterate < min_chi2_ext:
                    min_chi2_ext = chi2_iterate
                    min_ext = ext
        min_ext_arr.append(min_ext)

        plt.scatter(df["airmass"], model_intensity * np.exp(-np.array(df["airmass"])*min_ext_arr[i])/ np.max(model_intensity * np.exp(-np.array(df["airmass"])*min_ext_arr[i])), color = color)
        plt.plot(df["airmass"], data_arr, color = color)
    plt.savefig("test")
    plt.clf()

    for i, value in enumerate(wavelength[12:23]):
        color = cmap(norm(value))
        data_no_clouds = df[wavelength[i]]
        data_arr = data_no_clouds / np.max(data_no_clouds) 

        original = grass_data_no_cb_2[brightness_2[0]][()]
        model_intensity = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]

        min_chi2_ext = 10
        min_ext = 0
        for ext in float_range(0.09, 0.55, 50):
                I_ext = model_intensity * np.exp(-np.array(df["airmass"])*ext)
                ext_model = I_ext / np.max(I_ext)

                chi2_iterate = np.sum((data_arr - ext_model) ** 2 / ext_model)
                if chi2_iterate < min_chi2_ext:
                    min_chi2_ext = chi2_iterate
                    min_ext = ext
        min_ext_arr.append(min_ext)

    # print(min_ext_arr) #0.15571428571428572

    # Create a DataFrame
    df = pd.DataFrame({
        'Wavelength': wavelength,
        'ext': min_ext_arr
    })
    # Save to CSV file
    df.to_csv('NEID_ext_coeff_SSD_fine_time.csv', index=False)

def NEID_ext_SSD():
    df = NEID_flux(Fiber, path_october, airmass, wavelength, timestamps_october)[0:-28] 

    grass_data_no_cb_1 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD_3ext/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_ext_1.jld2", "r")
    brightness_1 = grass_data_no_cb_1["brightness"][()]
    grass_data_no_cb_2 = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD_3ext/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_ext_2.jld2", "r")
    brightness_2 = grass_data_no_cb_2["brightness"][()]

    fig = plt.figure()
    ax1 = fig.add_subplot()
    for i, value in enumerate(wavelength[0:12]):
        original = grass_data_no_cb_1[brightness_1[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        airmass_arr = df['airmass']
        comp_arr = ((df[wavelength[i]]/max(df[wavelength[i]]))/(array/max(array)))
        ax1.plot(airmass_arr, comp_arr, color=color)

    for i, value in enumerate(wavelength[12:23]):
        original = grass_data_no_cb_2[brightness_2[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        airmass_arr = df['airmass']
        comp_arr = ((df[wavelength[i]]/max(df[wavelength[i]]))/(array/max(array)))
        ax1.plot(airmass_arr, comp_arr, color=color)
    ax1.set_xlabel("airmass")
    ax1.set_ylabel("relative flux/intensity") 
    ax1.axvline(x = df['airmass'][80])
    plt.savefig("figures/ext/fine_time/coeff_eclipse_SSD_fine.png")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot()  
    for i, value in enumerate(wavelength[0:12]):
        original = grass_data_no_cb_1[brightness_1[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], array/max(array), label = wavelength[i], color=color, s = 1)
        ax1.plot(df['airmass'], df[wavelength[i]]/max(df[wavelength[i]]), color=color)

    for i, value in enumerate(wavelength[12:23]):
        original = grass_data_no_cb_2[brightness_2[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], array/max(array), label = wavelength[i], color=color, s = 1)
        ax1.plot(df['airmass'], df[wavelength[i]]/max(df[wavelength[i]]), color=color)
    rms = round(np.sqrt((np.nansum(((df[wavelength[8]]/max(df[wavelength[8]])) - array/max(array))**2))/len((df[wavelength[8]]/max(df[wavelength[8]])) - array/max(array))),4)
    ax1.text(np.array(df['airmass'])[-10], 0.5, "{} RMS {}".format(wavelength[8], rms))    
    ax1.set_xlabel("airmass")
    ax1.set_ylabel("relative flux") 
    plt.savefig("figures/ext/fine_time/flux_comp_SSD_fine.png")
    plt.clf()

    fig = plt.figure()
    ax1 = fig.add_subplot()  
    for i, value in enumerate(wavelength[0:12]):
        original = grass_data_no_cb_1[brightness_1[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], (df[wavelength[i]]/max(df[wavelength[i]])) - array/max(array), label = wavelength[i], color=color, s = 1)
    for i, value in enumerate(wavelength[12:23]):
        original = grass_data_no_cb_2[brightness_2[0]][()]
        array = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        color = cmap(norm(value))
        ax1.scatter(df['airmass'], (df[wavelength[i]]/max(df[wavelength[i]])) - array/max(array), label = wavelength[i], color=color, s = 1)

    ax1.set_xlabel("airmass")
    ax1.set_ylabel("relative flux") 
    ax1.axvline(x = df['airmass'][80])
    plt.savefig("figures/ext/fine_time/flux_comp_diff_SSD_fine.png")
    plt.clf()  

def flux_ext_no_comp():
    df = NEID_flux(Fiber, path_october, airmass, wavelength, timestamps_october)[0:-28] 
    grass_data_no_cb_no_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_1.jld2", "r")
    brightness_no_ext = grass_data_no_cb_no_ext["brightness"][()]
    grass_data_no_cb_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/SSD/Fine_Time/SSD_3ext/data/neid_all_lines_rv_off_SSD_gpu_fine_spectra_ext_1.jld2", "r")
    brightness_ext = grass_data_no_cb_ext["brightness"][()]

    fig, axs = plt.subplots(1, 2, sharex=False, sharey=True, figsize=(10, 5))
    for i, value in enumerate([wavelength[2]]):
        print([wavelength[2]])
        color = cmap(norm(value))
        data_no_clouds = df[wavelength[i]]
        data = data_no_clouds / np.max(data_no_clouds) #np.mean(data_no_clouds[115:len(data_no_clouds)])
        axs[0].plot(df['airmass'], data, color=color, label = "NEID")
        original = grass_data_no_cb_no_ext[brightness_no_ext[0]][()]
        no_ext = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        original = grass_data_no_cb_ext[brightness_ext[0]][()]
        ext = [np.mean(original[j:j+56]) for j in range(0, len(original), 56)]
        axs[0].scatter(df['airmass'], no_ext/max(no_ext), color=color, s = 1, label = "Model")
        axs[1].plot(df['airmass'], data, color=color)
        axs[1].scatter(df['airmass'], ext/max(ext), label = wavelength[i], color=color, s = 1)
    axs[0].set_xlabel("Airmass")
    axs[0].set_ylabel("Relative Intenstiy")
    axs[1].set_xlabel("Airmass")
    axs[0].set_title('SSD Limb Darkening Only')
    axs[1].set_title('SSD Limb Darkening + Extinction')
    axs[0].legend()
    plt.savefig("2nd_yr_talk.png")

# NEID_no_ext_SSD()
# determine_ext_SSD()
# NEID_ext_SSD()
flux_ext_no_comp()