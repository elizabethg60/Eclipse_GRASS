#add rms of all lines and get that rms (1 m/s)
#variability on vs off in grass
#check line by line rv 
#quantify photon noise of lines 

#weighted rms by rms or photon noise (error on line by line ccf) if it drops then limit is astrophysical or model
#look into line shape for blends / depth like eric said 
#eric's paper complex extinction model with linear weight to z function (Eric’s extinction + compare coef values reasonably match)
#Do for line by line rv - remove eclipse with model then do binning for entire eclipse time + do for multiple days 

#Other:
# do boulder october no ext
# Eric's intensity sensitive
# Eric's coverage area weighting
# same intensity pipeline on a different day

#Future: Boulder pattern sensitivity + ?  projected rv with NO cb wins out GRASS for some lines  - follow up ?
#future: calculate gravitational redshift of sun and a white dwarf + sensitivity test of LD + extinction coefficient
    # time weighting: 1. rerun (just my model for a single line) with (1) weighted flux midpoint and (2) 5 sec cadence - future analysis done with the best one
                    # 2. run with BEST timing option (N = 50 and subset = 10) for all lines (my model + GRASS) overnight (if not regular midpoint)
    # full resolution: 1. submit jobs to run updated code at full resolution for a single line (by then should have decided if 197 ok if 300 better)
                     # 2. submit jobs to run updated code at N = 197 and subgrid = 80 for a single line
    # extinction: try slope method again but with two different laws where extinction clearly changes slope
#plan b: consider how to validate model is correct in orientation / time, consider Reiners newer CB model, consider NEID SPF not a gaussian 
#barycentric correction is flux weighted by the counts-per-second recorded using exposure meter

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
            row.append(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[0]) #SNR
            # print("flux: {}".format(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[1]))
            # print("snr: {}".format(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[0]))
            # print("sqrt flux = snr : {}".format(np.sqrt(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[1])))
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
            row.append(calc_snr_at_wav(inputSpectrum,wavelength[wavl_ind],Fiber,calculate_at_center=True)[1]) #flux
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

def model_flux_3_ext(file_data):
    norm1 = pd.read_csv("/storage/home/efg5335/work/GRASS/data/NEID_three_extinction.csv")["ext1_norm"]
    norm2 = pd.read_csv("/storage/home/efg5335/work/GRASS/data/NEID_three_extinction.csv")["ext2_norm"]
    
    intensity_list = file_data["intensity_list"][()]
    intensity_array = [[]]*len(wavelength)
    for i in range(0,len(wavelength)):
        intensity = (file_data[intensity_list[i]][()])[0:-20]

        m1 = intensity[0:26]
        m2 = intensity[26:47]
        m3 = intensity[47:len(intensity)]
        
        m1_end = file_data_no_ext[intensity_list_no_ext[i]][()][0:-20][47:len(intensity)] * np.exp(-np.array(airmass[0:-20][47:len(intensity)])*norm1[i])
        m2_end = file_data_no_ext[intensity_list_no_ext[i]][()][0:-20][47:len(intensity)] * np.exp(-np.array(airmass[0:-20][47:len(intensity)])*norm2[i])
        
        m1 = list(np.array(m1) / np.mean(m1_end[68:len(intensity)]))
        m2 = list(np.array(m2) / np.mean(m2_end[68:len(intensity)]))
        m3 = list(np.array(m3) / np.mean(m3[68:len(intensity)]))

        intensity_array[i].append(np.array(m1+m2+m3))
    return intensity_list, intensity_array

Fiber = 'SCI'
path_october = "/storage/group/ebf11/default/pipeline/neid_solar/data/v1.3/L2/2023/10/14/"
timestamps_full_october = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/data/NEID_Data.csv")["filename"]
timestamps_october = timestamps_full_october[15:-150]
wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
airmass = [2.5632073517047638, 2.5361097748885557, 2.5093661650400723, 2.4832872029346307, 2.458152719768879, 2.4333280353218156, 2.4091026139254765, 2.3857380385003344, 2.362645575113234, 2.3400953849743926, 2.3183323449274487, 2.2968091303568645, 2.27577829301206, 2.2554693819371017, 2.2353725338593415, 2.2157241944955177, 2.1967398191699803, 2.1779435542209677, 2.159557077563393, 2.1415684495759097, 2.1241760200750517, 2.106944668512874, 2.0900782010965315, 2.0737632866181173, 2.05759228654005, 2.0417568347767414, 2.0264327849819934, 2.011237707679954, 1.9963519894425743, 1.9819414507550772, 1.9676468497988242, 1.9536381582317257, 1.9400718707467606, 1.926610129678612, 1.9134132425208898, 1.9006290194999498, 1.8879393729586693, 1.8754956248574794, 1.8634374833021545, 1.8514651818936574, 1.8397216717571738, 1.8282016018286706, 1.817034683533013, 1.8059435946574782, 1.7950609742097614, 1.7845095958663924, 1.7740276105655446, 1.763740511098882, 1.753764560828184, 1.743852354616247, 1.7341226871090905, 1.7246857022512947, 1.715307505275213, 1.7061006017579452, 1.697169374222985, 1.688292589769291, 1.6795768389575176, 1.6711210712224671, 1.6627159428484009, 1.6544624719613135, 1.646454476321036, 1.638493796560436, 1.630676192498115, 1.6230906273882857, 1.6155494821322123, 1.6081435466952902, 1.6008703361742513, 1.5938127284826231, 1.5867962627606653, 1.5799055052643187, 1.573219043936084, 1.5665716579824998, 1.5600435303724371, 1.5537091700100665, 1.5474121144362678, 1.5412283821572987, 1.5352285977595637, 1.5292646143956121, 1.523408487955313, 1.5177271385876636, 1.5120803277757258, 1.5065363351237604, 1.5011585426132867, 1.4958142450515572, 1.4905681178810164, 1.4854801586717268, 1.4804248500221022, 1.4754634222942877, 1.4705945147488082, 1.4658738276317298, 1.4611849448608079, 1.456584725994602, 1.4521258025734258, 1.447698240609205, 1.4433557809843884, 1.43914810522267, 1.43497149743696, 1.4308767024770255, 1.4269105616726943, 1.4229753336016377, 1.419118880539379, 1.4153853052103214, 1.4116826158282119, 1.408055897033766, 1.4045466068708823, 1.4010682957774412, 1.3976633677805828, 1.3943707230207358, 1.3911092629152901, 1.387918800523127, 1.3847985802902372, 1.381784212071226, 1.3788014670576347, 1.3758868280098038, 1.3730735241834369, 1.3702922880917996, 1.3675771947689528, 1.364959158943391, 1.3623737296001717, 1.3598526423845454, 1.3574245581270798, 1.3550297095145418, 1.3526975555563139, 1.350454558686985, 1.3482455138644898, 1.346097660961062, 1.3440353146548754, 1.34200772145745, 1.3400399546266697, 1.3381542271029911, 1.3363041363785326, 1.334512636897446, 1.3327998814765274, 1.331123727652466, 1.329505054361501, 1.3279435167148383, 1.3264565754953044, 1.3250076481928261, 1.323614903675427, 1.3222938236556199, 1.3210119084301528, 1.319785336067303, 1.3186276374092976, 1.3175103337752088, 1.3164476420975404, 1.3154511668285525, 1.3144963965600829, 1.3135956129740973, 1.3127585150081922, 1.3119645125498263, 1.3112239742717597, 1.310544710948751, 1.3099100153508805, 1.3093283618247469, 1.3088056868926414, 1.3083291351789732, 1.307905301914301, 1.3075341005197691]
#wavelength = list(np.arange(4292.473, 9300.471, 200))
# airmass_next_day = [1.766610154779519, 1.7566191871652905, 1.746790123204968, 1.7371414103451779, 1.7276693596938506, 1.7183702783550086, 1.7092404757642614, 1.7002766992439282, 1.6914753622423204, 1.6828335117028685, 1.6743475565571195, 1.6660147278317745, 1.6578319322153565, 1.6497962668339659, 1.6419049127001464, 1.6341550389718114, 1.6265442646697401, 1.6190696375641391, 1.6117289205340153, 1.604519582312987, 1.5974392513036884, 1.5904856221310038, 1.5836563716529657, 1.576949485805271, 1.5703627620709284, 1.5638941394660004, 1.5575416133744304, 1.5513031590298254, 1.5451770302626477, 1.5391613052270394, 1.5332541171772707, 1.527453930817201, 1.5217588329551917, 1.5161672390187295, 1.5106774719915244, 1.5052880962432893, 1.499997518005852, 1.4948042509119703, 1.4897068478761184, 1.484703840036097, 1.4797939765912973, 1.4749758607306593, 1.4702481918140997, 1.4656096477761906, 1.4610591075139228, 1.456595313440845, 1.4522170960806213, 1.447923316608493, 1.443712866055501, 1.4395846645406674, 1.4355376120675019, 1.4315707826229884, 1.4276831298091357, 1.423873682942278, 1.4201414969610042, 1.4164856080932275, 1.412905251967709, 1.4093993835768877, 1.4059672829497234, 1.4026080814598816, 1.3959381084413276, 1.392804034379087, 1.3897320542789846, 1.3867291698615183, 1.3837946448524314, 1.38092779997589, 1.3781279749759805, 1.3753944955614459, 1.372726804469675, 1.3701242627422294, 1.367586282409122, 1.3651122927011081, 1.36270173969934, 1.360354086033361, 1.3580687833094516, 1.3558453816115352, 1.3536833634005574, 1.3515822545387943, 1.3495415716912371, 1.3475609436760139, 1.3456398451405287, 1.343777931745373, 1.341974778935331, 1.3402299996988158, 1.3385432207639427, 1.336914063040604, 1.3353422195498081, 1.33382733701681, 1.332369095206067, 1.3309671866996038, 1.3296213167354674, 1.3283312030534131, 1.3270965611490668, 1.3259171631945352, 1.3247927483109057, 1.3237230828412432, 1.3227079330127351, 1.3217471132763106, 1.3208404125328952, 1.3199876438018885, 1.319188631518287, 1.318443202751978, 1.3177512225359715, 1.3171125397011638, 1.3165270235132573, 1.3159945542659572, 1.3155150177391568, 1.3150883325482992, 1.314714395291477, 1.314393131741419, 1.314124484311589, 1.3139083941907213, 1.3137448141726413, 1.3136337210397964, 1.3135750819355776, 1.3135688897745499, 1.3136151423722613, 1.3137138492374163, 1.3138650305503157, 1.3140687171729564, 1.3146508843676523, 1.315013810748451, 1.3154305552945205, 1.3159001247139817, 1.3164226149893354, 1.316998133008446, 1.3176268045440243, 1.3183087432551277, 1.3190440963728511, 1.319833015001365, 1.3206756614719402, 1.3215722206027718, 1.3225228557553914, 1.323527773990482, 1.3245871965910723, 1.3257013039655532, 1.3268703672028908, 1.3280946318919165, 1.3293743240710225, 1.3307097260573355, 1.3321011180265434, 1.333548792632386, 1.3350530551598503, 1.3366142429114671, 1.3382326491657885, 1.339908636631941, 1.3416425627938133, 1.343434798651181, 1.3452857516531398, 1.347195775648958, 1.3491652812378712, 1.3511947680767882, 1.3532846053103853, 1.355435273648915, 1.357647271829831, 1.3599210333055025, 1.362257087032262, 1.3646559514834933, 1.36711816132506, 1.3696442677190188, 1.3722348703643767, 1.374890491710917, 1.3776117652886417, 1.3803993115121311, 1.3832537689908864, 1.3861757948895377, 1.3891661019220714, 1.3922253132611158, 1.395354179930717, 1.3985534374202588, 1.4018238418835418, 1.4051661706042335, 1.4085812642268882, 1.412069861154891, 1.4156328458879297, 1.4192710852249688, 1.4229854696833204, 1.4267769604013154, 1.4306464053695533, 1.434594815180816, 1.443652244037057, 1.4478523219893995, 1.452145161976374, 1.4565223924481816, 1.4609851493925985, 1.4655348215999329, 1.4701722801398585, 1.4748990374443234, 1.4797161316399963, 1.4846250967396304, 1.4896270923593786, 1.4947236116835492, 1.4999158773517323, 1.5052055234645023, 1.5105940403429705, 1.516082824847498, 1.521673510376001, 1.527367708617811, 1.5331671461904635, 1.5390733841591562, 1.5450882402090875, 1.551213510834245, 1.5574510423297185, 1.5638027323458712, 1.570270531513536, 1.5768564450962952, 1.5835626165931016, 1.5903910037189721, 1.5973438668525093, 1.6044234480482498, 1.6116321414868295, 1.6189721452965287, 1.626445986819871, 1.634056270687161, 1.641805303214182, 1.6496960237383178, 1.6577310784209252, 1.6659132873930087, 1.6742456564270367, 1.682730974983207, 1.6913724239481458, 1.7001731776133666, 1.709136616342259, 1.7182658944975635, 1.7275645945624962, 1.7370362993674975, 1.74668470355897, 1.7565136178304737, 1.7665269733501452, 1.7767288263950174, 1.7871234901256565, 1.7977150343901336, 1.8085081771889717, 1.819507130785895, 1.8307171868173617, 1.8421428597377691, 1.8537893655981785, 1.8656619573003812, 1.8777660660016815, 1.8901073085770066, 1.916437987250106, 1.9295189076999746, 1.942887012496054, 1.9565238760263914, 1.9704364259498761, 1.9846318406777248, 1.9991175606773948, 2.013901480914375, 2.028991244782345, 2.0443953298814295, 2.0601223505892463, 2.0761812498315297, 2.092581515191484, 2.1093323979044714, 2.1264447479594546, 2.1439271055675717]

df = NEID_flux(Fiber, path_october, airmass, wavelength, timestamps_october)

# Normalize the values to the range [0, 1]
norm = mcolors.Normalize(vmin=np.min(wavelength), vmax=np.max(wavelength))
# Create a colormap from blue to red
cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

#----------------------------------------

# #SSD
file_data_no_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD/data/neid_october_N_50_KSSD.jld2", "r")
intensity_list_no_ext, intensity_array_no_ext = model_flux(file_data_no_ext)

# fig = plt.figure()
# ax1 = fig.add_subplot()
# for i, value in enumerate(wavelength):
#     color = cmap(norm(value))
#     airmass_arr = df['airmass']
#     comp_arr = ((df[wavelength[i]]/max(df[wavelength[i]]))/(intensity_array[0][i]))
#     ax1.plot(airmass_arr, comp_arr, color=color)
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux/intensity") 
# ax1.axvline(x = df['airmass'][46])
# ax1.axvline(x = df['airmass'][25])
# plt.savefig("figures/KSSD/coeff_eclipse.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()  
# for i, value in enumerate(wavelength):
#     color = cmap(norm(value))
#     ax1.scatter(df['airmass'], file_data[intensity_list[i]][()]/max(file_data[intensity_list[i]][()]), label = wavelength[i], color=color, s = 1)
#     ax1.plot(df['airmass'], df[wavelength[i]]/max(df[wavelength[i]]), color=color)
# rms = round(np.sqrt((np.nansum(((df[wavelength[1]]/max(df[wavelength[1]])) - file_data[intensity_list[1]][()]/max(file_data[intensity_list[1]][()]))**2))/len((df[wavelength[1]]/max(df[wavelength[1]])) - file_data[intensity_list[1]][()]/max(file_data[intensity_list[1]][()]))),4)
# ax1.text(np.array(df['airmass'])[-10], 0.5, "{} RMS {}".format(wavelength[1], rms))    
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux") 
# plt.savefig("figures/KSSD/flux_comp.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()  
# for i, value in enumerate(wavelength):
#     color = cmap(norm(value))
#     ax1.scatter(df['airmass'], df[wavelength[i]]/max(df[wavelength[i]]) - (file_data[intensity_list[i]][()]/max(file_data[intensity_list[i]][()])), label = wavelength[i], color=color, s = 1)
# rms = round(np.sqrt((np.nansum(((df[wavelength[1]]/max(df[wavelength[1]])) - file_data[intensity_list[1]][()]/max(file_data[intensity_list[1]][()]))**2))/len((df[wavelength[1]]/max(df[wavelength[1]])) - file_data[intensity_list[1]][()]/max(file_data[intensity_list[1]][()]))),4)
# ax1.text(np.array(df['airmass'])[-10], 0.5, "{} RMS {}".format(wavelength[1], rms))    
# ax1.set_xlabel("airmass")
# ax1.set_ylabel("relative flux") 
# ax1.axvline(x = df['airmass'][46])
# ax1.axvline(x = df['airmass'][25])
# plt.savefig("figures/KSSD/flux_comp_diff.png")
# plt.clf()

# def float_range(start, end, num_steps):
#     step_size = (end - start) / (num_steps - 1)
#     current = start
#     while current <= end:
#         yield current
#         current += step_size

# min_ext1_arr = []
# min_ext2_arr = []
# min_ext3_arr = []
# norm_ext1_arr = []
# norm_ext2_arr = []

# for i, value in enumerate(wavelength):
#     data_no_clouds = df[wavelength[i]][0:-20]
#     data_arr = data_no_clouds / np.mean(data_no_clouds[115:len(data_no_clouds)])
#     model_intensity = file_data[intensity_list[i]][()][0:-20]

#     min_chi2_ext3 = 10
#     min_ext3 = 0
#     for ext3 in float_range(0.09, 0.55, 50):
#             I_ext2 = model_intensity[47:len(data_arr)] * np.exp(-np.array(airmass[0:-20][47:len(data_arr)])*ext3)
#             ext_model = I_ext2 / np.mean(I_ext2[68:len(data_arr)])

#             chi2_iterate = np.sum((data_arr[47:len(data_arr)] - ext_model) ** 2 / ext_model)
#             if chi2_iterate < min_chi2_ext3:
#                 min_chi2_ext3 = chi2_iterate
#                 min_ext3 = ext3
#     min_ext3_arr.append(min_ext3)

#     min_chi2_ext1 = 10
#     min_ext1 = 0
#     ext1_norm = 0
#     for ext1 in float_range(0.0, 0.1, 200): 
#         for ext2 in float_range(0.0, 0.1, 200): 
#             #first extinction
#             I_ext1 = model_intensity[0:26] * np.exp(-np.array(airmass[0:26])*ext1)
#             #second extinction 
#             I_ext2 = model_intensity[47:len(data_arr)] * np.exp(-np.array(airmass[0:-20][47:len(data_arr)])*ext2)

#             ext_model = I_ext1 / np.mean(I_ext2[68:len(data_arr)])
#             chi2_iterate = np.sum((data_arr[0:26] - ext_model) ** 2 / ext_model)

#             if chi2_iterate < min_chi2_ext1:
#                 min_chi2_ext1 = chi2_iterate
#                 min_ext1 = ext1
#                 ext1_norm = ext2
#     min_ext1_arr.append(min_ext1)
#     norm_ext1_arr.append(ext1_norm)

    # min_chi2_ext2 = 10
    # min_ext2 = 0
    # ext2_norm = 0
    # for ext2 in float_range(0.0, 0.1, 200): 
    #     for ext3 in float_range(0.0, 0.1, 200): 
    #         #first extinction
    #         I_ext1 = model_intensity[26:47] * np.exp(-np.array(airmass[26:47])*ext2)
    #         #second extinction 
    #         I_ext2 = model_intensity[47:len(data_arr)] * np.exp(-np.array(airmass[0:-20][47:len(data_arr)])*ext3)

    #         ext_model = I_ext1 / np.mean(I_ext2[68:len(data_arr)])
    #         chi2_iterate = np.sum((data_arr[26:47] - ext_model) ** 2 / ext_model)

    #         if chi2_iterate < min_chi2_ext2:
    #             min_chi2_ext2 = chi2_iterate
    #             min_ext2 = ext2
    #             ext2_norm = ext3
    # min_ext2_arr.append(min_ext2)
    # norm_ext2_arr.append(ext2_norm)

# fig = plt.figure()
# ax1 = fig.add_subplot()  
# data_no_clouds = df[wavelength[0]][0:-20]
# data = data_no_clouds / np.mean(data_no_clouds[115:len(data_no_clouds)])
# ax1.plot(df['airmass'][0:-20], data)
# m1 = list(model_intensity[0:26] * np.exp(-np.array(airmass[0:26])*min_ext1_arr[0]))
# m2 = list(model_intensity[26:47] * np.exp(-np.array(airmass[26:47])*min_ext2_arr[0]))
# m3 = list(model_intensity[47:len(data_arr)] * np.exp(-np.array(airmass[0:-20][47:len(data_arr)])*min_ext3_arr[0]))
# m1_end = model_intensity[47:len(data_arr)] * np.exp(-np.array(airmass[0:-20][47:len(data_arr)])*norm_ext1_arr[0])
# m2_end = model_intensity[47:len(data_arr)] * np.exp(-np.array(airmass[0:-20][47:len(data_arr)])*norm_ext2_arr[0])
# m1 = list(np.array(m1) / np.mean(m1_end[68:len(data_arr)]))
# m2 = list(np.array(m2) / np.mean(m2_end[68:len(data_arr)]))
# m3 = list(np.array(m3) / np.mean(m3[68:len(data_arr)]))
# ax1.scatter(df['airmass'][0:-20], np.array(m1+m2+m3), s = 1, color = 'red')
# ax1.axvline(x = df['airmass'][47])
# plt.savefig("comp.png")
# plt.clf()

# fig = plt.figure()
# ax1 = fig.add_subplot()  
# ax1.scatter(df['airmass'][0:-20], data - np.array(m1+m2+m3))
# plt.savefig("diff.png")

# # Create a DataFrame
# df = pd.DataFrame({
#     'Wavelength': wavelength,
#     'ext1': min_ext1_arr,
#     'ext2': min_ext2_arr,
#     'ext3': min_ext3_arr,
#     'ext1_norm': norm_ext1_arr,
#     'ext3_norm': norm_ext2_arr
# })
# # Save to CSV file
# df.to_csv('NEID_three_extinctions.csv', index=False)

#----------------------------------------

# SSD + Extinction 
file_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/figures/NEID_October/KSSD/KSSD_3ext/data/neid_october_N_50_KSSD_3ext.jld2", "r")
intensity_list, intensity_array = model_flux_3_ext(file_data)

fig = plt.figure()
ax1 = fig.add_subplot()
for i, value in enumerate(wavelength):
    data_no_clouds = df[wavelength[i]][0:-20]
    data = data_no_clouds / np.mean(data_no_clouds[115:len(data_no_clouds)])
    color = cmap(norm(value))
    ax1.plot(np.delete(df['airmass'][0:-20],46), np.delete(data/(intensity_array[0][i]), 46), color=color)
ax1.set_xlabel("airmass")
ax1.set_ylabel("relative flux/intensity") 
plt.savefig("figures/KSSD_ext/coeff_eclipse.png")
plt.clf()

fig = plt.figure()
ax1 = fig.add_subplot()  
for i, value in enumerate(wavelength):
    data_no_clouds = df[wavelength[i]][0:-20]
    data = data_no_clouds / np.mean(data_no_clouds[115:len(data_no_clouds)])
    color = cmap(norm(value))
    ax1.scatter(np.delete(df['airmass'][0:-20],46), np.delete(intensity_array[0][i],46), label = wavelength[i], color=color, s = 1)
    ax1.plot(np.delete(df['airmass'][0:-20],46), np.delete(data, 46), color=color)
data_no_clouds = df[wavelength[0]][0:-20]
data = data_no_clouds / np.mean(data_no_clouds[115:len(data_no_clouds)])
rms = round(np.sqrt((np.nansum(np.delete(data - intensity_array[0][i],46)**2))/len((np.delete(data,46)))),4)
ax1.text(np.array(df['airmass'])[-10], 0.5, "{} RMS {}".format(wavelength[1], rms))    
ax1.set_xlabel("airmass")
ax1.set_ylabel("relative flux") 
plt.savefig("figures/KSSD_ext/flux_comp.png")
plt.clf()

fig = plt.figure()
ax1 = fig.add_subplot()  
for i, value in enumerate(wavelength):
    data_no_clouds = df[wavelength[i]][0:-20]
    data = data_no_clouds / np.mean(data_no_clouds[115:len(data_no_clouds)])
    color = cmap(norm(value))
    ax1.scatter(np.delete(df['airmass'][0:-20],46), np.delete(data - intensity_array[0][i],46), label = wavelength[i], color=color, s = 1)  
ax1.set_xlabel("airmass")
ax1.set_ylabel("relative flux") 
plt.savefig("figures/KSSD_ext/flux_comp_diff.png")
plt.clf()

#----------------------------------------

# df_pyrrheliometer = pd.read_csv("Data/neid_ljpyrohelio_chv0_20231009.tel.txt", sep=" ", names=["Date", "Photocell", "Irradiance"])
# df_pyrrheliometer["DateTime"] = [datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f") for date in df_pyrrheliometer['Date']]
# df_pyrrheliometer["Day"] = [str(datetime.strptime(date, "%Y-%m-%dT%H:%M:%S.%f").date()) for date in df_pyrrheliometer['Date']]
# df_pyrrheliometer["Hour"] = pd.to_datetime(df_pyrrheliometer['DateTime']).dt.hour 
# df = df_pyrrheliometer.loc[(df_pyrrheliometer['Day'] == "2023-10-14") & (df_pyrrheliometer['Hour'] >= 15)]

# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.plot(np.array(df["DateTime"]), np.array(df['Irradiance']))
# ax1.set_ylabel("Irradiance (W/m^2)") 
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("hour on 10/14")
# plt.savefig("figures/pyrrheliometer.png")
# plt.clf()

#----------------------------------------

# for file in os.listdir(directory):
#     filename = os.fsdecode(file)
#     if filename.endswith(".fits"): 
#         inputSpectrum = fits.open('Data/{}'.format(filename))
#         # inputSpectrum.info() # shows file extensions/info
#         # print(inputSpectrum[0].header)
#         # sciflux=inputSpectrum[1].data#[]
#         # scilambda=inputSpectrum[7].data#[]
#         # print(sci.shape)

#         # plt.plot(np.concatenate(scilambda),np.concatenate(sciflux),label="original")
#         # plt.legend()
#         # plt.show()

#         # plt.figure()
#         # plt.imshow(scilambda, aspect='auto', interpolation='nearest')
#         # plt.colorbar()
#         # plt.savefig('wavelength.png')

#         # plt.figure()
#         # plt.imshow(sciflux, aspect='auto', vmin=-10, vmax=10)
#         # plt.colorbar()
#         # plt.savefig('flux.png')
