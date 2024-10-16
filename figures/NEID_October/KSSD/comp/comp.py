import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
import matplotlib.colors as mcolors
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

#read in data
#GRASS
grass_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_all_lines_rv_regular_KSSD.jld2", "r")
lines = grass_data["name"][()]
GRASS_rv  = grass_data["rv"][()]
grass_data_no_cb = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_all_lines_rv_off_KSSD.jld2", "r")
lines_no_cb = grass_data_no_cb["name"][()]
GRASS_no_cb  = grass_data_no_cb["rv"][()]
#model
file = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_october_N_50_KSSD.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()]

#data 
data = pd.read_csv("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/data/NEID_Data.csv")
rv_obs = list(data["ccfrvmod"][15:-150]*1000 + 644.9)
UTC_time = []
time_julian = []
for i in data["obsdate"][15:-150]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)
#line by line data
line_data = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD/data/neid_RVlinebyline_KSSD.jld2", "r")
line_lines = line_data["name"][()]
line_rv  = line_data["rv"][()]

vb, warnings, flag = get_BC_vel(JDUTC=time_julian[0:-25], lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs[0:-25])
rv_obs -= rv_obs[-1]

UTC_time = UTC_time[0:-25]

#read in data with extinction 
#GRASS
grass_data_ext = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_all_lines_rv_regular_KSSD_2_ext.jld2", "r")
lines_ext  = grass_data_ext["name"][()]
GRASS_rv_ext   = grass_data_ext["rv"][()]
grass_data_no_cb_ext  = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_all_lines_rv_off_KSSD_2_ext.jld2", "r")
lines_no_cb_ext  = grass_data_no_cb_ext["name"][()]
GRASS_no_cb_ext   = grass_data_no_cb_ext["rv"][()]
#model
file_ext  = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_october_N_50_KSSD_2_ext.jld2", "r")
RV_list_no_cb_ext  = file_ext["RV_list_no_cb"][()]
#line by line data
line_data_ext  = h5py.File("/storage/home/efg5335/work/Eclipse_GRASS/src/plots/NEID_October/KSSD/KSSD_ext/data/neid_RVlinebyline_KSSD_2_ext.jld2", "r")
line_lines_ext  = line_data_ext["name"][()]
line_rv_ext   = line_data_ext["rv"][()]

def jld2_read(jld2_file, variable, vb, index):
    array = jld2_file[variable[index]][()][0:-25]
    array = np.array(array + vb)
    # array -= array[-1]
    return array

wavelength = [5250.2084, 5250.6453, 5379.5734, 5381.0216, 5382.2562, 5383.368, 5432.546, 5432.947, 5434.5232, 5435.8577, 5436.2945, 5436.5875, 5576.0881, 5578.718, 6149.246, 6151.617, 6169.042, 6169.563, 6170.5056, 6173.3344, 6301.5008, 6302.4932]
airmass = [2.572425367171159, 2.5448735167307444, 2.5180129852674447, 2.4918196512704975, 2.4662705139035777, 2.4413436288441934, 2.4170177578145795, 2.3932734823889623, 2.37009138693885, 2.3474531938592382, 2.3253414202609437, 2.3037393354077356, 2.2826306684738125, 2.2620005864480275, 2.2418338863836254, 2.222117195250764, 2.202836269389674, 2.1839783769544594, 2.165530844281854, 2.147482158803499, 2.1298206065409873, 2.1125351465601323, 2.0956149615094435, 2.079050251585693, 2.062830987750628, 2.046947716592267, 2.0313913365785403, 2.016153081856585, 2.0012245071345185, 1.986597298269381, 1.972263962461992, 1.958216754226753, 1.9444483742854706, 1.9309517792954787, 1.9177200124152, 1.9047468298365413, 1.8920257307863801, 1.8795505914120514, 1.8673154942382788, 1.8553147196941986, 1.8435427380566671, 1.831994063434159, 1.820663666748346, 1.809546809503011, 1.798638241102934, 1.7879334089554906, 1.777427647870953, 1.7671168159344612, 1.7569965267146577, 1.7470626546708954, 1.7373112019147001, 1.727738293459651, 1.718340060053852, 1.7091130864141855, 1.7000537251133014, 1.691158549369504, 1.6824242345351446, 1.6738474516746664, 1.665425277118639, 1.6571545668250862, 1.6490323699354348, 1.6410558210923736, 1.6332221375080378, 1.6255285239477424, 1.6179725405672654, 1.6105515420578271, 1.6032630490559034, 1.5961046518781352, 1.5890740082488395, 1.5821688411184676, 1.5753868552864425, 1.5687260619706278, 1.5621842847995455, 1.5557594874908791, 1.5494496136755136, 1.5432528890610906, 1.5371673619055677, 1.5311912080589878, 1.5253226524609174, 1.5195599676196931, 1.5139014722332405, 1.5083454632407058, 1.5028904820389442, 1.49753491009505, 1.4922772377184867, 1.487115995081263, 1.4820497510829715, 1.4770771122931756, 1.4721966633832764, 1.4674072011461599, 1.4627073793625198, 1.4580959450772295, 1.4535716239098995, 1.449133337013963, 1.444779871923669, 1.4405100505253015, 1.436322928737646, 1.432217235163834, 1.428192079430581, 1.4242463978391648, 1.4203792034339793, 1.4165895351965723, 1.412876457412079, 1.4092390590575627, 1.4056764105696757, 1.402187734730558, 1.398772147589441, 1.3954288311855378, 1.3921569503535707, 1.388955847927697, 1.3858246154153044, 1.3827626345546815, 1.3797691539915007, 1.376843480268806, 1.373984904963574, 1.3711928417826138, 1.368466618713325, 1.3658056167673367, 1.3632092346354667, 1.3606768883388554, 1.3582080108923182, 1.3558020519795388, 1.3534584496721638, 1.3511767427430967, 1.348956373841781, 1.3467969357634277, 1.3446978525698496, 1.342658745016126, 1.3406791455438265, 1.3387586268041718, 1.3368967756697012, 1.3350931930210945, 1.3333474727522334, 1.331659285322119, 1.3300282508874752, 1.32845402425539, 1.3269362732391417, 1.325474678483048, 1.324068916628486, 1.3227187274887562, 1.321423811884163, 1.320183900178243, 1.3189987347982757, 1.3178680701008785, 1.3167916722521957, 1.3157693070586747, 1.3148007886685817, 1.313885905273252, 1.3130244691189576, 1.312216303719761, 1.311461234956191, 1.3107591268798375, 1.3101098268427607, 1.309513202596138, 1.3089691327704955, 1.3084775011853784, 1.3080382199450036, 1.3076511937122919]

# Normalize the values to the range [0, 1]
norm = mcolors.Normalize(vmin=np.min(wavelength), vmax=np.max(wavelength))
# Create a colormap from blue to red
cmap = plt.get_cmap('coolwarm')  # or 'RdYlBu' or any other suitable colormap

def plot(variable, file, label, savfile, wavelength):
    fig = plt.figure()
    ax1 = fig.add_subplot()
    for i, value in enumerate(wavelength):
        if i == 5:
            continue
        color = cmap(norm(value))

        arr = jld2_read(file, variable, vb, i)

        ax1.plot(UTC_time, arr, color = color, linewidth = 2)

    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax1.set_xlabel("Time (UTC)", fontsize=12)
    ax1.set_ylabel("RV [m/s]", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(label)
    plt.savefig(savfile)
    plt.clf()

def plot_diff(variable1, file1, variable2, file2, label, savfile, wavelength):
    fig = plt.figure()
    ax1 = fig.add_subplot()
    for i, value in enumerate(wavelength):
        if i == 5:
            continue
        color = cmap(norm(value))

        arr1 = jld2_read(file1, variable1, vb, i)
        arr2 = jld2_read(file2, variable2, vb, i)

        ax1.plot(UTC_time, arr1 - arr2, color = color, linewidth = 2)

    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax1.set_xlabel("Time (UTC)", fontsize=12)
    ax1.set_ylabel("RV [m/s]", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(label)
    plt.savefig(savfile)
    plt.clf()

plot(RV_list_no_cb, file, "Projected - no cb (no ext)", "projected_no_ext.png", wavelength)
plot(RV_list_no_cb_ext, file_ext, "Projected - no cb (ext)", "projected_ext.png", wavelength)
plot_diff(RV_list_no_cb_ext, file_ext, RV_list_no_cb, file, "Projected no cb (ext - no ext)", "projected_ext_diff.png", wavelength)

plot(GRASS_rv, grass_data, "GRASS - cb (no ext)", "grass_cb_no_ext.png", wavelength)
plot(GRASS_rv_ext, grass_data_ext, "GRASS - cb (ext)", "grass_cb_ext.png", wavelength)
plot_diff(GRASS_rv_ext, grass_data_ext, GRASS_rv, grass_data, "GRASS cb (ext - no ext)", "grass_cb_ext_diff.png", wavelength)

plot(GRASS_no_cb, grass_data_no_cb, "GRASS - no cb (no ext)", "grass_no_cb_no_ext.png", wavelength)
plot(GRASS_no_cb_ext, grass_data_no_cb_ext, "GRASS - no cb (ext)", "grass_no_cb_ext.png", wavelength)
plot_diff(GRASS_no_cb_ext, grass_data_no_cb_ext, GRASS_no_cb, grass_data_no_cb, "GRASS no cb (ext - no ext)", "grass_no_cb_ext_diff.png", wavelength)

plot(line_rv, line_data, "Line (no ext)", "line_no_ext.png", wavelength)
plot(line_rv_ext, line_data_ext, "Line (ext)", "line_ext.png", wavelength)
