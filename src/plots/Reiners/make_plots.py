import h5py
import math
import csv
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt 
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from astropy.time import Time

#GRASS_rv = [185.22079487847145, 190.13037877392662, 188.27736197141184, 189.04543213216084, 190.308609257324, 190.3220785667663, 195.02263175690632, 193.50543868423844, 192.6093953203454, 192.91878039051295, 197.15646395874964, 195.10250588307508, 202.5011251287988, 203.83324093756943, 200.84227277585725, 201.81958558700902, 203.0075113726937, 204.38935499480843, 204.9464960461308, 207.884361511744, 210.07434834562278, 210.92553964989295, 209.750638178996, 213.5208525090023, 213.41851399145904, 214.43759863100624, 214.82496816247001, 219.02719939125268, 221.62858554856012, 219.30223044506508, 227.05684063263354, 223.54001004327458, 225.68435891039243, 224.57190374807692, 226.19316011580162, 229.6157731427279, 228.5190228310273, 233.65674452707674, 231.84067817639934, 237.87107663047982, 237.6508018026791, 241.56971892021983, 239.29422048123462, 242.48090631095158, 243.8544737372461, 245.25736976132444, 252.48540255847095, 250.25504966626204, 248.68877459457252, 241.1952085956524, 230.78872387295, 214.9600052670899, 198.5980502053981, 183.7163660941176, 164.92426559931957, 139.70585591435253, 110.88727107235715, 87.68317402203462, 56.02853894076872, 31.496414120834537, 4.274045744195973, -25.721896423155517, -60.34984154622469, -88.74425340141472, -123.54304941087554, -148.29718509677613, -181.19198030231374, -208.76540145279515, -241.6594113360246, -263.4729018962027, -288.851061596824, -316.46589010484263, -330.4599544178024, -337.07031027765913, -345.86929435433876, -347.6354530119901, -340.73407075863935, -320.01939443416757, -286.0115551309078, -229.8094030513495, -157.4517936834772, -77.63445170674396, 46.83606636786106, 186.12526490428377, 318.1896614159581, 457.8354948136906, 595.7363779461443, 720.2376374017119, 813.8208011250945, 892.9561613988197, 952.1581335034615, 983.7217001259189, 1016.4175210855004, 1026.2220770313857, 1037.536292642075, 1029.0574252012127, 1013.3781760075038, 996.981616691851, 981.6162722066422, 955.2982628392415, 921.5739337107699, 888.3753688830185, 864.7206645192479, 833.0333163739685, 802.8869607268857, 774.9724096635215, 742.6894133361378, 711.2255904851901, 679.5684027866535, 648.4307952569429, 619.1994876570518, 593.4614835498061, 566.5338037305227, 545.3938655996392, 521.1866227238324, 504.99294841342595, 479.4438235194204, 458.6724084340225, 443.38582111728203, 428.54379425150523, 421.91917987616324, 411.35504574319845, 406.58640205279744, 406.2324578379578, 406.13523160139533, 405.69229995061414, 411.3563383845905, 410.86151804337914, 414.1622631258724, 419.11761502951947, 421.8997907889057, 419.0813616066156, 426.7776398663327, 426.36843530131955, 429.19688415039434, 432.5813556474184, 436.02574119119316, 436.9185946080879, 437.200610376281, 438.7267903116193, 442.80831510027804, 441.26380557771313, 445.7835956357488, 450.0256831922874, 451.38848544063114, 452.66117927596633, 455.76649671724647, 459.2626773395419, 457.38216075453613, 461.75079799686534, 479.39094167146055, 478.7954812998366, 479.88859287046085, 483.3665254316129, 482.00916329581264, 485.8289692923192, 489.2550216534782, 486.1292794927257, 492.6894710114944]
#GRASS_rv = [82.24509865307951, 81.24498087661793, 80.84060920236968, 83.68609316922364, 85.37121891013193, 86.58045154340617, 84.90889322584879, 87.54506019062002, 90.48870060284352, 91.03382016026373, 95.53001818700436, 93.67310895607194, 95.72891657746763, 95.98454568883017, 94.44763815426579, 97.91474592702512, 99.0471377516749, 101.30086967519298, 101.40932264835833, 104.88050925278317, 104.46919568023509, 104.80374772206255, 111.01980166740047, 104.66637964183191, 110.22859144503289, 109.63348523803691, 113.45098876019122, 115.31602513017606, 113.20368575332726, 114.20382826633009, 112.7428937905931, 119.54166849133654, 120.24121982932881, 122.97539960375649, 123.60203172644584, 125.28046964184904, 123.14637084409807, 129.22927475938366, 128.64863359221857, 129.2204688403186, 132.725643645661, 133.94310474649092, 135.61915662936613, 138.51086913340453, 137.88062910984402, 134.74852476063108, 141.76252161570838, 145.56675425376142, 141.954573915278, 133.98344059747808, 121.68350670182113, 104.7662921724926, 90.60312834573092, 67.17823006835, 44.87335877563885, 18.153557444889596, -4.928364941894796, -34.30957936915593, -62.194636035224946, -89.56867817989314, -115.23671799321896, -146.08944767681538, -178.36587279159792, -206.06232722095626, -237.1708862380222, -261.703917635701, -293.5080777292452, -317.6246936987115, -342.33004039012036, -374.28332744712657, -390.45203818371954, -413.13903925649527, -437.8924074178846, -439.23516951718, -443.63405034972783, -445.12625617158636, -435.8420912796106, -412.09092477298174, -370.6163933932659, -321.36251531081615, -249.77176312494157, -164.48126841370274, -46.595995268706204, 86.9661224496431, 219.10184717668054, 356.67663066543037, 495.6485608905832, 616.8499345675807, 701.8396823719261, 776.8101827005252, 833.2339294186411, 877.0765801836652, 901.4887870253574, 916.5951377232944, 935.2346339617627, 924.5999100288129, 910.8245172950136, 892.6493261798875, 873.321514471984, 852.0626597975878, 814.7237005283047, 790.0151878033953, 762.0390618415083, 730.0393486641852, 706.8027548661705, 673.6000772204626, 643.779240480882, 612.9338988170279, 580.271070394082, 545.6216268375252, 526.8579672410083, 497.19902872050676, 474.4654097920989, 452.3270239677031, 417.86404183544397, 402.91358179973633, 379.8905881936049, 358.99842600461517, 343.4115213274059, 329.3689815215206, 313.04480326181556, 307.293061024825, 299.02293826899927, 295.1904476865396, 303.10241977697586, 299.9136192307739, 304.2597161959558, 306.55562344749256, 304.3176578555417, 311.28653809554703, 313.1506293425462, 318.7593351141312, 318.16070954027964, 320.40853636207623, 324.34872602742183, 323.7174559592192, 325.8986345969053, 330.5845115983359, 333.1679919650041, 332.9160718147737, 336.34347663881545, 336.6672015762774, 340.0113979038265, 342.3228975759106, 347.35242972251496, 346.39546719988596, 345.9472606410323, 349.8276951617625, 352.9219044590745, 355.2543419549156, 368.5768167083609, 371.3878936771033, 372.9373363633721, 373.03252476036124, 377.78257588740206, 378.68660703781904, 379.64025555469544, 384.40610075062636, 388.1399830866807]
GRASS_rv = [185.55524612668583, 184.1137494998501, 186.94954436150834, 186.6706675113734, 186.6919078765123, 187.83582651261756, 195.22408398226608, 192.38013894867728, 192.2506688802623, 193.02074588556403, 199.22091404587744, 196.197096121704, 194.76583480966892, 197.40028820435577, 202.53428702207077, 201.53555237575063, 201.34785900226504, 206.1405732503619, 204.6739142673757, 207.7991997779799, 207.41633783409137, 210.93746908047444, 210.08408008640788, 208.22248960444296, 211.69980295011558, 211.5730645645544, 216.33833834066468, 216.41914250719736, 219.37831656999145, 218.68409251612297, 220.59547903817761, 224.89306976864296, 224.569658598629, 221.5621161303031, 226.8754392253254, 232.91135831583807, 230.1376982885175, 230.40724005150176, 235.34403987756133, 233.31391203038936, 236.9195053012098, 239.48177289965432, 242.98663244322114, 239.4824441008731, 244.44137151406105, 242.4812867406585, 245.02755457054613, 249.61950930647924, 246.3728286999408, 237.98036164167047, 224.79142538883022, 202.5332377387986, 187.76749690113004, 169.52587229298243, 145.30107373468107, 116.61553566298451, 94.69873766437847, 66.13821115307663, 38.78061125750478, 7.472545898851981, -22.58777525147849, -53.114797295760695, -89.43167315587516, -120.71427508921722, -155.3269349099154, -182.08368635039926, -209.01350025254627, -241.74745941102717, -272.0721933561039, -296.5744070818533, -328.65163861289244, -341.7773048761867, -360.0010125471359, -372.43741827695817, -369.6979614020139, -373.3142777729079, -359.67749540798894, -333.0483848548124, -295.2872103404558, -232.75735696077612, -155.73528699465194, -58.074937811961384, 55.1673050366491, 176.6409457083123, 320.1516587559031, 457.76899733066574, 596.0265946113844, 708.6027152260364, 813.866620130112, 898.3920772367302, 961.3921703862617, 1005.9180150242154, 1038.8602397633322, 1050.723442219282, 1059.6541527226595, 1054.49202316837, 1050.715339495766, 1032.4912290978473, 1007.158164186087, 981.0520196372952, 950.9835723153774, 921.3209753291566, 893.3396975039103, 864.0543104743716, 830.3236419340526, 798.6690064172947, 764.4173163573632, 738.1597824994338, 707.6011467952809, 675.9463226712661, 648.0869738697758, 621.5245547858754, 588.8767973511049, 562.4617876882869, 538.0799435057301, 516.423206983711, 496.1791203163919, 472.6116090015251, 456.17393637281407, 440.2488163976377, 423.8090373176268, 414.9427485374574, 404.6826085707622, 404.72332488953293, 401.5277462455501, 405.847676575882, 410.8583729296204, 408.389088132633, 411.5195375177547, 415.38077773891837, 416.14342161191945, 422.2223808622882, 419.20675255037867, 422.5313874223972, 428.9980001359842, 432.1881607474027, 432.87469835334207, 433.76961552161004, 435.7868221558594, 437.6782273738152, 440.54528044188913, 440.59791965215874, 443.36683856322355, 447.0416038527263, 448.27170340535224, 450.90160778723754, 450.7569514301252, 457.3648178528346, 456.34813718159063, 459.2576859129441, 475.5930508733941, 474.51705739244926, 479.12204061426195, 477.1637410224978, 482.32501766715154, 483.3709681751527, 486.2059917645307, 488.6466047808169, 490.2355442971168]

#read in data
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()] 
RV_list_cb  = file["RV_list_cb"][()]
RV_list_cb_new  = file["RV_list_cb_new"][()]
intensity_list = file["intensity_list"][()]

#data 
f=open("data/Reiners_Data.txt","r")
lines=f.readlines()[1:]
time=[]
raw_rv = []
time_julian = []
UTC_time = []
for x in lines:
    dt = datetime.strptime("2015-03-20 {}".format((x.split()[1])), "%Y-%m-%d %H:%M:%S.%f")
    time.append(dt)
    time_julian.append((Time(dt)).jd)
    raw_rv.append(float(x.split()[4]))
for x in range(0,len(lines)):
    UTC_time.append((datetime.strptime("2015-03-20 {}".format((lines[x].split()[1])), "%Y-%m-%d %H:%M:%S.%f") + timedelta(seconds=float(50))))
f.close()
#reiners model residuals and intensity 
with open('data/model2_residuals.csv') as f:
    cf = csv.reader(f)
    time_residuals_model2 = []
    residual_model2 = []
    for row in cf:
        time_residuals_model2.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_model2.append(float(row[1]))
with open('data/new_residuals_model1.csv') as f:
    cf = csv.reader(f)
    time_residuals_new1 = []
    residual_new1 = []
    for row in cf:
        time_residuals_new1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        residual_new1.append(float(row[1]))
with open('data/flux_data.csv') as f:
    cf = csv.reader(f)
    time_flux_model1 = []
    flux_model1 = []
    for row in cf:
        time_flux_model1.append(datetime.strptime("2015-03-20 {}".format(row[0]), "%Y-%m-%d %H:%M"))
        flux_model1.append(float(row[1]))     

raw_rv = np.array(raw_rv)
raw_rv -= raw_rv[-1]

GRASS_rv = np.array(GRASS_rv)
GRASS_rv -= GRASS_rv[-1]

RV_list_no_cb = np.array(RV_list_no_cb) 
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb)
RV_list_cb -= RV_list_cb[-1]

RV_list_cb_new = np.array(RV_list_cb_new)
RV_list_cb_new -= RV_list_cb_new[-1]

residual_slope = ((raw_rv[46] - RV_list_cb[46])-(raw_rv[0] - RV_list_cb[0]))/46
linear_correction = residual_slope*range(0,46) + (raw_rv[0] - RV_list_cb[0])

corrected_rv = list((raw_rv[0:46]) - linear_correction) + list(raw_rv[46:,])

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_no_cb)**2))/len(corrected_rv - RV_list_no_cb)),2)
# axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals
# axs[1].scatter(UTC_time, corrected_rv - RV_list_no_cb, color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_new1, residual_new1, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals.png")
# plt.show()

#rm curve
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 2, label = "GRASS")
axs[0].plot(UTC_time, RV_list_no_cb, color = 'r', label = "Model - No CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_no_cb)**2))/len(corrected_rv - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((corrected_rv - GRASS_rv)**2))/len(corrected_rv - GRASS_rv)),2)
axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
axs[0].legend()
#residuals
axs[1].scatter(UTC_time, corrected_rv - RV_list_no_cb, color = 'r', marker = "x", s = 1)
axs[1].scatter(time_residuals_new1, residual_new1, color = 'k', marker = "x", s = 1)
axs[1].scatter(UTC_time, corrected_rv - GRASS_rv, color = 'b', marker = "x", s = 3)
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals")
plt.savefig("rm_and_residuals_no_cb2.png")
plt.show()

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_cb, color = 'r', label = "Model - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_cb)**2))/len(corrected_rv - RV_list_cb)),2)
# axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals       
# axs[1].scatter(UTC_time, corrected_rv - RV_list_cb, color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_cb.png")
# plt.show()

# #rm curve
# fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
# axs[0].scatter(UTC_time, corrected_rv, color = 'k', marker = "x", s = 15, label = "Reiners RVs")
# axs[0].plot(UTC_time, RV_list_cb_new, color = 'r', label = "Model - CB")
# axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# rms_model_no_cb = round(np.sqrt((np.nansum((corrected_rv - RV_list_cb_new)**2))/len(corrected_rv - RV_list_cb_new)),2)
# axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
# axs[0].set_xlabel("Time (UTC)")
# axs[0].set_ylabel("RV [m/s]")
# axs[0].legend()
# #residuals
# axs[1].scatter(UTC_time, corrected_rv - RV_list_cb_new, color = 'r', marker = "x", s = 1)
# axs[1].scatter(time_residuals_model2, residual_model2, color = 'k', marker = "x", s = 1)
# axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# axs[1].set_xlabel("Time (UTC)")
# axs[1].set_ylabel("Residuals")
# plt.savefig("rm_and_residuals_cb_new.png")
# plt.show()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.plot(UTC_time, intensity_list/intensity_list[129], color = 'r', label = "model")
# ax1.scatter(time_flux_model1, flux_model1, color = 'k', label = "Reiners")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity")
# plt.legend()
# plt.savefig("intensity.png")
# plt.show()

# #barycentric correction test: Sun center vs surface
# vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=51.560583 , longi=9.944333, alt=201, SolSystemTarget='Sun', predictive=False,zmeas=0.0)
# fig = plt.figure()
# ax1 = fig.add_subplot()
# plt.scatter(time, vb, color = 'g', label = "BC using center")
# plt.scatter(time, bc_list, color = 'r', label = "BC using surface")
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Barycentric Correction")
# plt.legend()
# plt.savefig("barycentric_correction.png")
# plt.show()

# #extiction test: differential vs absolute 
# f=open("Reiners_Data.txt","r")
# lines=f.readlines()[1:]
# airmass = []
# for x in range(0,len(lines)):
#    airmass.append(float((lines[x].split()[2])))
# f.close()
# ext_coefficient = list(np.linspace(0.4, 0.1, num = len(airmass)))
# fig = plt.figure()
# ax1 = fig.add_subplot()
# plt.scatter(time, airmass, color = 'b', label = "absolute airmass")
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Airmass")
# diff_airmass1 = []
# diff_airmass2 = []
# for i in range(0,len(airmass_matrix)):
#     diff_airmass1.append(1 / file[airmass_matrix[i]][()][1][1])
#     diff_airmass2.append(1 / file[airmass_matrix[i]][()][30][10])
# plt.scatter(time, diff_airmass1, color = 'r', label = "airmass for random patch")
# plt.scatter(time, diff_airmass2, color = 'g', label = "airmass for random patch")
# plt.legend()
# plt.savefig("airmass_test")
# plt.show()
