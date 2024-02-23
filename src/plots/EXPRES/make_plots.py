import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

GRASS_rv = [-975.0464667249638, -971.9646375548864, -971.9207030866073, -970.7827525481874, -966.6617539254989, -964.061117063012, -968.6577223159743, -964.026802613031, -958.8009420981863, -963.5119013948369, -959.3592886641252, -954.5530045408425, -954.5802757748658, -951.9502116679928, -951.208449942905, -949.8521199119386, -947.8525209566948, -947.5692811772909, -945.5888454840273, -943.7191628449236, -942.5023574914884, -941.2871269967877, -940.2035471553523, -939.1741741819316, -937.0012284326455, -934.0147598926922, -943.6190509438547, -954.0801411510315, -968.0902280862207, -984.7203882004258, -999.051321940244, -1023.6528364238034, -1047.7237595092463, -1067.8677449351137, -1093.37881589661, -1122.5588378125549, -1154.8762038346226, -1185.0331835947222, -1217.9833254972818, -1258.136171108366, -1293.8560091289614, -1325.7468955849472, -1371.3962882832923, -1409.4392353059432, -1466.3437178144386, -1513.9665472377853, -1557.1605772231756, -1596.5697148800045, -1619.8338861050709, -1625.5741408718886, -1588.9691043922346, -1448.0476758874997, -1147.61916987093, -583.3571981457955, -81.846312648811, 139.289368741028, 178.2070543969087, 164.06724503587162, 102.35471574301226, 41.85082818980176, -20.459623665617393, -80.97374997761743, -126.9879689843838, -184.69174581876663, -223.6224588511426, -263.26875590492307, -307.11145045348377, -344.09934812876986, -374.378718678369, -407.03278629050607, -441.6565994962974, -464.7620477462592, -491.0777091117972, -517.428979093915, -560.5052227256415, -586.6273050144432, -618.4122788033105, -646.525826115773, -664.882250498285, -687.3040548264875, -699.567489585619, -714.0865726138838, -716.8086792936016, -722.7356084360835, -720.8009173210739, -716.6196273202098, -710.9999042848183, -706.5873704459863, -701.8501902647864, -693.703905345073, -693.3751454988011, -690.5877627094304, -689.5247517540377, -677.4559543840229, -675.1709673966037, -675.2030371468998, -670.874293860989, -665.6631683630621, -661.6004939110472, -653.8529734337778, -650.7760392160151, -646.1065306654741, -643.7738686486118, -639.1894230869408, -635.722804964729, -633.6741240711237, -626.9833690337381, -626.3298956051486, -618.4756124648143, -614.9055208631536, -612.8948815399332, -608.0931549197711, -604.4355282295791, -599.0962928106646, -597.5381789663447, -590.9314725654688, -586.3078215510135, -581.6172537567734, -576.4938729112404, -577.9499884975191, -573.6549273843411, -565.1933619237309, -561.6898507224253, -560.0384189272769, -554.9547145723292, -553.2549138912394, -546.1428390965967, -542.1753527769359, -537.4645681451489, -540.7296636624422, -527.7535677782845, -526.1167669006932, -524.125624179013, -519.7622736377342, -513.4789910049094, -515.12940918166, -508.0555722641301, -504.931983867474, -498.2981864301512, -499.5246579284195, -490.63931716268036, -491.1645435183781, -482.7294628878588, -477.4809757309848, -479.4706198320252, -473.54198145225706, -472.5882060125817, -464.7799810259104, -462.16137715599035, -455.49034351768483, -455.02103593440194, -449.8141546320531, -444.21007737975975, -444.0432591510501, -442.64763763097875, -435.69727233751945, -434.6303732634703, -427.47795956749593, -426.0763788504863, -421.46109779490615, -419.2780437605392, -415.54468757770246, -414.1243735823258, -409.69922757554247, -405.9730373428257, -398.32207636026806, -400.7744359386828, -393.9090405754223, -389.5273185223373, -387.45945042724026, -384.7744728922968, -383.8177346319311, -377.4028459805479, -375.50501213607964, -373.71987496133823, -371.6498271377727, -366.0927570096379, -363.66567440823405, -358.13898896460796, -355.2507355037557, -349.95528246652907, -351.51623655140645, -343.02264855483185, -343.0059837811158, -340.73798684637694, -336.153550398437, -337.87361883850696, -329.4220797766438, -332.33144180510106, -321.4961031605616, -321.6426272233646, -317.56264639429185, -316.58520563067185, -313.53102057803375, -309.4371835241541, -307.57077871552394, -302.99281598513573, -303.00082521575854, -302.2703342261285, -296.80576662792794, -294.44315606094955, -291.1160328069923]

#read in data
#model
file = h5py.File("model_data.jld2", "r")
RV_list_no_cb = file["RV_list_no_cb"][()][20:-100]
RV_list_cb  = file["RV_list_cb"][()][20:-100]
intensity_list = file["intensity_list"][()][20:-100]
vel_no_cb = file["vel_no_cb"][()]
vel_cb = file["vel_cb"][()]
#data 
data = pd.read_csv("EXPRES_data.csv")
rv_obs = list(data["rv"][20:-100])
UTC_time = []
time_julian = []
for i in data["tobs"][20:-100]:
    dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S")
    UTC_time.append(dt)
    time_julian.append((Time(dt)).jd)

vb, warnings, flag = get_BC_vel(JDUTC=time_julian, lat=34.744444 , longi=-111.421944 , alt=235.9152, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

rv_obs = np.array(rv_obs)
rv_obs -= rv_obs[-1]

GRASS_rv = np.array(GRASS_rv[20:-100] + vb)
GRASS_rv -= GRASS_rv[-1]

RV_list_no_cb = np.array(RV_list_no_cb + vb)
RV_list_no_cb -= RV_list_no_cb[-1]

RV_list_cb = np.array(RV_list_cb + vb)
RV_list_cb -= RV_list_cb[-1]    

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 15, label = "EXPRES RVs") 
axs[0].plot(UTC_time, RV_list_no_cb , color = 'r',  label = "Model - No CB")
axs[0].plot(UTC_time, GRASS_rv, color = 'b', linewidth = 3, label = "GRASS")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
rms_model_no_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_no_cb)**2))/len(rv_obs - RV_list_no_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_no_cb))
rms_grass_no_cb = round(np.sqrt((np.nansum((rv_obs - GRASS_rv)**2))/len(rv_obs - GRASS_rv)),2)
axs[0].text(UTC_time[-40], -500, "GRASS RMS {}".format(rms_grass_no_cb))
axs[0].legend()
#residuals
axs[1].scatter(UTC_time, (rv_obs) - RV_list_no_cb, color = 'r', marker = "x", s = 1)   
axs[1].scatter(UTC_time, (GRASS_rv) - RV_list_no_cb, color = 'b', marker = "x", s = 3)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals") 
plt.savefig("rm_and_residuals_no_cb.png")
plt.show()

#rm curve 
fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axs[0].scatter(UTC_time, rv_obs, color = 'k', marker = "x", s = 15, label = "EXPRES RVs") 
axs[0].plot(UTC_time, RV_list_cb , color = 'r', label = "Model - CB")
axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[0].set_xlabel("Time (UTC)")
axs[0].set_ylabel("RV [m/s]")
rms_model_cb = round(np.sqrt((np.nansum((rv_obs - RV_list_cb)**2))/len(rv_obs - RV_list_cb)),2)
axs[0].text(UTC_time[-40], -400, "Model RMS {}".format(rms_model_cb))
axs[0].legend()
#residuals 
axs[1].scatter(UTC_time, (rv_obs) - RV_list_cb, color = 'r', marker = "x", s = 1)  
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
axs[1].set_xlabel("Time (UTC)")
axs[1].set_ylabel("Residuals") 
plt.savefig("rm_and_residuals_cb.png")
plt.show()

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter(UTC_time, intensity_list)  
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity") 
# plt.savefig("intensity.png")
# plt.show()