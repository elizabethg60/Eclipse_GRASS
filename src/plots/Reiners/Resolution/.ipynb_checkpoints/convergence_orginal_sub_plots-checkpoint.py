import h5py
import numpy as np
import matplotlib.pyplot as plt 

# min_bound = 4
# max_bound = 20

# lats_array = []
# mean_vel = []
# for i in range(min_bound, max_bound + 1, 2):
#     lats_array.append(i)
    
#     file = h5py.File("model_data_{}.jld2".format(i), "r")
#     vel_no_cb = file["vel_no_cb"][()]

#     mean_vel.append(np.mean(file[vel_no_cb[0]][()]))

# plt.scatter(lats_array, mean_vel, color = 'b')
# #plt.axhline(y = 0, color = 'r', linestyle = '-')
# plt.xlabel("number of subgrid latitude slices")
# plt.ylabel("mean projected velocity")
# plt.savefig("single_time_convergence_orginal_sub.png")
# plt.show()
# plt.clf()


path = "/storage/home/efg5335/work/Eclipse_GRASS/src/plots/Reiners/"
#data 
f=open(path + "Reiners_Data.txt","r")
lines=f.readlines()[1:]
raw_rv = []
for x in lines[46:-35]:
    raw_rv.append(float(x.split()[4]))
f.close()

min_bound = 4
max_bound = 16

rms_total = []
rms = []
for lat in range(min_bound, max_bound + 1):
    #model
    file = h5py.File("convergence_orginal_sub/model_data_{}.jld2".format(lat), "r")
    RV_list_no_cb = file["RV_list_no_cb"][()][46:-35]

    raw_rv = np.array(raw_rv)
    raw_rv -= raw_rv[-1]

    RV_list_no_cb = np.array(RV_list_no_cb) 
    RV_list_no_cb -= RV_list_no_cb[-1]

    rms.append(round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2))
    rms_total.append(round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2))

plt.scatter(range(min_bound, max_bound + 1), rms, color = 'b')

min_bound = 18
max_bound = 20

rms = []
for lat in range(min_bound, max_bound + 1):
    #model
    file = h5py.File("convergence_orginal_sub/model_data_{}.jld2".format(lat), "r")
    RV_list_no_cb = file["RV_list_no_cb"][()][46:-35]

    raw_rv = np.array(raw_rv)
    raw_rv -= raw_rv[-1]

    RV_list_no_cb = np.array(RV_list_no_cb) 
    RV_list_no_cb -= RV_list_no_cb[-1]

    rms.append(round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2))
    rms_total.append(round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2))

plt.scatter(range(min_bound, max_bound + 1), rms, color = 'b')
    
tol = 10**(-16)
converge = False 
for i in range(1,len(rms_total)):
    if rms_total[i]-rms_total[i-1] <= tol:
        converge = True
        converge_rms = rms_total[i]
    
plt.axhline(y = converge_rms, color = 'r', linestyle = '-')
plt.axvline(x = 19, color = 'r', linestyle = '-')

plt.xlabel("number of subgrid latitude slices")
plt.ylabel("rms")
plt.savefig("eclipse_convergence_original_sub.png")
plt.show()
plt.close()