import h5py
import numpy as np
import matplotlib.pyplot as plt 

min_bound = 80
max_bound = 219

# lats_array = []
# mean_vel = []
# for i in range(min_bound, max_bound + 1, 5):
#     lats_array.append(i)
    
#     file = h5py.File("convergence_orginal_large/model_data_{}.jld2".format(i), "r")
#     vel_no_cb = file["vel_no_cb"][()]

#     mean_vel.append(np.mean(file[vel_no_cb[0]][()]))

# plt.scatter(lats_array, mean_vel, color = 'b')
# plt.axhline(y = 0, color = 'r', linestyle = '-')
# plt.xlabel("number of latitude slices")
# plt.ylabel("mean projected velocity")
# plt.savefig("single_time_convergence_original_large.png")
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


rms = []
for lat in range(min_bound, max_bound + 1):
    #model
    file = h5py.File("convergence_orginal_large/model_data_{}.jld2".format(lat), "r")
    RV_list_no_cb = file["RV_list_no_cb"][()][46:-35]

    raw_rv = np.array(raw_rv)
    raw_rv -= raw_rv[-1]

    RV_list_no_cb = np.array(RV_list_no_cb) 
    RV_list_no_cb -= RV_list_no_cb[-1]

    rms.append(round(np.sqrt((np.nansum((raw_rv - RV_list_no_cb)**2))/len(raw_rv - RV_list_no_cb)),2))

tol = 10**(-16)
converge = False 
for i in range(1,len(rms)):
    if rms[i]-rms[i-1] <= tol:
        converge = True
        converge_rms = rms[i]
    
plt.axhline(y = converge_rms, color = 'r', linestyle = '-')
plt.axvline(x = 200, color = 'r', linestyle = '-')
plt.scatter(range(min_bound, max_bound + 1), rms, color = 'b')
plt.xlabel("number of latitude slices")
plt.ylabel("rms")
plt.savefig("eclipse_convergence_original_large.png")
plt.show()
plt.close()