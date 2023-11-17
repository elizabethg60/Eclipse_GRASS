import h5py
import matplotlib.pyplot as plt

serial = h5py.File("grid_serial.jld2", "r")
serial_N = serial["N_steps"][()]
serial_time = serial["time_vector"][()]

pa = h5py.File("grid_parallel.jld2", "r")
pa_N = pa["N_steps"][()]
pa_time = pa["time_vector"][()]
  
plt.scatter(serial_N, serial_time, label = "Serial Code")
plt.scatter(pa_N, pa_time, label = "Parallel Code")
plt.xlabel("Grid Size - N")
plt.ylabel("Time (seconds)")
plt.legend()
plt.savefig("gridvstime.png")
plt.show()

parallel_v2 = h5py.File("scaling_v2.jld2", "r")
num_workers_all = parallel_v2["num_workers_all"][()] 
wall_time = parallel_v2["wall_time"][()] 
wall_time_mid = parallel_v2["wall_time_mid"][()] 
wall_time_weak = parallel_v2["wall_time_weak"][()] 

plt.scatter(range(1,num_workers_all+1), wall_time, label = "Strong Scaling, 50x100 grid")
plt.scatter(range(1,num_workers_all+1), wall_time_mid, label = "Strong Scaling, 250x500 grid")
plt.scatter(range(1,num_workers_all+1), wall_time_weak, label = "Weaking Scaling")
plt.xlabel("number of workers")
plt.ylabel("Time (seconds)")
plt.legend()
plt.savefig("scaling_v2.png")
plt.show()