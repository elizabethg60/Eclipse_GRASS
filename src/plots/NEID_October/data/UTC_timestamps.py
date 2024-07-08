import pandas as pd
import numpy as np
from datetime import datetime, timedelta

#print out timestamps for NEID data
data = pd.read_csv("NEID_Data.csv")
UTC_time = []
for i in data["obsdate"][15:-150]:
    UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S"))

    # if i == data["obsdate"][15]:
    #     inner = []
    #     for sec in np.arange(0, 55, 5):
    #         inner.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=float(sec))).strftime("%Y-%m-%dT%H:%M:%S"))
    #     print(inner)
print(UTC_time)

# import matplotlib.pyplot as plt
# mpl = plt.matplotlib 
# import matplotlib.dates as mdates

# first_time_fine = ["2023-10-14T15:26:18", "2023-10-14T15:26:23", "2023-10-14T15:26:28", "2023-10-14T15:26:33", "2023-10-14T15:26:38", "2023-10-14T15:26:43", "2023-10-14T15:26:48", "2023-10-14T15:26:53", "2023-10-14T15:26:58", "2023-10-14T15:27:03", "2023-10-14T15:27:08"]
# first_time_intensity = [8.438065653169185e11, 8.433680132962926e11, 8.428748370503627e11, 8.424467423005983e11, 8.419269260639641e11, 8.414955491635146e11, 8.410779998526726e11, 8.405270565811194e11, 8.401363338521866e11, 8.395848489158354e11, 8.391811487038573e11]

# UTC_time = []
# for i in first_time_fine:
#     UTC_time.append((datetime.strptime(i, "%Y-%m-%dT%H:%M:%S")))

# #intensity
# fig = plt.figure()
# ax1 = fig.add_subplot()
# ax1.scatter(UTC_time, np.array(first_time_intensity)/max(first_time_intensity), label = "Model")  
# ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
# ax1.set_xlabel("Time (UTC)")
# ax1.set_ylabel("Relative Intensity") 
# #plt.legend()
# plt.savefig("intensity_single_exposure.png")
# plt.show()