import pandas as pd
import numpy as np
from astropy.time import Time
from datetime import datetime, timedelta

#print out timestamps for NEID data
data = pd.read_csv("neid_april_data.csv")
UTC_time = []
fine_sample_UTC_time = []
fine_sample_UTC_time_jd = []
for i in data["obsdate"][100:-10]:
    UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S"))

    inner = []
    inner_jd = []
    for sec in np.arange(0, 55, 5):
        inner.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=float(sec))).strftime("%Y-%m-%dT%H:%M:%S"))
        dt = datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=float(sec))
        inner_jd.append((Time(dt)).jd)
    fine_sample_UTC_time.append(inner)
    fine_sample_UTC_time_jd.append(inner_jd)
# UTC_time.remove('2024-04-08T20:08:41')
# UTC_time.remove('2024-04-08T20:05:55')
# UTC_time.remove('2024-04-08T19:56:16')
# UTC_time.remove('2024-04-08T19:50:46')
# print(UTC_time)

# print(fine_sample_UTC_time)
print(fine_sample_UTC_time_jd)
