import numpy as np
from datetime import datetime, timedelta

#print timestamps for binned data
data_time = np.loadtxt("Boulder_Data_bin.txt")[:, 0]
UTC_time = []
for i in data_time:
    UTC_time.append((datetime.fromtimestamp(i) + timedelta(hours=4)).strftime("%Y-%m-%dT%H:%M:%S"))
print(UTC_time)