import pandas as pd
import numpy as np
from datetime import datetime, timedelta

#print out timestamps for NEID data
data = pd.read_csv("NEID_Data.csv")
UTC_time = []
fine_sample_UTC_time = []
for i in data["obsdate"][15:-178]: #-150
#     UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S:%f"))

#     inner = []
#     for sec in np.arange(0.0, 56.0, 1.0):
#         inner.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=sec)).strftime("%Y-%m-%dT%H:%M:%S.%f"))
#     fine_sample_UTC_time.append(inner)

    inner = []
    for sec in [0.0, 15.0, 30.0, 45.0, 55.0, 70.0]:
        inner.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=sec)).strftime("%Y-%m-%dT%H:%M:%S.%f"))
    fine_sample_UTC_time.append(inner)

# print(UTC_time)

# print(len(data["obsdate"][15:-178]) * len(np.arange(0.0, 56.0, 1.0)))
# print(len(sum(fine_sample_UTC_time, [])))
# df = pd.DataFrame(sum(fine_sample_UTC_time, []), columns=['FineTime'])

# # Save the DataFrame to CSV
# df.to_csv("NEID_October_fine_time.csv", index=False)

print(data["obsdate"][15])
print(sum(fine_sample_UTC_time, [])[0:7])
print(data["obsdate"][16])
df = pd.DataFrame(sum(fine_sample_UTC_time, []), columns=['FineTime'])

# Save the DataFrame to CSV
df.to_csv("NEID_October_15s_time.csv", index=False)