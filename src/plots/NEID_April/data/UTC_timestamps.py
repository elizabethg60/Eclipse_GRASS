import pandas as pd
from datetime import datetime, timedelta

#print out timestamps for NEID data
data = pd.read_csv("neid_april_data.csv")
UTC_time = []
for i in data["obsdate"][100:-10]:
    UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S"))
UTC_time.remove('2024-04-08T20:08:41')
UTC_time.remove('2024-04-08T20:05:55')
UTC_time.remove('2024-04-08T19:56:16')
UTC_time.remove('2024-04-08T19:50:46')
print(UTC_time)