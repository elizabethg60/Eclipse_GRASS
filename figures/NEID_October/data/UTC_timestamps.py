import pandas as pd
import numpy as np
from datetime import datetime, timedelta

#print out timestamps for NEID data
data = pd.read_csv("NEID_Data.csv")
UTC_time = []
for i in data["obsdate"]:
    UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S:%f"))

print(UTC_time)