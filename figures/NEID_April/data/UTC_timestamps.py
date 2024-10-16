import pandas as pd
import numpy as np
from astropy.time import Time
from datetime import datetime, timedelta

#print out timestamps for NEID data
data = pd.read_csv("neid_april_data.csv")
UTC_time = []
for i in data["obsdate"][100:-10]:
    UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S"))

print(UTC_time)

