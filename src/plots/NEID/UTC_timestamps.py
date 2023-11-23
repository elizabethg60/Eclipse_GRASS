import pandas as pd
from datetime import datetime, timedelta

data = pd.read_csv("NEID_Data.csv")
UTC_time = []
for i in data["obsdate"][15:-150]:
    UTC_time.append((datetime.strptime(i, "%Y-%m-%d %H:%M:%S") + timedelta(seconds=27.5)).strftime("%Y-%m-%dT%H:%M:%S"))
print(UTC_time)