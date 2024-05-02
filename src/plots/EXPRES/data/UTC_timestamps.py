import pandas as pd
from datetime import datetime

#print out timestamps for EXPRES data 
data = pd.read_csv("EXPRES_data.csv")
UTC_time = []
for i in data["tobs"]:
    UTC_time.append(datetime.strptime(i, "%Y-%m-%d %H:%M:%S").strftime("%Y-%m-%dT%H:%M:%S"))
print(UTC_time)