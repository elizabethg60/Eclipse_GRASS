from datetime import datetime, timedelta
import numpy as np

f=open("Reiners_Data.txt","r")
lines=f.readlines()[1:]

UTC_time = []
#airmass = []
fine_sample_UTC_time = []
for x in range(0,len(lines)):
   #UTC_time.append("2015-03-20T{}".format((lines[x].split()[1])))
   UTC_time.append((datetime.strptime("2015-03-20 {}".format((lines[x].split()[1])), "%Y-%m-%d %H:%M:%S.%f") + timedelta(seconds=float(50))).strftime("%Y-%m-%dT%H:%M:%S.%f"))
   #airmass.append((lines[x].split()[2]))

   inner = []
   for sec in np.arange(-30,30, 5):
      inner.append((datetime.strptime("2015-03-20 {}".format((lines[x].split()[1])), "%Y-%m-%d %H:%M:%S.%f") + timedelta(seconds=float(50)) + timedelta(seconds=float(sec))).strftime("%Y-%m-%dT%H:%M:%S.%f"))
   fine_sample_UTC_time.append(inner)
f.close()
# print(UTC_time)
# print(airmass)
print(fine_sample_UTC_time)

#print(list(np.linspace(0.4, 0.1, num = len(airmass))))