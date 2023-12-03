f=open("Reiners_Data.txt","r")
lines=f.readlines()[1:]
UTC_time = []
airmass = []
for x in lines:
   UTC_time.append("2015-03-20T{}".format((x.split()[1])))
   airmass.append((x.split()[2]))
f.close()
#print(UTC_time)
#print(airmass)

import numpy as np
print(list(np.linspace(0.4, 0.1, num = len(airmass))))