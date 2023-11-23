f=open("Reiners_Data.txt","r")
lines=f.readlines()[1:]
UTC_time = []
for x in lines:
   UTC_time.append("2015-03-20T{}".format((x.split()[1])))
f.close()
print(UTC_time)
