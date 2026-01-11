from sunpy.net import Fido, attrs as a

# Search HMI continuum image at a given time
result = Fido.search(a.Time("2023-10-14 15:26:45", "2023-10-14 19:03:06"),
                     a.jsoc.Series("hmi.Ic_45s"), a.jsoc.Segment("continuum"), a.jsoc.Notify("efg5335@psu.edu"))

# Download
files = Fido.fetch(result)
print(files)
