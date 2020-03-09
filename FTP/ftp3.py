import urllib.request
import os

ftpIP = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2017/04/PXD006131"
response = urllib.request.urlopen(ftpIP)
a = response.read().decode('utf-8').split("\n")

for line in a:
    name = line.split(" ")[-1].strip()
    if "HCD" in name and name.endswith(".raw"):
        print(name)
        path = os.path.join(ftpIP, name)
        print(path)
        cmd = "wget " + path
        print(cmd)
        #os.system(cmd)
