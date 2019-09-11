import urllib.request
import os

ftpIP = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/07/PXD012546"
response = urllib.request.urlopen(ftpIP)
a = response.read().decode('utf-8').split("\n")

for line in a:
    name = line.split(" ")[-1].strip()
    if name.startswith("R1") and name.endswith("raw"):
        print(name)
        path = os.path.join(ftpIP, name)
        print(path)
        cmd = "wget " + path
        if "A10" not in cmd:
            os.system(cmd)