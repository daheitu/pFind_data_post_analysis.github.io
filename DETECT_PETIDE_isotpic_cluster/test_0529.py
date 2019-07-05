import os
from detect_isotopic import detectIsotopic

os.chdir(r"E:\Script\test_ZMQ_20190529")
file_list = os.listdir(os.getcwd())
for name in file_list:                
    if name[-4:] == ".ms1": 
        all_inf = open(name, 'r').readlines()
    else:
        continue
read_dict = {}

i = 0
while i < len(all_inf):
    mz = float(all_inf[i].split(" ")[0].strip())
    intensity = float(all_inf[i].split(" ")[1].strip())
    read_dict[mz] = intensity
    i += 1

found_isotopic = detectIsotopic(read_dict)

print(found_isotopic)