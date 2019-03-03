# import os

mass_find = -1
delta_ppm = 100000


def generate_ion_mass_range(num):
    deta = num * delta_ppm / 1000000
    return num - deta, num + deta


f = open(r"D:\Program Files\pFind3\bin\modification.ini", 'r').readlines()
modifi_dic = {}
for line in f[1:]:
    if "name" not in line and "=" in line:
        name = line.strip().split("=")[0]
        name_info_list = line.strip().split("=")[1].split(" ")
        site = name_info_list[0]
        com_state = name_info_list[1]
        mass = name_info_list[2]
        neutral_loss = bool(name_info_list[4])
        chem_formal = name_info_list[-1]
        modifi_dic[name] = [site, com_state, mass, chem_formal]
    else:
        continue

mass_low, mass_up = generate_ion_mass_range(mass_find)

find_bool = False
for name in modifi_dic:
    if float(modifi_dic[name][2]) > mass_low and float(
            modifi_dic[name][2]) < mass_up:
        find_bool = True
        print(name, modifi_dic[name])
    else:
        continue
if find_bool is False:
    print("None")
else:
    pass