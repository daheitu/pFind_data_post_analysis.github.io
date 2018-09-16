import os

os.chdir(r"D:\Program Files (x86)\pFind_build_20170814\GMDSD")  # path
file = r"pFind.protein"  # filename
delta_mass = "295.12"  # mass


report_file_name = "report_delta_mass_" + delta_mass + ".txt"
b = open(report_file_name, 'w')
f = open(file, "r").readlines()
b.write(f[1])

for line in f[3:]:
    if line[:3] == "---":
        break
    else:
        line_list = line.strip().split("\t")
        if len(line_list) == 25:
            modificaiton = line_list[6]
            check_str = "PFIND_DELTA_" + delta_mass
            if check_str in modificaiton:
                b.write(line)
            else:
                continue
        else:
            continue

print("Done")
