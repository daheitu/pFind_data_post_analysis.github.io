import os
import re
os.chdir(
    r"C:\Users\Yong\Desktop\sgc\sheyang_list"
)
spec_cutoff = 2  # spectra number cut-off
E_value_cutoff = 0.01


def get_report_file_name():
    path = os.getcwd()
    path_d = os.path.dirname(path)
    os.chdir(path_d)
    file_list = os.listdir(path_d)
    for fl in file_list:
        if fl[-5:] == "plink":
            para = open(fl).readlines()
        else:
            continue
    for line in para:
        if line[:10] == "spec_title":
            spec_title = line.rstrip("\n").split("=")[1].strip()
        elif line[:11] == "enzyme_name":
            enzyme = line.rstrip("\n").split("=")[1].strip()
        elif line[:7] == "linker1":
            linker = line.rstrip("\n").split("=")[1].strip()
        else:
            continue
    report_file_name = spec_title + "_" + enzyme + "_" + linker + ".txt"
    os.chdir(path)
    return report_file_name


# print(get_report_file_name())


def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    if position1 == position2:
        return "Self cross-link"
    else:
        if int(position1) < 500 and int(position2) > 900:
            link_type = "Inter"
        else:
            link_type = "Intra"
        return link_type


f = open("list.txt", 'r').readlines()
b = open("repo.csv", 'w')
for line in f[0:]:
    b.write(",".join([line.rstrip("\n"), get_linked_site_inform(line.rstrip("\n"))]))
    b.write("\n")

b.close()
print("done")
