import os
import re
os.chdir(r"D:\E\Collabaration\TC\DIUB_5\Quant\B\results")
file_list = os.listdir(os.getcwd())


def site_correct(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    if "/" in linked_site:
        return linked_site
    else:
        if int(position1) <= int(position2):
            return linked_site
        else:
            a = linked_site.split("-")[0]
            b = linked_site.split("-")[1]
            return b + "-" + a


tmp = open('tmp', 'w')
for file in file_list:
    if file[-4:] ==".xls":
        f = open(file, 'r').readlines()
        for line in f:
            line_list = line.strip().split("\t")
            line_list_length = len(line_list)
            if line_list_length == 11:
                if line_list[2] != "Score" and float(line_list[2]) < 0.001:
                    line_list[-1] = line_list[-1]
                    del line_list[0]
                    tmp.write("\t".join(line_list))
                    tmp.write("\n")
tmp.close()

t = open("tmp", "r").readlines()
raw_list = []
link_pair_list = []
for line in t:
    tmp_line_list = line.strip().split("\t")
    title = tmp_line_list[0]
    raw = title[:title.find(".")]
    link_pair = tmp_line_list[-1]
    if raw not in raw_list:
        raw_list.append(raw)
    if link_pair not in link_pair_list:
        link_pair_list.append(link_pair)
raw_list.sort()
print(raw_list, link_pair_list)

r = open("report.txt", 'w')
column =["Linked_site", "totel spectra", "Best evalue"]
for raw in raw_list:
        column.append(raw+"_specta")
        column.append(raw + "_Best evale")
r.write("\t".join(column))
r.write("\n")
for link_pair in link_pair_list:
    spectra = 0
    E_value_list = []
    raw_spectra_dic = {}
    raw_evalue_dic = {}
    for raw in raw_list:
        raw_spectra_dic[raw] = 0
        raw_evalue_dic[raw] = []
    for line in t:
        l_list = line.strip().split("\t")
        title = l_list[0]
        raw_name = title[:title.find(".")]
        E_value = float(l_list[1])
        if line.strip().split("\t")[-1] == link_pair:
            raw_spectra_dic[raw_name] += 1
            raw_evalue_dic[raw_name].append(E_value)
    for line in t:
        if line.strip().split("\t")[-1] == link_pair:
            spectra += 1
            E_value_list.append(float(line.strip().split("\t")[1]))

    l_f_w = [link_pair, str(spectra), str(min(E_value_list))]
    for raw in raw_list:
        if raw_spectra_dic[raw] != 0:
            l_f_w.append(str(raw_spectra_dic[raw]))
            l_f_w.append(str(min(raw_evalue_dic[raw])))
        else:
            l_f_w.append("")
            l_f_w.append("")

    r.write("\t".join(l_f_w))
    r.write("\n")
r.close()