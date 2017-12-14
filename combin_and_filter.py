import os
import re


def combine_and_filter_spectra_file(path):
    os.system("cd " + path)
    os.chdir(path)
    comb = open("combine.spectra", "w")
    comb.write(template_2[0])
    comb.write(template_2[1])
    print os.listdir(path)
    file_dic = {}
    for name in os.listdir(path):
        if name != "combine.spectra":
            file_dic[name] = open(name, 'r').readlines()
            for i in range(2, len(file_dic[name]) - 1, 2):
                if float(
                        file_dic[name][i].rstrip('\n').split('\t')[5]) < 0.001:
                    comb.write(file_dic[name][i])
                    comb.write(file_dic[name][i + 1])
    comb.close()
    return


template = open(r"/nibs/home/ycao/Collaboration/TC/YJ/quant/pQuant_cfg.txt",
                'r').readlines()

template_2 = open(
    r"/nibs/home/ycao/Collaboration/TC/TC_diUB_20161101/20161107/BS3/pQuant/cy_combine.spectra.xls_filter.txt",
    'r').readlines()

for i in range(18):
    path = "/nibs/home/ycao/Collaboration/TC/YJ/quant/YJ_" + str(
        i + 1) + "QUANT"
    combine_and_filter_spectra_file(path)

    tmp = open("pQuant.cfg", 'w')
    for line in template:
        if line[:16] == "PATH_INI_RESIDUE":
            if i % 3 == 0:
                tmp.write(line[:63] + "aa-BS3.ini;VIP")
                tmp.write("\n")
            elif i % 3 == 1:
                tmp.write(line[:63] + "aa-BS2G.ini;VIP")
                tmp.write("\n")
            else:
                tmp.write(line[:63] + "aa-DST.ini;VIP")
                tmp.write("\n")
        elif line[:24] == "PATH_IDENTIFICATION_FILE":
            tmp.write(line[:25] + path + "/combine.spectra")
            tmp.write("\n")
        elif line[:11] == "DIR_EXPORT=":
            tmp.write(line[:11] + path)
            tmp.write("\n")
        else:
            tmp.write(line)
    tmp.close()

    os.system(
        "java -jar /nibs/home/yhding/program/pQuant20160120/pQuant_64bit.jar ./pQuant.cfg"
    )
