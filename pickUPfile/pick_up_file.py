import os
import re

for i in range(18):
    os.system("mkdir -p YJ_" + str(i + 1) + "QUANT")

for word in ["BS3", "BS2G", "DST"]:
    os.system("cd /nibs/home/ycao/Collaboration/TC/YJ/MGF/" + word +
              "/2.report/sample1")
    os.chdir("/nibs/home/ycao/Collaboration/TC/YJ/MGF/" + word +
             "/2.report/sample1")
    os.system("rename cy_ cy_0 cy_[0-9].spectra.xls ")
    filename_list = os.listdir(os.getcwd())
    for name in filename_list:
        if re.match('cy_[0-1][0-9].spectra.xls', name) is not None:
            tmp = open(name, 'r').readlines()
            title = tmp[2].strip().split("\t")[1]
            if "_R2" in title:
                raw_name = title[:title.find('_R2')]
            else:
                raw_name = title[:title.find('.')]
            os.system(
                "cp " + name + " /nibs/home/ycao/Collaboration/TC/YJ/quant/" +
                raw_name + "QUANT/" + name)