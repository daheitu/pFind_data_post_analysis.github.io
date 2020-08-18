# coding = utf-8
# author Yong CAO

import os
import sys
from numpy import median
sys.path.append(r"C:\Users\Yong Cao\Documents\github\pFind_data_post_analysis.github.io\dsso_ms3_ions")
from readMS3info import getLinageInfo, readms3Info

os.chdir(r"F:\data_from_paper\maxlinker_MCP\JinLiang10366422")
raw_name = "JinLiang10366422_XL_SCX_v3_Fr_18"
# ms3_fl = raw_name + ".ms3"
ms3InfoDic = readms3Info("./JinLiang10366422_XL_SCX_v3_Fr_18.ms3")
ms2TriggerMS3dic = getLinageInfo(ms3InfoDic)
# print(ms2TriggerMS3dic)

title_num_dic = {}

f = open(r"F:\data_from_paper\maxlinker_MCP\JinLiang10366422\raw\JinLiang10366422_XL_SCX_v3_Fr_18_CIDFT_output\tmps\doublets.txt", 'r').readlines()
i = 0
while i < len(f):
    if not f[i].startswith(raw_name):
        i += 1
        # print("line %d is wrong" % i)
    else:
        title = f[i].strip()
        scan_num = int(title.split('.')[1])
        p = i + 1
        if scan_num in list(ms2TriggerMS3dic.keys()):
            delta_32_2plus = 0
            paired_num = 0
            while p < len(f):
                if f[p].startswith("#masspairs"):
                    
                    break
                else:
                    sub_list = f[p].strip().split("\t")
                    if sub_list[0] == "all_doublet":
                        charge = int(sub_list[4])
                        if charge > 1:
                            delta_32_2plus += 1
                    if sub_list[0] == "found":
                        found_type = sub_list[3]
                        if found_type == "1":
                            print(sub_list)

                            paired_num += 1
                    p += 1
            
            
            if paired_num > 0:
                # num_paired += 1
                title_num_dic[title] = [delta_32_2plus, True]
            else:
                title_num_dic[title] = [delta_32_2plus, False]
        i = p + 1

info_list = list(title_num_dic.values())
all_num_list = [x[0] for x in info_list]
true_num_list = [x[0] for x in info_list if x[1]]
print(len(all_num_list), median(all_num_list))
print(len(true_num_list), median(true_num_list))

# print(median(list(title_num_dic.values())))      
# for scan in ms2TriggerMS3dic.keys():
