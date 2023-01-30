# coding: utf-8

"""
此脚本用于将plabel文件中的部分按需求过滤出来
"""


import os

wk_dir = r"Z:\pLink_workspace\pLink_task_2022.10.09.13.31.22\BSA_1mM_BS3_R1_T1_HCDFT.cross-linked.DSS_pepNC.plabel"

f = open(wk_dir).readlines()


def is_tail_head(pep_info_list):
    p1 = pep_info_list[3]
    p2 = pep_info_list[5]
    s1 = pep_info_list[1]
    s2 = pep_info_list[2]
    if s1 == "1" or s2 == "1":
        if len(p1) == int(s1) or len(p2) == int(s2):
            return True
        else:
            return False
    else:
        return False


wlist = []
for line in f[:8]:
    print(line)
    if line.startswith("xlink"):
        wlist.append(line.replace("DSS_pepNC", "DSS_pepNC_label"))
    else:
        wlist.append(line)

i = 0
for i in range(9, len(f), 3):
    print(f[i])
    pep_info_list = f[i+2].split(" ")
    if not is_tail_head(pep_info_list):
        i += 1
        wlist.append("[Spectrum%d]\n" % i)
        wlist.append(f[i+1])
        wlist.append(f[i+2])
wlist.insert(8, "total=%d\n" % i)


