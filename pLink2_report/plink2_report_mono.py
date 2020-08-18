# coding: utf-8
"""
This script can help you to summary the plink2 report file
"""

import os
import re

wk_dir = r"G:\msData\SLu\G0"
# reports_path = r"G:\msData\SLu\G9\LuS_NP_output_BS3\reports"

spec_cutoff = 2  # spectra number cut-off
Best_evalue_cutoff = 2 # 交联位点对层次最好的e-value cutoff
E_value_cutoff_SpecLvl = 2 # 谱图层次的e-value cutoff


def find_monoPeptides_File(reports_path):
    for fl in os.listdir(reports_path):
        if fl.endswith("mono-linked_sites.csv"):
            return os.path.join(reports_path, fl)
    return ""

def main(reports_path):
    b = open(os.path.join(reports_path, "Mono_link.csv"), 'w')
    b.write("site,spec,bestE-value,bestSVM\n")
    mono_fl = find_monoPeptides_File(reports_path)
    print(mono_fl)
    f = open(mono_fl).readlines()
    i = 2
    while i < len(f):
        linelist = f[i].split(',')
        if not linelist[0].isdigit():
            print("wrong")
        else:
            site = linelist[1]
            spec = linelist[3].strip()
            p = i + 1
            bestE = 2
            pep_list = []
            while p < len(f):
                sub_list = f[p].split(',')
                if sub_list[0].isdigit():
                    break
                else:
                    bestSVM = sub_list[9]
                    if float(sub_list[8]) < bestE:
                        bestE = float(sub_list[8])
                    pep = sub_list[5]
                    if pep not in pep_list:
                        pep_list.append(pep)
                    p += 1
            b.write(",".join([site, spec, str(bestE), bestSVM, ";".join(pep_list)])+"\n")
            i =  p
    b.close()


if __name__ == "__main__":
    for root, dirs, fls in os.walk(wk_dir):
        root_list = root.split("\\")
        if root_list[-1] == "reports":
            main(root)
# main(r"G:\msData\SLu\G9\LuS_NP_output_BS3\reports")