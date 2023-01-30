# coding : utf-8

import os

wk_dir = r"G:\acetylation_R\urea\AC\VALA"
pep = "GVAAAKAAAAR"

# 是否已经导出plabel文件
def is_plabeled(wk_dir, mgf_name):
    plabel_name = mgf_name[:-4] + ".regular.plabel"
    plbl_path = os.path.join(wk_dir, plabel_name)
    if os.path.exists(plabel_name):
        return True, plbl_path
    else:
        return False, plbl_path

def write_pable(b, mgf_path):
    f = open(mgf_path).readlines()
    b.write("[FilePath]\n")
    b.write("File_Path=%s\n" % mgf_path)
    b.write("[Modification]\n")
    b.write("1=Delta_H(4)C(2)[AnyN-term]\n")
    b.write("2=Delta_H(4)C(2)[K]\n")
    b.write("3=Amidated[AnyC-term]\n")
    b.write("[xlink]\nxlink=NULL\n[Total]\n")
    write_lines = []
    spec_num = 0
    for line in f:
        if line.startswith("TITLE="):
            title = line.rstrip().split("=")[1]
            spec_num += 1
            write_lines.append("[Spectrum%d]\n" % spec_num)
            write_lines.append("name=%s\n" % (title[:-4]+".DTA"))
            write_lines.append("pep1=0 %s 1 0,1 11,2 12,3\n" % pep)
    write_lines.insert(0, "total=%d\n" % spec_num)
    b.writelines(write_lines)


def main(wk_dir):
    for fl in os.listdir(wk_dir):
        if fl.endswith(".mgf"):
            is_plabel, pl_path = is_plabeled(wk_dir, fl)
            if not is_plabel:
                b = open(pl_path, 'w')
                mgf_path = os.path.join(wk_dir, fl)
                write_pable(b, mgf_path)
                b.close()

main(wk_dir)                
