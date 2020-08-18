import os

rep_path = r"G:\msData\SLu"

para_demo = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2020.08.11.08.13.10\pQuant_cfg.txt"

plink_bin_path = r"E:\pFindStudio\pLink2.3.9_0415\bin"

# rp_path = r"G:\msData\SLu\G1\LuS_NP_output_BS3\reports"


def change_file_format(reporter_path):
    for fl in os.listdir(reporter_path):
        if fl.endswith("cross-linked_spectra.csv"):
            f = open(os.path.join(reporter_path, fl)).readlines()
            b = open(os.path.join(reporter_path, fl[:-4] + "2.csv"), 'w')
            b.write(f[0])
            for line in f[1:]:
                linelist = line.split(',')
                linelist[14] = linelist[13]
                b.write(','.join(linelist))
            b.close()

# change_file_format(reporter_path)


def get_pf1_path(reporter_path):
    pf_path_list = []
    pathpf = "\\".join(reporter_path.split("\\")[:-2])
    print(pathpf)
    for fl in os.listdir(pathpf):
        if fl.endswith(".pf1"):
            print(fl)
            pf_path_list.append(os.path.join(pathpf, fl))
    print(pf_path_list)
    pf_path_list.append(";")
    return pf_path_list
# get_pf1_path(reporter_path)

def get_xl_path(reporter_path):
    for fl in os.listdir(reporter_path):
        if fl.endswith("cross-linked_spectra2.csv"):
            print(os.path.join(reporter_path, fl))
            return os.path.join(reporter_path, fl)


def write_pQuant_para(reporter_path):
    f = open(para_demo).readlines()
    pf1_path_list = get_pf1_path(reporter_path)
    id_file_path = get_xl_path(reporter_path)
    pq_path = os.path.join(reporter_path, "pQuant_xlink_d4.cfg")
    b = open(pq_path, 'w')
    for line in f:
        if line.startswith("PATH_MS1"):
            b.write("PATH_MS1=" + "|".join(pf1_path_list)+"\n")
        elif line.startswith("PATH_IDENTIFICATION_FILE"):
            b.write("PATH_IDENTIFICATION_FILE=" + id_file_path + "|;\n")
        elif line.startswith("DIR_EXPORT"):
            b.write("DIR_EXPORT="+reporter_path + "\\pQuant_xlink;\n")
        else:
            b.write(line)
    b.close()
    os.chdir(plink_bin_path)
    os.system("pQuant.exe " + pq_path)



def main():
    for root, dirs, fls in os.walk(rep_path):
        if root.split("\\")[-1] == 'reports' and "heavy" in root:
            change_file_format(root)
            write_pQuant_para(root)


main()