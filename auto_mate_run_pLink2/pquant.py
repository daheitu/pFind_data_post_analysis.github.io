import os
# from automate_pLink2 import makedir


def makedir(rootpath, dir_name):
    path = os.path.join(rootpath, dir_name)
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path)


def add_plink_path(line, plink_bin_path):
    header, file = line.split("=")
    file_name = file.split('\\')[-1]
    file_path = os.path.join(plink_bin_path, file_name)
    return "=".join([header, file_path])


def get_ms1_files(raw_path):
    ms1_fl_list = []
    for fl in os.listdir(raw_path):
        if fl.endswith("pf1"):
            ms1_fl_list.append(os.path.join(raw_path, fl))
    return "|".join(ms1_fl_list)


def get_id_file(raw_p):
    for root, drs, fls in os.walk(raw_p):
        if "score" in root:
            for fl in fls:
                if fl.endswith("cross-linked_spectra.csv"):
                    return os.path.join(root, fl)+'|'


def gene_run_pquant(plink_bin_path, raw_p):
    os.chdir(plink_bin_path)
    xl_spectra_path = get_id_file(raw_p)
    export_path = "\\".join(xl_spectra_path.split('\\')[:-2])
    print(export_path)
    f = open(os.path.join(plink_bin_path, 'pQuant_CXMS_demo.txt')).readlines()
    pquant_para_path = os.path.join(export_path, 'pquant_para.txt')
    b = open(pquant_para_path, 'w')
    for line in f:
        header = line.split('=')[0]
        if header in ["PATH_INI_ELEMENT", 'PATH_INI_MODIFICATION', 'PATH_INI_RESIDUE', 'PATH_INI_GLYCO', 'PATH_INI_LINKER']:
            b.write(add_plink_path(line, plink_bin_path))
        elif header == "PATH_MS1":
            b.write('='.join([header, get_ms1_files(raw_p)])+';\n')
        elif header == "EXTENSION_TEXT_MS1":
            b.write("=".join([header, 'pf1'])+';\n')
        elif header == "PATH_IDENTIFICATION_FILE":
            b.write('='.join([header, xl_spectra_path])+';\n')
        elif header == "NUMBER_SCANS_HALF_CMTG":
            b.write('='.join([header, "1000"])+';\n')
        elif header == "DIR_EXPORT":
            makedir(export_path, 'pQuant')
            out_path = os.path.join(export_path, 'pQuant')
            b.write('='.join([header, out_path])+';\n')
        else:
            b.write(line)
            
    b.close()
    pquant_cmd = 'pQuant.exe ' + pquant_para_path
    os.system(pquant_cmd)

if __name__ == "__main__":
    plink_bin_path = r"E:\pFindStudio\pLink2.3.9_0415\bin"
    raw_p = r"K:\20200820_PRM" # raw文件所在目录
    gene_run_pquant(plink_bin_path, raw_p)