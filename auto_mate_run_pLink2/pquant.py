import os

plink_bin_path = r"E:\pFindStudio\pLink2.3.9_0415\bin"
raw_p = r"G:\msData\20200419\BSA"

def add_plink_path(line):
    header, file = line.split("=")
    file_path = os.path.join(plink_bin_path, file)
    return "=".join([header, file_path])


def get_ms1_files(raw_path):
    ms1_fl_list = []
    for fl in os.listdir(raw_path):
        if fl.endswith("pf1"):
            ms1_fl_list.append(os.path.join(raw_path, fl))
    return "|".join(ms1_fl_list)


def get_id_file():
    for root, drs, fls in os.walk(raw_p):
        for fl in fls:
            if fl.endswith("cross-linkedsites.csv"):
                return os.path.join(root, fl)+'|'


os.chdir(plink_bin_path)
f = open('./pQuant_cfg.txt').readlines()
b = open('pquant_para.txt', 'w')
for line in f:
    header = line.split('=')[0]
    if header in ["PATH_INI_ELEMENT", 'PATH_INI_MODIFICATION', 'PATH_INI_RESIDUE', 'PATH_INI_GLYCO', 'PATH_INI_LINKER']:
    # if line.startswith("PATH_INI_ELEMENT"):
        b.write(add_plink_path(line))
    # elif line.startswith('PATH_INI_MODIFICATION'):
    elif header == "PATH_MS1":
        b.write('='.join([header, get_ms1_files(raw_path)])+';\n')
    elif header == "EXTENSION_TEXT_MS1":
        b.write("=".join([header, 'pf1'])+';\n')
    elif header == "PATH_IDENTIFICATION_FILE":
        b.write('='.join([header, get_id_file()])+';\n')
    elif header == "NUMBER_SCANS_HALF_CMTG":
        b.write('='.join([header, "1000"])+';\n')
    elif header == "DIR_EXPORT":
        export_path = "\\".join(get_id_file().split('\\')[:-2])
        b.write('='.join([header, export_path])+';\n')
    else:
        b.write(line)
        
b.close()