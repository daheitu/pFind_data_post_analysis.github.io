import os
os.chdir(r"E:\workspace\pFindTask97\result")
#pro_looking = "sp|P35527|K1C9_HUMAN"  # important! enter the protein name you want search
modify = "hp3m[C]"  # important! enter the modification type you want search


def filter_modification_from_pfind(path):
    os.chdir(path)
    os.chdir(os.getcwd())
    tab = open("pFind.protein", 'r').readlines()
    b = open("report.txt", 'w')
    b.write("\t".join([
        "ID", "Sequence", "Protein", "Modification", "Modi_site",
        "Spectra Number", "Final score", "Spectra title"
    ]))
    b.write("\n")
    pms_total = 0
    for line in tab[2:]:
        line_list = line.rstrip("\n").split("\t")
        if len(line_list) < 11:
            continue
        else:
            mod_tab = line_list[8]
            psm = line_list[-1]
            if modify in mod_tab:
                pms_total += int(psm)
            else:
                continue
    return pms_total


def main():
    path = os.getcwd()
    print(filter_modification_from_pfind(path))


if __name__ == "__main__":
    main()
