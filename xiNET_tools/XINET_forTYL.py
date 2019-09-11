import os
import re

os.chdir(r"C:\Users\Administrator\Documents\pLink\pLink_task_2019.02.14.16.56.36\reports")
input_file = "RTT105_RPA_2019.02.14_Trypsin_EDC-DEv3.txt"

def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    if protein1 != protein2:
        link_type = "Inter"
    else:
        if abs(int(position1)-int(position2)) < 5:
            link_type = "Inter"
        else:
            link_type = "Intra" 
    return protein1, protein2, position1, position2, link_type


def main():
    file_list = os.listdir(os.getcwd())
    for file in file_list:
        if input_file == file:
            f = open(file).readlines()
            rep_name = file[:-4]+"_xinet.csv"
            b = open(rep_name, "w")
            b.write(",".join(["Score","SVM", "Protein1", "Protein2", "LinkPos1", "LinkPos2"]))
            b.write("\n")

            for line in f[1:]:
                line_list = line.rstrip("\n").split("\t")
                site = line_list[0]
                spectra = line_list[1]
                best_svm = line_list[3]
                write_list = [spectra, best_svm]
                if site.isdigit():
                    print(line)
                    continue
                else:                
                    if "/" in site:
                        continue
                    else:
                        if "Molecular" in site:
                            continue
                        else:
                            write_list.extend(get_linked_site_inform(site)[:-1])
                            b.write(",".join(write_list))
                            print(write_list)
                            b.write("\n")
            b.close()
            print("Done")
        else:
            continue


if __name__ == "__main__":
    main()
