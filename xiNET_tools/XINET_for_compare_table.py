import os
import re

wk_path = r"D:\wechat_files\WeChat Files\yunzaitianya2010\FileStorage\File\2021-11" # 比较文件所在的路径

input_file = "YangSS_RPAvsRPA-105_30NT.csv" #比较文件的命字


####################Don't change the following lines#########################
os.chdir(wk_path)

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


def get_delta_spec(spec1, spec2):
    if "" in [spec1, spec2]:
        other = float([x for x in [spec1, spec2] if x != ""][0])
        return str(round(other, 4))
    else:
        return str(round(abs(float(spec1) - float(spec2)), 4))


def main():
    file_list = os.listdir(os.getcwd())
    for file in file_list:
        if input_file == file:
            f = open(file).readlines()
            rep_name_up = file[:-4]+"_up_xinet.csv"
            rep_name_down = file[:-4]+"_down_xinet.csv"
            b = open(rep_name_up, "w")
            c = open(rep_name_down, 'w')
            title_line = "Score, Protein1, Protein2, LinkPos1, LinkPos2\n"
            b.write(title_line)
            c.write(title_line)
            for line in f[1:-1]:
                line_list = line.rstrip("\n").split(",")
                site = line_list[0]
                spectra1 = line_list[1]
                spectra2 = line_list[4]
                spectra = get_delta_spec(spectra1, spectra2)
                write_list = [spectra]
                up_or_d = line_list[-1]

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
                            site_extrac = get_linked_site_inform(site)
                            if site_extrac[-1] == "Inter" or site_extrac[-1] == "Intra" :
                                write_list.extend(site_extrac[:-1])
                                if up_or_d == "up":
                                    b.write(",".join(write_list)+"\n")
                                    print(write_list)
                                elif up_or_d == "down":
                                    c.write(",".join(write_list)+"\n")
                                    print(write_list)
            b.close()
            c.close()
            print("Done")
        else:
            continue


if __name__ == "__main__":
    main()
