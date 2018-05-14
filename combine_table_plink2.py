import os
import re

spec_cutoff = 2
evaule_cutoff = 0.01
os.chdir(r"C:\Users\Yong\Desktop\results\plink2")

# file_list = os.listdir(os.getcwd())


def site_correct(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    if int(position1) <= int(position2):
        return linked_site
    else:
        a = linked_site.split("-")[0]
        b = linked_site.split("-")[1]
        return b + "-" + a


def site_list_correction(linked_site):
    if "/" in linked_site:
        site_list = linked_site.split("/")
        i = 0
        while i < len(site_list):
            if "REVERSE" in site_list[i]:
                site_list.pop(i)
                i -= 1
            elif "gi|CON" in site_list[i]:
                site_list.pop(i)
                i -= 1
            else:
                pass
            i += 1

        if len(site_list) == 0:
            return ""
        elif len(site_list) == 1:
            return site_correct(site_list[0])
        else:
            for i in range(len(site_list)):
                site_list[i] = site_correct(site_list[i])
            site_list.sort()
            return "/".join(site_list)
    else:
        if "gi|CON" in linked_site:
            return ""
        elif "REVERSE" in linked_site:
            return ""
        else:
            return site_correct(linked_site)


def combine_peptide_file(path):
    file_list = os.listdir(path)
    tmp = open('tmp', 'w')
    for file in file_list:
        if "cross-linked_spectra" in file:
            f = open(file, 'r').readlines()
            for line in f[1:]:
                line_list = line.strip().split(",")
                line_list[13] = site_list_correction(line_list[13][:-1])

                tmp.write("\t".join(line_list[:14]))
                tmp.write("\n")
        else:
            continue
    tmp.close()
    return


def handle_tmp(file):
    t = open(file, "r").readlines()
    pair_dic = {}
    for line in t:
        line_list = line.strip().split("\t")
        site = line_list[-1]
        title = line_list[1]
        evalue = line_list[9]
        svm = line_list[10]
        pep = line_list[4]
        if site not in pair_dic:
            pair_dic[site] = [
                site, 1, [float(svm)], [float(evalue)], [pep], [title]
            ]
        else:
            pair_dic[site][1] += 1
            pair_dic[site][2].append(float(svm))
            pair_dic[site][3].append(float(evalue))
            pair_dic[site][4].append(pep)
            pair_dic[site][5].append(title)

    for site in pair_dic:
        svm_list = pair_dic[site][2]
        evalue_list = pair_dic[site][3]
        pep_list = pair_dic[site][4]
        index = svm_list.index(min(svm_list))
        evalue = evalue_list[index]
        pep = pep_list[index]
        title = pair_dic[site][5][index]
        pair_dic[site][1] = str(pair_dic[site][1])
        pair_dic[site][2] = str(min(svm_list))
        pair_dic[site][3] = str(evalue)
        pair_dic[site][4] = pep
        pair_dic[site][5] = title

    return pair_dic


def main():
    combine_peptide_file(os.getcwd())
    pair_dic = handle_tmp("tmp")
    b = open("plink2.txt", "w")
    for site in pair_dic:
        if int(pair_dic[site][1]) > spec_cutoff and float(
                pair_dic[site][3]) < evaule_cutoff:
            b.write("\t".join(pair_dic[site]))
            b.write("\n")
        else:
            continue


if __name__ == "__main__":
    main()
