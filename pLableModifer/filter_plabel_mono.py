
from distutils.log import info
import os
plable_path = r"Z:\STY_PROJ\low_comp_sample\BS3_BSA_CY\xiSearch_results\bsa_BS3SmallScale_output\RES5"

fasta_path = r"Z:\STY_PROJ\low_comp_sample\bsa.fasta"


def get_fasta_dic(fasta_path):
    fasta_dic = {}
    f = open(fasta_path).readlines()
    i = 0
    while i < len(f):
        if not f[i].startswith(">"):
            i += 1
        else:
            name = f[i][1:-1]
            seq = ""
            p = i+ 1
            while p < len(f):
                if f[p].startswith(">"):
                    break
                else:
                    seq += f[p].strip()
                    p += 1
            fasta_dic[name] = seq
            i = p
    print(len(fasta_dic))
    return fasta_dic 


def is_mono_link_pep(p1,p2, fasta_dic):
    is_mono_list = []
    for pro in fasta_dic:
        seq = fasta_dic[pro]
        if p1 in seq and p2 in seq:
            if p1+p2 in seq or p2+p1 in seq:#is_poss_mono(p1, p2 ,seq) or is_poss_mono(p2, p1, seq):
                is_mono_list.append("T")
            else:
                is_mono_list.append("F")
    if is_mono_list == []:
        return False
    elif is_mono_list.count("T") >= 1:#is_mono_list.count("F"):
        # print(is_mono_list)
        return True
    else:
        return False


def not_kk_links(pep1, pep2, site1, site2):
    aa1 = pep1[int(site1) - 1]
    aa2 = pep2[int(site2) - 1]
    if aa1 == "K" and aa2 == "K":
        return False
    else:
        return True


fasta_dic = get_fasta_dic(fasta_path)
for fl in os.listdir(plable_path):
    if fl.endswith(".plabel"):
        print(fl)
        f = open(os.path.join(plable_path, fl)).readlines()
        tgt_path = os.path.join(plable_path, fl[:-7]+"filter_spe.plabel")
        b = open(tgt_path, 'w')
        for i in range(20):
            if f[i].startswith("[Total]"):
                b.write(f[i])
                break
            else:
                b.write(f[i])

        p = i + 2
        wlist = []
        n = 0
        for i in range(p, len(f), 3):
            info_list = f[i+2].split(" ")
            print(info_list)
            pep1 = info_list[3]
            pep2 = info_list[5]
            site1 = info_list[1]
            site2 = info_list[2]
            if not_kk_links(pep1, pep2, site1, site2) and is_mono_link_pep(pep1, pep2, fasta_dic):
                n += 1
                wlist.append("[Spectrum%d]\n" % n)
                wlist.append(f[i+1])
                wlist.append(f[i+2])
        wlist.insert(0, "total=%d\n" % n)
        b.writelines(wlist)
        b.close()

