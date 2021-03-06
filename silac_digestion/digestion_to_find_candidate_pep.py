# coding = utf-8
import os
import gc


gc.enable()

fasta_path = r"./bsa.fasta"
max_misclavage = 3
max_pep_len = 30
min_pep_len = 6

sites_dic = {"N":"D", "C":"KR"}  # important, 


mpMassTable = {
        'A': 71.037114, 'R': 156.101111,
        'N':114.042927, 'D': 115.026943,
        'C': 103.009185, 'E': 129.042593,
        'Q': 128.058578, 'G': 57.021464,
        'H': 137.058912, 'I': 113.084064,
        'L': 113.084064, 'K': 128.094963,
        'M': 131.040485, 'F': 147.068414,
        'P': 97.052764, 'S': 87.032028,
        'T': 101.047679, 'U': 150.95363,
        'W': 186.079313, 'Y': 163.06332,
        'V': 99.068414, 'H2O': 18.01056,
        'H1': 1.00782
        }

mpModMass = {
        'Carbamidomethyl[C]': 57.021464,
        'Oxidation[M]': 15.994915
    }


def pretreatment_fasta(fasta_path):
    fasta = open(fasta_path).readlines()
    FastaDic = {}
    Pro_position_list = []
    i = 0 
    while i < len(fasta):
        if fasta[i][0] != ">":
            i += 1
        else:
            proLine = fasta[i]
            if " " in proLine:
                proName = proLine[1:proLine.find(" ")]
            else:
                proName = proLine[1:-1]
            seq = ""
            p = i + 1
            while p < len(fasta):
                if fasta[p][0] == ">":
                    break
                else:
                    seq += fasta[p].strip()
                p += 1

            FastaDic[proName] = seq
            i = p
    return FastaDic


def split_sequence(seq, sites_dic, start_pos, i):
    if seq[i] in sites_dic["N"] and seq[i] in sites_dic["C"]:
        left = seq[:i]
        mid = seq[i]
        right = seq[i+1:]
        return left, mid, right, start_pos + len(left) + 1
    elif seq[i] in sites_dic["N"]:
        left = seq[:i]
        mid = ""
        right = seq[i:]
        return left, mid, right, start_pos + len(left)
    else:
        left = seq[:i+1]
        mid = ""
        right = seq[i+1:]
        return left, mid, right, start_pos + len(left)
    

def judge_char(seq, i, sites_dic):
    sites_list = "".join(list(sites_dic.values()))
    if seq[i] in sites_dic["N"] and seq[i] in sites_dic["C"]:
        return True, i
    else:
        if i == 0 and seq[i] in sites_dic["N"]:
            return False
        elif i == len(seq)-1 and seq[i] in sites_dic["C"]:
            return False
        else:
            if seq[i] in sites_list:
                return True, i
            else:
                return False


def isContain_sites(seq, sites_dic):
    for i in range(len(seq)):
        if judge_char(seq, i, sites_dic):
            return True, i
    else:
        return False


def find_all_pep(seq, sites_dic):
    right = seq
    start_pos = 0
    rep_list = []
    start_pos_list = []
    while isContain_sites(right, sites_dic):
        start_pos_list.append(start_pos)    
        cleave_site = isContain_sites(right, sites_dic)[1]
        left, mid, right, start_pos = split_sequence(right, sites_dic, start_pos, cleave_site)
        print(left, mid, right, start_pos)
        if mid != "" and left != "":
            rep_list.extend([left, mid])
            start_pos_list.append(start_pos_list[-1]+len(left))
        else:
            rep_list.append("".join([left, mid]))
    rep_list.append(right)
    start_pos_list.append(start_pos)
    return rep_list, start_pos_list


def get_peps_idx(piecelist, miss_cleav, index_list):
    m = miss_cleav
    final_pep = []
    final_index_list = []
    for j in range(len(piecelist)-m):
        pep = "".join(piecelist[j:j+m+1])
        final_pep.append(pep)
        final_index_list.append(index_list[j])
    return final_pep, final_index_list


def getall_peps_idx(piecelist, max_mis, index_list):
    rep_list = [piecelist]
    rep_index_list = [index_list]
    for i in range(1, max_mis+1):
        cur_list, cur_index_list = get_peps_idx(piecelist, i, index_list)
        rep_list.append(cur_list)
        rep_index_list.append(cur_index_list)
    return rep_list, rep_index_list


def writePep2file_idx(rep_list, rep_index_list, pname, b):
    for i in range(len(rep_list)):
        for j in range(len(rep_list[i])):
            pep = rep_list[i][j]
            if min_pep_len <= len(pep) <= max_pep_len:
                start_index = rep_index_list[i][j]+1
                end_index = start_index+len(pep)
                pep_name = "_".join([pname, str(start_index), str(end_index)])
                b.write(">"+pep_name+"\n")
                b.write(pep+"\n")


def main_idx():
    fasta_dic = pretreatment_fasta(fasta_path)
    b = open(os.path.basename(fasta_path)+"_enzyme", 'w')
    print(">>>>>>>写入文件<<<<<<<<<<<<<<<<<")
    for pname in fasta_dic:
        seq = fasta_dic[pname]
        piece_peps, piece_idx = find_all_pep(seq, sites_dic)
        rep_list, rep_idx_list = getall_peps_idx(piece_peps, max_misclavage, piece_idx)
        writePep2file_idx(rep_list, rep_idx_list,pname,b)
    b.close()
    print(">>>>>>>写入完成<<<<<<<<<<<<<<<<<")




if __name__ == "__main__":
    main_idx()