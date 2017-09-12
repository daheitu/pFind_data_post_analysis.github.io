import math
import os
from numpy import *

#计算PDB里两个原子之间的距离，需给出chain name 和 atom number
def cal_pdb_dis((chain_a, num1), (chain_b, num2)):
    # pdb=open(pdb_name,"r").readlines()
    x_1 = 0
    x_2 = 0
    y_1 = 0
    y_2 = 0
    z_1 = 0
    z_2 = 0
    i = 1000
    for i in range(len(pdb)):
        if pdb[i][:4] == "ATOM" and str(pdb[i][13:15]) == 'CA':
            if str(pdb[i][21]) == chain_a and int(pdb[i][22:26]) == num1:
                x_1 = float(pdb[i][30:38])
                y_1 = float(pdb[i][38:46])
                z_1 = float(pdb[i][46:54])
    for i in range(len(pdb)):
        if pdb[i][:4] == "ATOM" and str(pdb[i][13:15]) == 'CA':
            if str(pdb[i][21]) == chain_b and int(pdb[i][22:26]) == num2:
                x_2 = float(pdb[i][30:38])
                y_2 = float(pdb[i][38:46])
                z_2 = float(pdb[i][46:54])
    sd = (x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1) + (z_2 - z_1) * (z_2 - z_1)
    return round(math.sqrt(sd), 2)

# 找到一个矩阵的最大值，以及所在的位置
def find_maxvaule_poistion(matrix):
    max = matrix[0][0]
    for i in range(len(matrix[:, 0])):
        for j in range(len(matrix[0, :])):
            if matrix[i][j] > max:
                max = matrix[i][j]
                poistion = (i, j)
    return poistion, max

#计算动态规划矩阵
def longest_common_subsequence(lhs, rhs):
    l_len = len(lhs)
    r_len = len(rhs)
    matrix = zeros((l_len, r_len))
    for i in arange(l_len):
        for j in arange(r_len):
            if lhs[i] == rhs[j]:
                if i != 0 and j != 0:
                    matrix[i][j] = matrix[i - 1][j - 1] + 1
                else:
                    matrix[i][j] = 1
            elif i != 0 and j != 0:
                matrix[i][j] = max(matrix[i - 1][j], matrix[i][j - 1])
    return matrix


def longest_over_substrate(lhs, rhs):
    l_len = len(lhs)
    r_len = len(rhs)
    matrix = zeros((l_len, r_len))
    for i in arange(l_len):
        for j in arange(r_len):
            if lhs[i] == rhs[j]:
                if i != 0 and j != 0:
                    matrix[i][j] = matrix[i - 1][j - 1] + 1
                else:
                    matrix[i][j] = 1
    return matrix


def print_longest_over_substrate(lhs, rhs):
    matrix = longest_over_substrate(lhs, rhs)
    substrate_length = find_maxvaule_poistion(matrix)[1]
    return lhs[int(find_maxvaule_poistion(matrix)[0][0] - substrate_length) + 1:int(
        find_maxvaule_poistion(matrix)[0][0]) + 1]

#找到
def print_longest_common_subsequence(lhs, rhs):
    matrix = longest_common_subsequence(lhs, rhs)
    l_len = len(lhs)
    r_len = len(rhs)
    i = l_len - 1
    j = r_len - 1
    rst = []
    while j > 0 and i > 0:
        if matrix[i][j] != matrix[i - 1][j] and matrix[i][j] != matrix[i][j - 1]:
            rst.insert(0, lhs[i])
            j -= 1
            i -= 1
        elif matrix[i][j] == matrix[i - 1][j]:
            i -= 1
        else:
            j -= 1
    if j != 0 and i == 0:
        if matrix[i][j - 1] != matrix[i][j]:
            rst.insert(0, lhs[i])
    elif j == 0 and 0:
        if matrix[i][j] != matrix[i - 1][j]:
            rst.insert(0, lhs[j])
    return rst


AA_dict = dict(HIS="H", MET="M", THR="T", PHE="F", PRO="P", SER="S", TRP="W", 
               TYR="Y", VAL="V", M3L="m3K", GLY="G", ILE="I", ARG="R", LYS="K",
               LEU="L", ALA="A", CYS="C", ASN="N", GLN="Q", ASP="D", GLU='E')

def link_fasta(num):
    s = ""
    while num < len(fasta) - 1 and fasta[num + 1][0] != ">":
        s += fasta[num + 1].strip()
        num += 1
    return s

def get_linked_site_inform(linked_site):
    m = linked_site.find("(")
    n = linked_site.find(")")
    p = linked_site.find("-")
    x = linked_site.find("(", p)
    y = linked_site.find(")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:x].strip()
    position1 = int(linked_site[m + 1:n])
    position2 = int(linked_site[x + 1:y])
    return protein1, protein2, position1, position2

def get_pdb_distance(cross_link_pair):
    (protein1, protein2, position1, position2) = get_linked_site_inform(cross_link_pair)

    def judge_cross_link_type(protein1, protein2):
        if protein1 == protein2:
            return "intra"
        else:
            return "inter"

    if not chain_name_dict:
        return "PDB is not consistent with FASTA ", judge_cross_link_type(protein1, protein2)
    else:
        if protein1 not in chain_name_dict.values() or protein2 not in chain_name_dict.values():
            return "no structure infor", judge_cross_link_type(protein1, protein2)
        else:
            protein1_chain = []
            protein2_chain = []
            for chain in chain_name_dict.keys():
                if chain_name_dict[chain] == protein1 and (position1 - fasta_pdb_seq_delta_dict[chain]) in dic_list_dic[
                    chain]:
                    if dic_list_dic[chain][position1 - fasta_pdb_seq_delta_dict[chain]] == site_1 or \
                                    dic_list_dic[chain][position1 - fasta_pdb_seq_delta_dict[chain]] == site_2 or \
                                    position1 == 1:
                        protein1_chain.append((chain, position1 - fasta_pdb_seq_delta_dict[chain]))
                    else:
                        return "wrong", "wrong"
            for chain in chain_name_dict.keys():
                if chain_name_dict[chain] == protein2 and (position2 - fasta_pdb_seq_delta_dict[chain]) in dic_list_dic[
                    chain]:
                    if dic_list_dic[chain][position2 - fasta_pdb_seq_delta_dict[chain]] == site_1 or \
                                    dic_list_dic[chain][position2 - fasta_pdb_seq_delta_dict[chain]] == site_2 or \
                                    position1 == 1:
                        protein2_chain.append((chain, position2 - fasta_pdb_seq_delta_dict[chain]))
                    else:
                        return "wrong", "wrong"
            print protein1_chain, protein2_chain
            if len(protein1_chain) and len(protein2_chain):
                all_distance = []
                all_distance_dic = {}
                for i in range(len(protein1_chain)):
                    for j in range(len(protein2_chain)):
                        # all_distance.append(cal_pdb_dis(protein1_chain[i],protein2_chain[j]))
                        all_distance_dic[protein1_chain[i], protein2_chain[j]] = cal_pdb_dis(protein1_chain[i],
                                                                                             protein2_chain[j])
                print all_distance_dic
                all_distance = all_distance_dic.values()
                if len(all_distance) == 1:
                    if int(min(all_distance)) == 0:
                        return "self cross link", "inter"
                    else:
                        return round(min(all_distance), 2), judge_cross_link_type(protein1, protein2)
                else:
                    while min(all_distance) == 0.0:
                        all_distance.remove(0.0)
                    for word in all_distance_dic:
                        if all_distance_dic[word] == min(all_distance):
                            chain_1 = word[0][0]
                            chain_2 = word[1][0]
                    if chain_1 != chain_2:
                        return round(min(all_distance), 2), "inter"
                    else:
                        return round(min(all_distance), 2), "intra"
            else:
                return "no structure information", judge_cross_link_type(protein1, protein2)

os.chdir(r"C:\Users\Yong\Desktop\Cal_Pdb_Distance")
filename = os.listdir(os.getcwd())
for file in filename:
    if file[-4:] == '.pdb':
        pdb = open(file, 'r').readlines()
    elif file[-6:] == '.fasta':
        fasta = open(file, 'r').readlines()
    elif file[-22:-4] == "cross-linked_sites":
        site_table = open(file, 'r').readlines()



if not fasta:
    print "no fasta file"
else:
    if not pdb:
        print "no pdb"
    else:
        name_seq_in_fasta_dic = {}
        for line in fasta:
            if line[0] == ">":
                if " " in line:
                    name = line[1:line.strip().find(" ")]
                elif "|" in line:
                    name = line[1:line.strip().find("|")]
                else:
                    name = line[1:len(line.strip())]
                name_seq_in_fasta_dic[name] = link_fasta(fasta.index(line))
        print name_seq_in_fasta_dic

chains = []
for line in pdb:
    if line[0:6] == "COMPND":
        if line[11:17] == "CHAIN:":
            chain_length = len(line[17:len(line.strip()) - 1])
            if chain_length == 2:
                chains.append(line[17:len(line.strip()) - 1])
            else:
                chain_num = len(line[17:len(line.strip()) - 1].split(","))
                for i in range(chain_num):
                    chains.append(line[17:len(line.strip()) - 1].split(",")[i])
    elif line[0:6] == "SOURCE":
        break
for chain in chains:
    chains[chains.index(chain)] = chain.lstrip()
print chains


chain_name_dict = {}
chain_seq_pdb_dic = {}
for chain in chains:
    vars()["chain" + chain + "_seq_list"] = []
    for line in pdb:
        if line[:6] == "SEQRES" and line[11] == chain:
            sequence = line[19:len(line) - 1].strip()
            if len(sequence) == 3:
                vars()["chain" + chain + "_seq_list"].append(AA_dict[sequence])
            else:
                seq_length = len(sequence.split(" "))
                for i in range(seq_length):
                    vars()["chain" + chain + "_seq_list"].append(AA_dict[sequence.split(" ")[i]])
        elif line[:4] == "ATOM":
            break
    chain_seq_pdb_dic[chain] = "".join(vars()["chain" + chain + "_seq_list"])

    for name in name_seq_in_fasta_dic:
        overlap = "".join(print_longest_common_subsequence(name_seq_in_fasta_dic[name], chain_seq_pdb_dic[chain]))
        if len(overlap) > len(chain_seq_pdb_dic[chain]) * 0.8 and len(overlap) > len(name_seq_in_fasta_dic[name]) * 0.4:
            chain_name_dict[chain] = name
print chain_name_dict

fasta_pdb_seq_delta_dict = {}
dic_list_dic = {}
seq2_list_dic = {}
for chain in chains:
    vars()["chain" + chain + "_list"] = {}
    for line in pdb:
        if line[:4] == "ATOM" and line[21] == chain:
            vars()["chain" + chain + "_list"][int(line[22:26].lstrip())] = AA_dict[line[17:20]]
    dic_list_dic[chain] = vars()["chain" + chain + "_list"]
    seq2_list_dic[chain] = "".join(vars()["chain" + chain + "_list"].values())
    max_substr_of_fasta_pdb = print_longest_over_substrate(seq2_list_dic[chain],
                                                           name_seq_in_fasta_dic[chain_name_dict[chain]])
    sub_str_index_in_fasta = name_seq_in_fasta_dic[chain_name_dict[chain]].find(max_substr_of_fasta_pdb)
    sub_str_index_in_seq = chain_seq_pdb_dic[chain].find(max_substr_of_fasta_pdb)
    sub_str_index_in_seq2 = seq2_list_dic[chain].find(max_substr_of_fasta_pdb)
    series_num_in_pdb = dic_list_dic[chain].keys()[sub_str_index_in_seq2]
    fasta_pdb_seq_delta_dict[chain] = sub_str_index_in_fasta + 1 - series_num_in_pdb
print fasta_pdb_seq_delta_dict

site_1="K"
site_2="K"




list=open("list.txt").readlines()
B=open("report.txt",'w')
for pairs in list:
    Prot1 = get_linked_site_inform(pairs)[0]
    Prot2 = get_linked_site_inform(pairs)[1]
    if Prot1 == "Nhp2" or Prot2 == "Nhp2":
        distance="no"
    else:
        distance = get_pdb_distance(pairs.strip())[0]
    B.write("\t".join([pairs.strip(),str(distance)]))
    B.write("\n")
B.close()


