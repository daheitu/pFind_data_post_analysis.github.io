import math
import os
import re
from numpy import *


path = r"G:\DSSO_1008\structure"
XL_sites_list = ["R", "K"]
os.chdir(path)

AA_dict = dict(
    HIS="H", MET="M", THR="T",
    PHE="F", PRO="P", SER="S",
    TRP="W", TYR="Y", VAL="V",
    M3L="m3K", GLY="G", ILE="I",
    ARG="R", LYS="K", LEU="L",
    ALA="A", CYS="C", ASN="N",
    GLN="Q", ASP="D", GLU='E')


def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    return protein1, protein2, position1, position2


def find_maxvaule_poistion(matrix):
    max = matrix[0][0]
    for i in range(len(matrix[:, 0])):
        for j in range(len(matrix[0, :])):
            if matrix[i][j] > max:
                max = matrix[i][j]
                poistion = (i, j)
    return poistion, max


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
    return lhs[int(find_maxvaule_poistion(matrix)[0][0] - substrate_length) +
               1:int(find_maxvaule_poistion(matrix)[0][0]) + 1]


def print_longest_common_subsequence(lhs, rhs):
    matrix = longest_common_subsequence(lhs, rhs)
    l_len = len(lhs)
    r_len = len(rhs)
    i = l_len - 1
    j = r_len - 1
    rst = []
    while j > 0 and i > 0:
        if matrix[i][j] != matrix[i - 1][j] and matrix[i][j] != matrix[i][j
                                                                          - 1]:
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


def cal_pdb_dis(chain_a, num1, chain_b, num2, pdb):
    x_1 = 0
    x_2 = 0
    y_1 = 0
    y_2 = 0
    z_1 = 0
    z_2 = 0
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
    sd = (x_2 - x_1) * (x_2 - x_1) + (y_2 - y_1) * (y_2 - y_1) + (
        z_2 - z_1) * (z_2 - z_1)
    return round(math.sqrt(sd), 2)


def pretreatment_fasta(fasta):
    FastaDic = {}
    Pro_position_list = []
    for line in fasta:
        if line[0] == ">":
            Pro_position_list.append(fasta.index(line))
    Pro_position_list.append(len(fasta))

    for position in Pro_position_list:
        if position == len(fasta):
            pass
        else:
            Pro_line = fasta[position].rstrip("\n")
            if " " in Pro_line:
                Pro_name = Pro_line[1:Pro_line.find(" ")]
            else:
                Pro_name = Pro_line[1:]

            fasta_seq = ""
            for i in range(
                    position + 1,
                    Pro_position_list[Pro_position_list.index(position) + 1]):
                fasta_seq += fasta[i].strip()

            FastaDic[Pro_name] = fasta_seq
    return FastaDic


def pretreatment_pdb(pdb):
    pdb_chains = []
    for line in pdb:
        if line[0:6] == "SEQRES":
            pdb_chains.append(line[11:13])
        pdb_chains = list(set(pdb_chains))
        for chain in pdb_chains:
            pdb_chains[pdb_chains.index(chain)] = chain.rstrip()

    pdb_chain_To_seq_dic = {}
    for chain in pdb_chains:
        vars()["chain" + chain + "_seq_list"] = []
        for line in pdb:
            if line[:6] == "SEQRES" and line[11] == chain:
                sequence = line[19:len(line) - 1].strip()
                if len(sequence) == 3:
                    vars()["chain" + chain + "_seq_list"].append(
                        AA_dict[sequence])
                else:
                    seq_length = len(sequence.split(" "))
                    for i in range(seq_length):
                        vars()["chain" + chain + "_seq_list"].append(
                            AA_dict[sequence.split(" ")[i]])
            elif line[:4] == "ATOM":
                break
        pdb_chain_To_seq_dic[chain] = "".join(
            vars()["chain" + chain + "_seq_list"])

    structuredCA_To_chain_dic = {}
    for chain in pdb_chains:
        positionCA_To_Num = {}
        for line in pdb:
            if line[:4] == "ATOM" and line[21] == chain and line[13:
                                                                 15] == "CA":
                structured_CA_number = int(line[22:26])
                positionCA_To_Num[structured_CA_number] = AA_dict[line[17:20]]

        structuredCA_To_chain_dic[chain] = positionCA_To_Num

    structuredChain_To_seq = {}
    for chain in pdb_chains:
        structuredChain_To_seq[chain] = "".join(
            structuredCA_To_chain_dic[chain].values())

    return pdb_chains, pdb_chain_To_seq_dic, structuredCA_To_chain_dic, structuredChain_To_seq


def get_fasta_pdb_infor(fasta, pdb):
    FastaDic = pretreatment_fasta(fasta)
    pdb_chains = pretreatment_pdb(pdb)[0]
    pdb_chain_To_seq_dic = pretreatment_pdb(pdb)[1]
    structuredCA_To_chain_dic = pretreatment_pdb(pdb)[2]
    structuredChain_To_seq = pretreatment_pdb(pdb)[3]

    PdbChain_To_ProName_dict = {}
    for chain in pdb_chain_To_seq_dic:
        for name in FastaDic:
            overlap = "".join(
                print_longest_common_subsequence(FastaDic[name],
                                                 pdb_chain_To_seq_dic[chain]))
            if len(overlap) > len(pdb_chain_To_seq_dic[chain]) * 0.8 and len(
                    overlap) > len(FastaDic[name]) * 0.4:
                PdbChain_To_ProName_dict[chain] = name

    Delta_PdbNum_To_FastaNum = {}
    for chain in pdb_chains:
        pdb_seq = structuredChain_To_seq[chain]  # 有结构的PDB 的序列
        fasta_seq = FastaDic[PdbChain_To_ProName_dict[chain]]  # fasta 序列
        MaxSubstr_fastaTOpdb = print_longest_over_substrate(pdb_seq, fasta_seq)
        Substr_index_in_fasta = fasta_seq.index(MaxSubstr_fastaTOpdb)
        Substr_index_in_pdb = pdb_seq.index(MaxSubstr_fastaTOpdb)
        Substr_SeriesNum_in_Pdb = list(
            structuredCA_To_chain_dic[chain].keys())[Substr_index_in_pdb]
        Delta = Substr_index_in_fasta + 1 - Substr_SeriesNum_in_Pdb
        Delta_PdbNum_To_FastaNum[chain] = Delta

    return PdbChain_To_ProName_dict, Delta_PdbNum_To_FastaNum


def get_pdb_distance(cross_link_pair, PdbChain_To_ProName_dict,
                     Delta_PdbNum_To_FastaNum, structuredCA_To_chain_dic, pdb):
    protein1, protein2, position1, position2 = get_linked_site_inform(
        cross_link_pair)
    site_1_chain = []
    site_2_chain = []
    for chain in PdbChain_To_ProName_dict:
        if PdbChain_To_ProName_dict[chain] == protein1:
            Correct_Posi1 = int(position1) - Delta_PdbNum_To_FastaNum[chain]
            if Correct_Posi1 in structuredCA_To_chain_dic[chain]:
                if structuredCA_To_chain_dic[chain][Correct_Posi1] in XL_sites_list:
                    site_1_chain.append((chain, Correct_Posi1))
                else:
                    return "not cross-linked AA"
            else:
                return "No struc info"

        if PdbChain_To_ProName_dict[chain] == protein2:
            Correct_Posi2 = int(position2) - Delta_PdbNum_To_FastaNum[chain]
            if Correct_Posi2 in structuredCA_To_chain_dic[chain]:
                if structuredCA_To_chain_dic[chain][Correct_Posi2] in XL_sites_list:
                    site_2_chain.append((chain, Correct_Posi2))
                else:
                    return "not cross-linked AA"
            else:
                return "No struc info"
    print(site_1_chain, site_2_chain)
    if len(site_1_chain) and len(site_2_chain):
        all_distance = []
        all_distance_dic = {}
        for i in range(len(site_1_chain)):
            for j in range(len(site_2_chain)):
                all_distance_dic[
                    site_1_chain[i], site_2_chain[j]] = cal_pdb_dis(
                        site_1_chain[i][0], site_1_chain[i][1],
                        site_2_chain[j][0], site_2_chain[j][1], pdb)
        all_distance = list(all_distance_dic.values())
        if len(all_distance) == 1:
            return round(min(all_distance), 2)
        else:
            while min(all_distance) == 0.0:
                all_distance.remove(0.0)
            return round(min(all_distance), 2)
    else:
        return "no structure information"


def main():
    f = open("DSS.txt").readlines()
    pair_list = []
    for line in f[1:]:
        line_list = line.strip().split("\t")
        if line_list[0].isdigit():
            break
        else:
            pair_list.append(line_list[0])
    
    B = open("DSS_DISTANCE.txt", 'w')
    file_list = os.listdir(os.getcwd())
    for fl in file_list:
        if fl[-6:] == ".fasta":
            fasta = open(fl, 'r').readlines()
            print("The fasta file is " + fl)
        elif fl[-4:] == ".pdb":
            pdb_name = fl[:-4]
            pdb = open(fl, 'r').readlines()
            print("The pdb file is " + pdb_name)
        else:
            continue

    PdbChain_To_ProName_dict = get_fasta_pdb_infor(fasta, pdb)[0]
    Delta_PdbNum_To_FastaNum = get_fasta_pdb_infor(fasta, pdb)[1]
    structuredCA_To_chain_dic = pretreatment_pdb(pdb)[2]
    print(PdbChain_To_ProName_dict)
    print(Delta_PdbNum_To_FastaNum)
    for pairs in pair_list:
        ptn1, ptn2, pos1, pos2 = get_linked_site_inform(pairs)
        site_1_chain = []
        site_2_chain = []
        chain_list = list(PdbChain_To_ProName_dict.keys())
        chain_list.sort()
        for chain in chain_list:
            if PdbChain_To_ProName_dict[chain] == ptn1:
                pdb_Posi1 = int(pos1) - Delta_PdbNum_To_FastaNum[chain]
                if pdb_Posi1 in structuredCA_To_chain_dic[chain]:
                    if structuredCA_To_chain_dic[chain][pdb_Posi1] in XL_sites_list:
                        site_1_chain.append((chain, pdb_Posi1))
                    else:
                        pass
                else:
                    pass
            else:
                continue

        for chain in chain_list:
            if PdbChain_To_ProName_dict[chain] == ptn2:
                pdb_Posi2 = int(pos2) - Delta_PdbNum_To_FastaNum[chain]
                if pdb_Posi2 in structuredCA_To_chain_dic[chain]:
                    if structuredCA_To_chain_dic[chain][pdb_Posi2] in XL_sites_list:
                        site_2_chain.append((chain, pdb_Posi2))
                    else:
                        pass
                else:
                    pass
            else:
                continue
        print(site_1_chain, site_2_chain)
        
        if len(site_1_chain) and len(site_2_chain):
            all_distance = []
            all_dist_dic = {}
            jwalkCMDlist = []
            for i in range(len(site_1_chain)):
                chain1, site1 = site_1_chain[i]
                ori_sit1 = int(site1) + Delta_PdbNum_To_FastaNum[chain1]
                for j in range(len(site_2_chain)):            
                    chain2, site2 = site_2_chain[j]
                    ori_site2 = int(site2) + Delta_PdbNum_To_FastaNum[chain2]
                    dis_pdb_ED = cal_pdb_dis(chain1, site1, chain2, site2, pdb)
                    all_distance.append(dis_pdb_ED)
                    jwalkCMD = "|".join([site1,chain1, site2, chain2, ""])
                    jwalkCMDlist.append(jwalkCMD)
                    all_dist_dic[site_1_chain[i], site_2_chain[j]] = [str(dis_pdb), jwalk_position]
            n = 0
            while n < len(all_distance):
                if all_distance[n] == 0.00:
                    all_distance.remove(0.00)
                else:
                    n += 1
            if all_distance == []:
                min_dis = "self-cross-link"
            else:
                min_dis = str(min(all_distance))
            line_write = [pairs, min_dis]
            if all_dist_dic:
                for key in all_dist_dic:
                    line_write.extend(all_dist_dic[key])
        else:
            line_write = [pairs, "no Stru"]    

        B.write("\t".join(line_write))
        B.write("\n")

    B.close()


if __name__ == "__main__":
    main()
