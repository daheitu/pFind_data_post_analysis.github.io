import math
import os
import re
from numpy import *


path = r"F:\Script\pdbDatabase\3u28_CNG"
XL_sites_list = ["R", "K"]
os.chdir(path)
input_file = "BDG.txt"
output_file = input_file + "_out.txt"


AA_dict = dict(
    HIS="H", MET="M", THR="T",
    PHE="F", PRO="P", SER="S",
    TRP="W", TYR="Y", VAL="V",
    M3L="m3K", GLY="G", ILE="I",
    ARG="R", LYS="K", LEU="L",
    ALA="A", CYS="C", ASN="N",
    GLN="Q", ASP="D", GLU='E')


def getLinkedSiteInfo(linked_site):
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


def findCAcoord(chain_a, num1, pdb):
    for line in pdb:
        if line[:4] == "ATOM" and str(line[13:15]) == 'CA':
            if str(line[21]) == chain_a and int(line[22:26]) == num1:
                x_1 = float(line[30:38])
                y_1 = float(line[38:46])
                z_1 = float(line[46:54])
    return x_1, y_1, z_1


def cal_pdb_dis(chain_a, num1, chain_b, num2, pdb):
    xA, yA, zA = findCAcoord(chain_a, num1, pdb)
    xB, yB, zB = findCAcoord(chain_b, num2, pdb)
    sd = (xA - xB) ** 2 + (yA - yB) ** 2 + (zA - zB) ** 2
    return round(math.sqrt(sd), 2)


def pretreatment_fasta(fasta):
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


def pretreatment_pdb(pdb):
    pdb_chains = []
    for line in pdb:
        if line[0:6] == "ATOM  ":
            chain = line[21].strip()
            if chain not in pdb_chains:
                pdb_chains.append(chain)
        else:
            continue

    pdbChainToSeqDic = {}
    for chain in pdb_chains:
        pdbChainToSeqDic[chain] = []
        for line in pdb:
            if line[:6] == "SEQRES" and line[11] == chain:
                sequence = line[19:len(line) - 1].strip()
                if len(sequence) == 3:
                    pdbChainToSeqDic[chain].append(
                        AA_dict[sequence])
                else:
                    aaList = sequence.split(" ")
                    for aa3 in aaList:
                        pdbChainToSeqDic[chain].append(AA_dict[aa3])
            elif line[:4] == "ATOM":
                break
        pdbChainToSeqDic[chain] = "".join(pdbChainToSeqDic[chain])

    strucCAToChainDic = {}
    for chain in pdb_chains:
        positionCA_To_Num = {}
        for line in pdb:
            if line[:4] == "ATOM":
                if line[21] == chain and line[13:15] == "CA":
                    struc_CA_num = int(line[22:26])
                    positionCA_To_Num[struc_CA_num] = AA_dict[line[17:20]]

        strucCAToChainDic[chain] = positionCA_To_Num

    structChainToSeq = {}
    for chain in pdb_chains:
        structChainToSeq[chain] = "".join(list(strucCAToChainDic[chain].values()))

    return pdb_chains, pdbChainToSeqDic, strucCAToChainDic, structChainToSeq


def correlateFastaPdb(fasta, pdb):
    FastaDic = pretreatment_fasta(fasta)
    pdb_chains = pretreatment_pdb(pdb)[0]
    chainToSeqInpdbDic = pretreatment_pdb(pdb)[1]
    structCAtoChainDic = pretreatment_pdb(pdb)[2]
    structChainToSeq = pretreatment_pdb(pdb)[3]
    print(FastaDic, pdb_chains)
    structCAtoChainDic
    chainInpdbToPronameDic = {}
    for chain in chainToSeqInpdbDic:
        for name in FastaDic:
            overlap = "".join(
                print_longest_common_subsequence(FastaDic[name],
                                                 chainToSeqInpdbDic[chain]))
            if len(overlap) > len(chainToSeqInpdbDic[chain]) * 0.8 and len(
                    overlap) > len(FastaDic[name]) * 0.4:
                chainInpdbToPronameDic[chain] = name

    deltaPdbNumToFastaNum = {}
    for chain in pdb_chains:
        pdb_seq = structChainToSeq[chain]  # 有结构的PDB 的序列
        fasta_seq = FastaDic[chainInpdbToPronameDic[chain]]  # fasta 序列
        MaxSubstr_fastaTOpdb = print_longest_over_substrate(pdb_seq, fasta_seq)
        Substr_index_in_fasta = fasta_seq.index(MaxSubstr_fastaTOpdb)
        Substr_index_in_pdb = pdb_seq.index(MaxSubstr_fastaTOpdb)
        Substr_SeriesNum_in_Pdb = list(
            structCAtoChainDic[chain].keys())[Substr_index_in_pdb]
        Delta = Substr_index_in_fasta + 1 - Substr_SeriesNum_in_Pdb
        deltaPdbNumToFastaNum[chain] = Delta

    return chainInpdbToPronameDic, deltaPdbNumToFastaNum


def calPdbDistance(linkXpair, chainInpdbToPronameDic,
                     deltaPdbNumToFastaNum, structCAtoChainDic, pdb):
    ptn1, ptn2, pos1, pos2 = getLinkedSiteInfo(linkXpair)
    chaProN = chainInpdbToPronameDic
    deltaPdbFasta = deltaPdbNumToFastaNum
    strCaCha = structCAtoChainDic
    site_1_chain = calSiteToPdbPos(ptn1, pos1, chaProN, deltaPdbFasta, strCaCha)
    site_2_chain = calSiteToPdbPos(ptn2, pos2,  chaProN, deltaPdbFasta, strCaCha)
    print(site_1_chain, site_2_chain)
    if len(site_1_chain) and len(site_2_chain):
        all_distance = []
        all_dist_dic = {}
        for i in range(len(site_1_chain)):
            chain1, site1 = site_1_chain[i]
            for j in range(len(site_2_chain)):            
                chain2, site2 = site_2_chain[j]
                dis_pdb = cal_pdb_dis(chain1, site1, chain2, site2, pdb)
                all_dist_dic[site_1_chain[i], site_2_chain[j]] = [str(dis_pdb)]
                all_distance.append(dis_pdb)
        
        all_distance = delEleFromList(all_distance, 0.00)
        if all_distance == []:
            min_dis = "Self_XL"
        else:
            min_dis = str(min(all_distance))
        return min_dis
    else:
        return "No_Stru"


def calSiteToPdbPos(ptn1, pos1, chainInpdbToPronameDic, deltaPdbNumToFastaNum, structCAtoChainDic):
    site_1_chain = []
    chain_list = list(chainInpdbToPronameDic.keys())
    chain_list.sort()
    for chain in chain_list:
        if chainInpdbToPronameDic[chain] == ptn1:
            pdb_Posi1 = int(pos1) - deltaPdbNumToFastaNum[chain]
            if pdb_Posi1 in structCAtoChainDic[chain]:
                if structCAtoChainDic[chain][pdb_Posi1] in XL_sites_list:
                    site_1_chain.append((chain, pdb_Posi1))

    return site_1_chain


def delEleFromList(eList, ele):
    if ele not in eList:
        return eList
    else:
        n = 0
        while n < len(eList):
            if eList[n] == ele:
                del eList[n]
            else:
                n += 1
        return eList


def geneXLjwalkCMD():
    f = open(output_file, 'r').readlines()
    c = open(output_file + "_jwalk", 'w')
    for line in f:
        lineList = line.split("\t")
        jwalkCMD = lineList[-1]
        if jwalkCMD.strip() != "No_Stru":
            c.write(jwalkCMD)
    c.close()


def main():
    f = open(input_file, 'r').readlines()
    pair_list = []
    for line in f:
        line_list = line.strip().split("\t")
        if line_list[0].isdigit():
            break
        else:
            pair_list.append(line_list[0])
    
    B = open(output_file, 'w')
    flList = os.listdir(os.getcwd())
    for fl in flList:
        if fl[-6:] == ".fasta":
            fasta = open(fl, 'r').readlines()
            print("The fasta file is " + fl)
        elif fl[-4:] == ".pdb":
            pdb_name = fl[:-4]
            pdb = open(fl, 'r').readlines()
            print("The pdb file is " + pdb_name)
        else:
            continue

    chainInpdbToPronameDic = correlateFastaPdb(fasta, pdb)[0]
    deltaPdbNumToFastaNum = correlateFastaPdb(fasta, pdb)[1]
    structCAtoChainDic = pretreatment_pdb(pdb)[2]
    print(chainInpdbToPronameDic)
    print(deltaPdbNumToFastaNum)
    for pairs in pair_list:
        ptn1, ptn2, pos1, pos2 = getLinkedSiteInfo(pairs)
        site_1_chain = calSiteToPdbPos(ptn1, pos1, chainInpdbToPronameDic, deltaPdbNumToFastaNum, structCAtoChainDic)
        site_2_chain = calSiteToPdbPos(ptn2, pos2, chainInpdbToPronameDic, deltaPdbNumToFastaNum, structCAtoChainDic)
        print(site_1_chain, site_2_chain)
        
        if len(site_1_chain) and len(site_2_chain):
            all_distance = []
            all_dist_dic = {}
            for i in range(len(site_1_chain)):
                chain1, site1 = site_1_chain[i]
                for j in range(len(site_2_chain)):            
                    chain2, site2 = site_2_chain[j]
                    dis_pdb = cal_pdb_dis(chain1, site1, chain2, site2, pdb)
                    obj_name = "dist " + chain1 + pos1 + "_" + chain2 + pos2
                    selec1 = " /" + pdb_name + "//" + chain1 + "/" + str(site1) + "/"+"CA"
                    selec2 = " /" + pdb_name + "//" + chain2 + "/" + str(site2) + "/"+"CA"
                    pym_cmd = obj_name + "," + selec1 + "," + selec2
                    jwalkPos = "|".join([str(site1), chain1,  str(site2), chain2, ""])
                    all_dist_dic[site_1_chain[i], site_2_chain[j]] = [str(dis_pdb), jwalkPos]
                    all_distance.append(dis_pdb)

            all_distance = delEleFromList(all_distance, 0.00)
            if all_distance == []:
                min_dis = "Self_XL"
            else:
                min_dis = str(min(all_distance))
            line_write = [pairs, min_dis]
            if all_dist_dic:
                for key in all_dist_dic:
                    line_write.extend(all_dist_dic[key])
        else:
            line_write = [pairs, "No_Stru"]    

        B.write("\t".join(line_write))
        B.write("\n")

    B.close()
    


if __name__ == "__main__":
    main()
    geneXLjwalkCMD()