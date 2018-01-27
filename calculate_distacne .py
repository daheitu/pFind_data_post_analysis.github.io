import math
import os
from numpy import *

#os.chdir(r"C:\Users\Yong\Documents\pLink\search_task_2018.01.18.09.19.12_PUD12_ARGO2_4RAW\reports")

AA_dict = dict(HIS="H", MET="M", THR="T", PHE="F", PRO="P", SER="S", TRP="W",
               TYR="Y", VAL="V", M3L="m3K", GLY="G", ILE="I", ARG="R", LYS="K",
               LEU="L", ALA="A", CYS="C", ASN="N", GLN="Q", ASP="D", GLU='E')

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

def cal_pdb_dis((chain_a, num1), (chain_b, num2)):
    #pdb=open(pdb_name,"r").readlines()
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
    return math.sqrt(sd)


def pretreatment_fasta(fastafile):
    fasta_dic = {}
    Pro_position_list =[]
    fasta = open(fastafile, 'r').readlines()
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
            for i in range(position+1, Pro_position_list[Pro_position_list.index(position)+1]):
                fasta_seq += fasta[i].strip()

            fasta_dic[Pro_name] = fasta_seq
    return fasta_dic




def pretreatment_pdb(pdbfile):
    pdb = open(pdbfile, 'r').readlines()
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
                    vars()["chain" + chain + "_seq_list"].append(AA_dict[sequence])
                else:
                    seq_length = len(sequence.split(" "))
                    for i in range(seq_length):
                        vars()["chain" + chain + "_seq_list"].append(
                            AA_dict[sequence.split(" ")[i]])
            elif line[:4] == "ATOM":
                break
        pdb_chain_To_seq_dic[chain] = "".join(vars()["chain" + chain + "_seq_list"])

    structuredCA_To_chain_dic = {}
    for chain in pdb_chains:
        positionCA_To_Num = {}
        for line in pdb:
            if  line[:4] == "ATOM" and line[21] == chain and line[13:15] == "CA":
                structured_CA_number = int(line[22:26])
                positionCA_To_Num[structured_CA_number] = AA_dict[line[17:20]]
        
        structuredCA_To_chain_dic[chain] = positionCA_To_Num

    structuredChain_To_seq ={}
    for chain in pdb_chains:
        structuredChain_To_seq[chain] = "".join(structuredCA_To_chain_dic[chain].values())

    return pdb, pdb_chains, pdb_chain_To_seq_dic, structuredCA_To_chain_dic, structuredChain_To_seq




def get_fasta_pdb_infor(path):
    filename = os.listdir(path)
    for file in filename:
        if file[-4:] == '.pdb':
            print "The pdb is" + file
            pdb, pdb_chains, pdb_chain_To_seq_dic, structuredCA_To_chain_dic, structuredChain_To_seq = pretreatment_pdb(file)
        elif file[-6:] == '.fasta':
            print "The fasta is" + file
            fasta_dic = pretreatment_fasta(file)
    
    PdbChain_To_ProName_dict = {}
    for chain in pdb_chain_To_seq_dic:
        for name in fasta_dic:
            overlap = "".join(
                print_longest_common_subsequence(fasta_dic[name], pdb_chain_To_seq_dic[chain]))
            if len(overlap) > len(pdb_chain_To_seq_dic[chain]) * 0.8 and len(
                    overlap) > len(fasta_dic[name]) * 0.4:
                PdbChain_To_ProName_dict[chain] = name

    Delta_PdbNum_To_FastaNum = {}
    for chain in pdb_chains:
        pdb_seq = structuredChain_To_seq[chain]
        fasta_seq = fasta_dic[PdbChain_To_ProName_dict[chain]]
        MaxSubstr_fastaTOpdb = print_longest_over_substrate(pdb_seq, fasta_seq)
        Substr_index_in_fasta = fasta_seq.index(MaxSubstr_fastaTOpdb)
        Substr_index_in_pdb = pdb_seq.index(MaxSubstr_fastaTOpdb)
        Substr_SeriesNum_in_Pdb = structuredCA_To_chain_dic[chain].keys()[Substr_index_in_pdb]
        Delta = Substr_index_in_fasta + 1 - Substr_SeriesNum_in_Pdb
        Delta_PdbNum_To_FastaNum[chain] = Delta
    
    return pdb, PdbChain_To_ProName_dict, Delta_PdbNum_To_FastaNum, structuredCA_To_chain_dic


pdb, PdbChain_To_ProName_dict, Delta_PdbNum_To_FastaNum, structuredCA_To_chain_dic = get_fasta_pdb_infor(os.getcwd())

print Delta_PdbNum_To_FastaNum, PdbChain_To_ProName_dict
Cross_link_site_list =["K", "R"]


def get_pdb_distance(cross_link_pair):
    protein1, protein2, position1, position2 = get_linked_site_inform(cross_link_pair)
    site_1_chain = []
    site_2_chain = []
    for chain in PdbChain_To_ProName_dict:
        if PdbChain_To_ProName_dict[chain] == protein1:
            Correct_Posi1 = position1 - Delta_PdbNum_To_FastaNum[chain]
            if Correct_Posi1 in structuredCA_To_chain_dic[chain] and structuredCA_To_chain_dic[chain][Correct_Posi1] in Cross_link_site_list:
                site_1_chain.append((chain, Correct_Posi1))
            else:
                return "W"

        if PdbChain_To_ProName_dict[chain] == protein2 :
            Correct_Posi2 = position2 - Delta_PdbNum_To_FastaNum[chain]
            if Correct_Posi2 in structuredCA_To_chain_dic[chain] and structuredCA_To_chain_dic[chain][Correct_Posi2] in Cross_link_site_list:
                site_2_chain.append((chain, Correct_Posi2))
            else:
                return "W"
    print site_1_chain, site_2_chain
    if len(site_1_chain) and len(site_2_chain):
        all_distance = []
        all_distance_dic = {}
        for i in range(len(site_1_chain)):
            for j in range(len(site_2_chain)):
                all_distance_dic[site_1_chain[i], site_2_chain[j]] = cal_pdb_dis(site_1_chain[i], site_2_chain[j])
        all_distance = all_distance_dic.values()
        if len(all_distance) == 1:
            return round(min(all_distance), 2)
        else:
            while min(all_distance) == 0.0:
                all_distance.remove(0.0)
            return round(min(all_distance), 2)
    else:
        return "no structure information"

def main():
    list = open("list.txt").readlines()
    B = open("report.txt", 'w')
    for pairs in list:
        distance = get_pdb_distance(pairs.strip())
        B.write("\t".join([pairs.strip(), str(distance)]))
        B.write("\n")
    B.close()
