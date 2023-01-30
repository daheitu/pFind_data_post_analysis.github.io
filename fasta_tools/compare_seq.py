from numpy import *


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


def get_name_seq_dic(fasta_file):
    name_seq_dic = {}
    f = open(fasta_file).readlines()
    
    i = 0
    while i < len(f):
        if not f[i].startswith(">"):
            i += 1
        else:
            pro_name = f[i].strip()
            seq = ""
            p = i + 1
            while p < len(f):
                if f[p].startswith(">"):
                    break
                else:
                    seq += f[p].strip()
                    p += 1
            name_seq_dic[pro_name] = seq
            i = p
    return name_seq_dic        


def get_overlap_pro(demo_dic, target_dic):
    over_dic = {}
    for pro in demo_dic.keys():
        over_dic[pro] = []
        for pro_tgt in target_dic:
            over = print_longest_common_subsequence(demo_dic[pro], target_dic[pro_tgt])
            if len(over)/ len(demo_dic[pro]) > 0.8 or len(over)/len(target_dic[pro_tgt]) > 0.8:
                # print(pro_tgt)
                over_dic[pro].append(pro_tgt)
    return over_dic            


def main():
    demo = r"D:\wechat_files\WeChat Files\yunzaitianya2010\FileStorage\File\2021-12\20211224_MS_Protein List(20-tag).txt"
    demo_dic = get_name_seq_dic(demo)
    targe_fasta = r"\\192.168.54.221\cy_data\pFind_work_space\yz_FSD_20210107\result\pBuild_tmp\translocon_FDS.fasta"
    target_dic = get_name_seq_dic(targe_fasta)
    
    over_dic = get_overlap_pro(demo_dic, target_dic)
    print(over_dic)
    overlap_proteins = []
    for pro in over_dic:
        overlap_proteins.extend(over_dic[pro])
    final_fasta = "./translocon_FDS_20220107.fasta"
    b = open(final_fasta, 'w')
    for pro in demo_dic:
        b.write(pro+"\n")
        b.write(demo_dic[pro] + "\n")

    for pro in target_dic:
        if pro not in overlap_proteins:
            b.write(pro+"\n")
            b.write(target_dic[pro] + "\n")
            

main()