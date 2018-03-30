import os
os.chdir(r"D:\fasta database")


def statstic_aa_abundance(fasta_file):
    aa_num_dic = {}
    aa_freq_dic = {}
    f = open(fasta_file).readlines()
    for line in f:
        if line:
            if line[0] == ">":
                continue
            else:
                for char in line.strip():
                    if char not in aa_num_dic:
                        aa_num_dic[char] = 1
                    else:
                        aa_num_dic[char] += 1
        else:
            continue
    total_aa = sum(aa_num_dic.values())
    for aa in aa_num_dic:
        aa_freq_dic[aa] = aa_num_dic[aa] / total_aa
    return sorted(aa_freq_dic.items(), key=lambda item: item[1]) # 按照value排序


print(statstic_aa_abundance("uniprot-proteome_UP000000589_Mus musculus.fasta"))
