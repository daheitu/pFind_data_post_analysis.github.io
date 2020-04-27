# /usr/bin/python
# coding = utf-8
import os, random


def extractFasta(fastaFL):
    f = open(fastaFL, 'r').readlines()
    num_fasta_dic = {}
    i = 0
    m = 0
    while i < len(f):
        if f[i][0] != ">":
            i =+ 1
        else:
            m += 1
            lns = [f[i]]
            p = i + 1
            while p < len(f):
                if f[p][0] == ">":
                    break
                else:
                    lns.append(f[p])
                    p += 1
            num_fasta_dic[m] = lns
        i = p
    return num_fasta_dic


def del_ele_seq(num_seq, sub_seq):
    for num in sub_seq:
        num_seq.remove(num)


def random_num_seq(num_fasta_dic, k_pro):
    num_pro = len(num_fasta_dic)
    num_seq = list(range(1, num_pro + 1))
    a1 = sorted(random.sample(num_seq, k_pro))
    print(sorted(a1))
    del_ele_seq(num_seq, a1)
    a2 = sorted(random.sample(num_seq, k_pro))
    print(sorted(a2))
    del_ele_seq(num_seq, a2)
    a3 = sorted(random.sample(num_seq, k_pro))
    print(sorted(a3))
    return a1, a2, a3


def add_BSA2fasta(comFasta):
    bsa = open(comFasta, 'r').readlines()
    for fl in os.listdir():
        if "E.coli_random600" in fl:
            b = open(fl, 'a')
            for line in bsa:
                b.write(line)
            b.close()


def combinefasta(fastaFL_list, comb_name):
    b = open(comb_name, 'w')
    for fl in fastaFL_list:
        f = open(fl, 'r').readlines()
        for line in f:
            b.write(line)
    b.close()


def main():
    num_fasta_dic = extractFasta("./uniprot_Ecoli_K12.fasta")
    a1, a2, a3 = random_num_seq(num_fasta_dic, 600)
    for ls in [a1, a2, a3]:
        b = open("E.coli_random600_"+str([a1, a2, a3].index(ls)+ 1)+".fasta", 'w')
        for num in ls:
            wlist = num_fasta_dic[num]
            for line in wlist:
                b.write(line)
        b.close()

if __name__ == "__main__":
    #main()
    #add_BSA2fasta('./synthetic_pep.fasta')
    #add_BSA2fasta("./contaminant.fasta")
    combinefasta(['./synthetic_pep.fasta', './contaminant.fasta', './uniprot_Ecoli_K12.fasta'], "synth_pep_cont_E.coli.fasta")
    combinefasta(['./synthetic_pep.fasta', './contaminant.fasta'], "synth_pep_contaminant.fasta")