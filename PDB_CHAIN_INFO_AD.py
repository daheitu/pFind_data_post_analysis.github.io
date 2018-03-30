import os
os.chdir(r"C:\Users\Yong\Documents\pro-rna_pdb")
AA_dict = dict(
    HIS="H",
    MET="M",
    THR="T",
    PHE="F",
    PRO="P",
    SER="S",
    TRP="W",
    TYR="Y",
    VAL="V",
    M3L="m3K",
    GLY="G",
    ILE="I",
    ARG="R",
    LYS="K",
    LEU="L",
    ALA="A",
    CYS="C",
    ASN="N",
    GLN="Q",
    ASP="D",
    GLU='E',
    DPR="Dpr",
    MSE="Mse",
    MLZ="Mlz",
    TPO="Tpo",
    UNK="Uuk",
    SMC="Smc")
NA_list = ["A", "T", "C", "G", "U", "PSU", "5MU", "H2U"]
filename = os.listdir(os.getcwd())
out = open("output_1GST.txt", 'w')
for name in filename:
    if name[-4:] == ".pdb":
        print(name)
        pdb = open(name, 'r').readlines()
        chains = []
        for line in pdb:
            if line[0:6] == "COMPND":
                if line[11:17] == "CHAIN:":
                    print(line[18])
                    chains.append(line[18])
            elif line[0:6] == "SOURCE":
                break
        for chain in chains:
            chains[chains.index(chain)] = chain.lstrip()
        chains.sort()
        print(chains)

        chain_seq_pdb_dic = {}
        chain_length_pdb_dic = {}
        chain_type_pdb_dic = {}
        pdb_info = [name[:name.find(".")], str(len(chains))]
        for chain in chains:
            vars()["chain" + chain + "_seq_list"] = []
            chain_type = ""
            chain_length = ""
            for line in pdb:
                if line[:6] == "SEQRES" and line[11] == chain:
                    sequence = line[18:70].rstrip()
                    # print(sequence)
                    chain_length = line[13:17].strip()
                    raw_seq = []
                    if sequence[:5].strip() in AA_dict:
                        chain_type = "Protein"
                    else:
                        chain_type = "RNA"
                    element_num = int(len(sequence) / 4)
                    if chain_type == "Protein":
                        for i in range(element_num):
                            raw_seq.append(
                                AA_dict[sequence[i * 4:(i + 1) * 4].strip()])
                    else:
                        for i in range(element_num):
                            if sequence[i * 4:(i + 1) * 4].strip():
                                raw_seq.append(
                                    sequence[i * 4:(i + 1) * 4].strip())
                    vars()["chain" + chain + "_seq_list"].append(
                        "".join(raw_seq))
                elif line[:4] == "ATOM":
                    break
            chain_seq_pdb_dic[chain] = "".join(
                vars()["chain" + chain + "_seq_list"])
            chain_length_pdb_dic[chain] = chain_length
            chain_type_pdb_dic[chain] = chain_type
            chain_info = ",".join([
                chain,
                str(chain_length_pdb_dic[chain]), chain_type_pdb_dic[chain]
            ])
            pdb_info.append(chain_info)

        out.write("\t".join(pdb_info))
        out.write("\n")

out.close()