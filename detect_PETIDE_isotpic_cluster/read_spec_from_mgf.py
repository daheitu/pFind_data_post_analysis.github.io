def read_spec(mgf_file):
    mgf_info_dic = {}
    f = open(mgf_file, 'r').readlines()
    i = 0
    while i < len(f):
        # print("current is " + str(i))
        if f[i].strip() == "BEGIN IONS":
            # begin_idx = i
            # total_spec += 1
            pass
        else:
            print("wrong")
        title = f[i+1].strip().split("=")[1]
        chrg = int([i for i in f[i+2] if i.isdigit()][0])
        # print(chrg)
        mz = float(f[i+4][:-1].split("=")[1])
        # mass = mz * chrg
        # print(chrg, mz)
        p = i + 5
        
        ms2_dic = {}
        while p < len(f) and f[p][0].isdigit():
            if f[p].strip() == "END IONS":                        
                break
            else:
                mass = float(f[p].strip().split(" ")[0])
                insy = float(f[p].strip().split(" ")[1])
                if mass not in ms2_dic:
                    ms2_dic[mass] = insy
            p += 1
        # end_idx = p
        mgf_info_dic[title] = ms2_dic
        i = p + 1
    return mgf_info_dic

# print(read_spec("./test_spec.mgf"))