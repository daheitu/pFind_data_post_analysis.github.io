import os

os.chdir(r"G:\Huronggui\test\test")
deta_ppm = 20
feature_ion_mass_list = [147.1128, 187.0713, 204.1342, 275.1714, 422.2398] # 126.1282 114.128
insy_cutoff = 0.15


def generate_ion_mass_range(num):
    deta = num * deta_ppm / 1000000
    return num - deta, num + deta


def judge_spectra(feature_list, ms2_dic):
    count_list = [0] * len(feature_list)
    mass_list = list(ms2_dic.keys())
    max_insy = max(list(ms2_dic.values()))
    high_feature_list = []
    mass_range_dic = {}
    for mass in feature_list:
        low_mass, high_mass = generate_ion_mass_range(mass)
        high_feature_list.append(high_mass)
        mass_range_dic[mass] = [low_mass, high_mass]
    
    ft_idx = 0
    ms2_idx = 0
    while ft_idx < len(feature_list) and ms2_idx < len(mass_list):
        if mass_list[ms2_idx] > high_feature_list[ft_idx]:
            ft_idx += 1
        elif mass_list[ms2_idx] >= mass_range_dic[feature_list[ft_idx]][0]:
            ft_idx += 1
            ms2_idx += 1
            
            if ms2_dic[mass_list[ms2_idx -1]]/max_insy >= insy_cutoff:
                count_list[ft_idx -1] = 1
            else:
                pass
        else:
            ms2_idx += 1
    if count_list.count(0):
        return False
    else:
        return True


def main():
    fl_list = os.listdir(os.getcwd())
    for fl in fl_list:
        if fl[-9:] != "HCDFT.mgf":
            continue
        else:
            print(fl)
            total_spec = 0
            ft_spec = 0
            f = open(fl, 'r').readlines()
            rep_name = fl[:-4] + "_filter.mgf"
            b = open(rep_name, 'w')
            i = 0
            while i < len(f):
                print("current is " + str(i))
                if f[i].strip() == "BEGIN IONS":
                    begin_idx = i
                    total_spec += 1
                else:
                    print("wrong")
                chrg = int([i for i in f[i+2] if i.isdigit()][0])
                # print(chrg)
                mz = float(f[i+4][:-1].split("=")[1])
                # mass = mz * chrg
                print(chrg, mz)
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
                end_idx = p

                if judge_spectra(feature_ion_mass_list, ms2_dic):
                    ft_spec += 1
                    for m in range(begin_idx, end_idx + 1):
                        b.write(f[m])
                else:
                    pass
                
                i = p + 1
            print(round(float(ft_spec)/total_spec, 10))
            b.close()


if __name__ == "__main__":
    main()           
