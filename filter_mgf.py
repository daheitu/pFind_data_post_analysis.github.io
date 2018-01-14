import os


def generate_ion_mass_range(num):
    deta = num * 0.5 / 100000
    return num - deta, num + deta


def words_frequcecy(list):
    mass_frequcy_dic = {}
    for word in list:
        if word not in mass_frequcy_dic:
            mass_frequcy_dic[word] = 0
        mass_frequcy_dic[word] += 1
    #x,y : operator.lt(x[1],y[1]), reverse=True)
    return sorted(mass_frequcy_dic.items(), key=lambda d: d[1], reverse=True)


os.chdir(r"E:\TEST_AB+N\B_N")
file_list = os.listdir(os.getcwd())
for name in file_list:
    if name[-3:] == "mgf":
        f = open(name)
        all = f.readlines()
        ion_be = []
        ion_end = []
        for i in range(len(all)):
            if all[i].strip() == "BEGIN IONS":
                ion_be.append(i)
            elif all[i].strip() == "END IONS":
                ion_end.append(i)
        print(len(ion_be))
        if len(ion_be) != len(ion_end):
            print("erro")
        else:
            report_file_name = name[:name.find(".mgf")] + ".txt"
            b = open(report_file_name, 'w')
            b.write("\t".join(["title", "charge", "m/z", "mass"]))
            b.write("\n")
            for k in ion_be:
                title = all[k + 1].strip()
                charge = all[k + 2][7:9]
                mass_over_z = all[k + 4][all[k + 4].find("=") + 1:all[k + 4]
                                         .find("\n")]
                mass = float(mass_over_z) * int(charge[0]) - int(charge[0])
                w_list = [title, charge, mass_over_z, str(mass)]
                ms2_ion_list = []
                #print k+5,ion_end[ion_be.index(k)]
                for m in range(k + 5, ion_end[ion_be.index(k)]):
                    #if all[m].strip().split(" ")[0].isdigit():
                    ms2_ion_list.append(all[m].strip().split(" ")[0])
                #print ms2_ion_list
                [y1_low, y1_up] = generate_ion_mass_range(147.1128)
                [y2_low, y2_up] = generate_ion_mass_range(244.1656)
                [y3_low, y3_up] = generate_ion_mass_range(391.2340)
                [b2_low, b2_up] = generate_ion_mass_range(201.1234)
                include_y1 = False
                include_y2 = False
                include_y3 = False
                include_b2 = False
                for ion in ms2_ion_list:
                    if float(ion) > y1_low and float(ion) < y1_up:
                        include_y1 = True
                #print include_y1
                    elif float(ion) > b2_low and float(ion) < b2_up:
                        include_b2 = True
                    elif float(ion) > y2_low and float(ion) < y2_up:
                        include_y2 = True
                    elif float(ion) > y3_low and float(ion) < y3_up:
                        include_y3 = True
                if include_y1 and include_y2 and include_y3:
                    b.write("\t".join(w_list))
                    b.write("\n")
            b.close()
        b = open(report_file_name, 'a+')
        b.seek(0)
        print(len(b.readlines()))
        mass_list = []
        table = open(report_file_name).readlines()
        for i in range(1, len(table)):
            mass_list.append(round(float(table[i].split("\t")[3]), 1))

        mass_frequcy_dic = words_frequcecy(mass_list)
        print(mass_frequcy_dic)
        for key in mass_frequcy_dic:
            b.write("\t".join([str(key[0]), str(key[1])]))
            b.write("\n")
        b.close()
        f.close()
