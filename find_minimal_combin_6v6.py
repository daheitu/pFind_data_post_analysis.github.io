all_AA = "STYHRKEDQNGIVLAPMWF"
rep = open('min_combin_4.txt', 'w')


AA_resi_dic = {'A':71.037114,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185, \
'E':129.042593,'Q':128.058578,'G':57.021464,'H':137.058912,'I':113.084064, \
'L':113.084064,'K':128.094963,'M':131.040485,'F':147.068414,'P':97.052764, \
'S':87.032028,'T':101.047679,'U':150.95363,'W':186.079313,'Y':163.06332,'V':99.068414, \
'H2O':18.01056,'Proton':1.0072766}


def cal_pep_mass(seq):
    mass = AA_resi_dic["H2O"]
    for resi in list(seq):
        mass += AA_resi_dic[resi]
    return round(mass, 5)


# rep_dic = {}
for m in range(len(all_AA)):
    for n in range(m + 1, len(all_AA)):
        for o in range(n+1, len(all_AA)):
            site1_list = [all_AA[m], all_AA[n], all_AA[o]]
            for i in range(len(all_AA)):
                for j in range(i + 1, len(all_AA)):
                    site2_list = [all_AA[i], all_AA[j]]
                    for a in range(len(all_AA)):
                        for b in range(a + 1, len(all_AA)):
                            for c in range(b+1, len(all_AA)):
                                site3_list = [all_AA[a], all_AA[b], all_AA[c]]
                                for g in range(len(all_AA)):
                                    for h in range(g + 1, len(all_AA)):
                                        site4_list = [all_AA[g], all_AA[h]]
                                        mass_list = []
                                        for res1 in site1_list:
                                            for res2 in site2_list:
                                                for res3 in site3_list:
                                                    for res4 in site4_list:
                                                        mass_list.append(
                                                            cal_pep_mass(res1 + res2 + res3 + res4))
                                        mass_list.sort()
                                        delta_list = []
                                        for x in range(len(mass_list) - 1):
                                            delta_list.append(mass_list[x + 1] -
                                                            mass_list[x])
                                        pair = all_AA[m] + "/" +all_AA[n] + "/" + all_AA[o] + ", " + all_AA[i] + "/" +all_AA[j] +", "\
                                            + all_AA[a] + "/" + all_AA[b] + "/" + all_AA[c]+ ", "+ all_AA[g] + "/" + all_AA[h]
                                        min_delta = round(min(delta_list), 4)
                                        # print(pair, type(str(min(delta_list))))
                                        # rep_dic[pair] = str(min(delta_list))
                                        if min(delta_list) > 3:
                                            print(pair, str(min_delta))
                                            rep.write(pair + "\t" + str(min(delta_list)) + "\n")
                                   

rep.close()
