import os
import re
os.chdir(r"D:\Onedriver\OneDrive\文档")
all_AA = "STYHRKEDQNGIVLAPMWF"
f = open('./peptides_sythensis_sgc.csv', 'r').readlines()


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


def find_fix2_list(site1_list, site2_list,rep):
    report_dic = {}
    max_combine = []
    for m in range(len(all_AA)):
        for n in range(m + 1, len(all_AA)):
            site3_list = [all_AA[m], all_AA[n]]
            for a in range(len(all_AA)):
                for b in range(len(all_AA)):
                    site4_list = [all_AA[a], all_AA[b]]
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
                        delta_list.append(mass_list[x + 1] - mass_list[x])
                    pair = all_AA[m] + "/" + all_AA[n] + ", " + all_AA[
                        a] + "/" + all_AA[b]
                    min_delta = round(min(delta_list), 4)
                    if max_combine == []:
                        max_combine = [min_delta, [pair, min_delta]]
                    else:
                        if min_delta > max_combine[0]:
                            max_combine = [min_delta, [pair, min_delta]]
                    if min(delta_list) > 3:
                        # print(pair, str(min_delta))
                        report_dic[pair] = min_delta
                        #rep.write(pair + "\t" + str(min(delta_list)) + "\n")
    print(max_combine)
    new_list = sorted(report_dic.items(), key=lambda d: d[1], reverse=True)
    if len(new_list) < 20:
        for sit in new_list:
            rep.write("\t".join([str(ele) for ele in list(sit)]) + "\n")
    else:
        for sit in new_list[:20]:
            rep.write("\t".join([str(ele) for ele in list(sit)]) + "\n")


def find_fix1_list(site1_list, rep):
    report_dic = {}
    max_combine = []
    for m in range(len(all_AA)):
        for n in range(m + 1, len(all_AA)):
            site3_list = [all_AA[m], all_AA[n]]
            for a in range(len(all_AA)):
                for b in range(len(all_AA)):
                    site4_list = [all_AA[a], all_AA[b]]
                    mass_list = []
                    for res1 in site1_list:
                        for res3 in site3_list:
                            for res4 in site4_list:
                                mass_list.append(cal_pep_mass(res1 + res3 + res4))
                    mass_list.sort()
                    delta_list = []
                    for x in range(len(mass_list) - 1):
                        delta_list.append(mass_list[x + 1] - mass_list[x])
                    pair = all_AA[m] + "/" + all_AA[n] + ", " + all_AA[
                        a] + "/" + all_AA[b]
                    min_delta = round(min(delta_list), 4)
                    if max_combine == []:
                        max_combine = [min_delta, [pair, min_delta]]
                    else:
                        if min_delta > max_combine[0]:
                            max_combine = [min_delta, [pair, min_delta]]
                    if min(delta_list) > 3:
                        report_dic[pair] = min_delta
    print(max_combine)
    new_list = sorted(report_dic.items(), key=lambda d: d[1], reverse=True)
    if len(new_list) < 20:
        for sit in new_list:
            rep.write("\t".join([str(ele) for ele in list(sit)]) + "\n")
    else:
        for sit in new_list[:20]:
            rep.write("\t".join([str(ele) for ele in list(sit)]) + "\n")


def main():
    rep = open("report_best.txt", 'w')
    for line in f[1:]:
        line_list = line.strip().split(',')
        pep = line_list[4]
        name = line_list[1]
        rep.write("\t".join(line_list) + "\n")
        var_group_list = re.findall("\((.*?)\)", pep)
        if len(var_group_list) == 1:
            site1_list = var_group_list[0].split("/")
            find_fix1_list(site1_list, rep)
        elif len(var_group_list) == 2:
            site1_list = var_group_list[0].split("/")
            site2_list = var_group_list[1].split("/")
            find_fix2_list(site1_list, site2_list, rep)
    rep.close()


if __name__ == "__main__":
    main()