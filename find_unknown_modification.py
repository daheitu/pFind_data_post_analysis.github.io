import os
from math import log


def compare_num(num1, num2):
    num1 = int(num1)
    num2 = int(num2)
    if num1 < 10 and num2 < 10:
        return False
    else:
        if num1 == 0:
            return True
        elif num2 == 0:
            return True
        else:
            if abs(log(num1/num2, 2)) > 2:
                return True
            else:
                return False


def get_info(opened_file):
    modi_dic = {}
    for line in opened_file[3:]:
        if line[:3] == "---":
            break
        else:
            line_list = line.strip().split("\t")
            if len(line_list) == 25:
                peptide = line_list[1]
                modificaiton = line_list[6]
                spec_num = int(line_list[-1])
                if "PFIND_DELTA" in modificaiton:
                    modi_list = modificaiton.split(";")[:-1]
                    for modi in modi_list:
                        if "PFIND_DELTA" in modi:
                            modi_pos = modi.split(",")[0]
                            delta_mass = modi.split(",")[1].split("_")[-1]
                            modi_site = peptide[int(modi_pos) - 1]
                            if modi_site not in modi_dic:
                                modi_dic[modi_site] = {}
                                if delta_mass not in modi_dic[modi_site]:
                                    modi_dic[modi_site][delta_mass] = spec_num
                                else:
                                    pass
                            else:
                                if delta_mass not in modi_dic[modi_site]:
                                    modi_dic[modi_site][delta_mass] = spec_num
                                else:
                                    modi_dic[modi_site][delta_mass] += spec_num
                        else:
                            continue
                else:
                    continue
            else:
                continue

    print(modi_dic)

    b = open("report.txt", 'w')
    for AA in modi_dic:
        b.write(AA)
        b.write("\n")
        for mass in modi_dic[AA]:
            if modi_dic[AA][mass] > 4:
                b.write("\t".join([mass, str(modi_dic[AA][mass])]))
                b.write("\n")
            else:
                continue
    b.close()

    return modi_dic


os.chdir(r"D:\Program Files (x86)\pFind_build_20170814\GMDSD")
f1 = open("pFind.protein", 'r').readlines()
GSDMS_dic = get_info(f1)
os.chdir(r"D:\Program Files (x86)\pFind_build_20170814\CONTROL")
f2 = open(r"pFind.protein", 'r').readlines()
control_dic = get_info(f2)

print(len(GSDMS_dic), len(control_dic))

new_dic = {}
for AA in control_dic:
    comb_dic = {}
    mass_all_list = [] 
    for mass in control_dic[AA]:
        if mass not in mass_all_list:
            mass_all_list.append(mass)
        else:
            continue
    for mass in GSDMS_dic:
        if mass not in mass_all_list:
            mass_all_list.append(mass)
        else:
            continue
    for mass in mass_all_list:
        if mass in control_dic[AA]:
            if mass in GSDMS_dic[AA]:
                w_list = ["", mass, str(control_dic[AA][mass]), str(GSDMS_dic[AA][mass])]
            else:
                w_list = ["", mass, str(control_dic[AA][mass]), "0"]
        else:
            if mass in GSDMS_dic[AA]:
                w_list = ["", mass, "0", str(GSDMS_dic[AA][mass])]
        comb_dic[mass] = w_list
    new_dic[AA] = comb_dic
print(new_dic)

c = open("sig_site.txt", 'w')
for AA in new_dic:
    c.write(AA)
    c.write("\n")
    for mass in new_dic[AA]:
        sig_bool = compare_num(new_dic[AA][mass][-2], new_dic[AA][mass][-1])
        if sig_bool:
            # wr_list = new_dic[AA][mass].insert(0, "")
            c.write("\t".join(new_dic[AA][mass]))   
            c.write("\n")
    
        # c.write("\t".join(new_dic[AA][mass]))
        # c.write("\n")
    # sig_bool = compare_num(new_dic[AA][-2], new_dic[AA][-1])
    # if sig_bool:
    #     c.write("\t".join(new_dic[AA]))
    #     c.write("\n")
c.close()
print("done")
