import os

os.chdir(r"E:\MASS SEPTRA DATA\20170310\H2")
f = open("H2.txt", 'r').readlines()
delta_list = [1.0, 2.0, 3.0]
deta_Da = 0.05
b = open("rep.txt", 'w')


def judge_num(num1, num2):
    if num1 > num2 and num1 - num2 < deta_Da:
        return True
    else:
        return False


title_pre = f[1].strip().split("\t")[0][:-4]
mass_pre = float(f[1].strip().split("\t")[3])
k = 2
while k < len(f):
    line_list = f[k].strip().split("\t")
    if len(line_list) == 4:
        title_post = line_list[0][:-5]
        mass_post = float(line_list[3])
        if title_post == title_pre:
            isotop_bool = False
            for ms in delta_list:
                if judge_num(mass_post - ms, mass_pre):
                    isotop_bool = True
                    break
                else:
                    continue
            print(isotop_bool)
            if isotop_bool is False:
                b.write(f[k])
            else:
                pass
        else:
            b.write(f[k])
        title_pre = title_post
        mass_pre = mass_post
        k += 1
    else:
        break

b.close()