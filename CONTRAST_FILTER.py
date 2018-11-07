import os

os.chdir(r"D:\E\Collabaration\Jingji20180318\Group6")


def caompare_two_value(str1, str2):
    num1 = int(str1)
    num2 = int(str2)
    if num1 < num2:
        num1, num2 = num2, num1
    else:
        pass

    if num1 == 0:
        return False
    else:
        if num2 == 0:
            if num1 - num2 > 9:
                return True
            else:
                return False
        else:
            if num1 > 10 and num1 / num2 >= 2.0:
                return True
            else:
                return False


f = open("result.txt", 'r').readlines()
B = open("G6_report.txt", "w")
B.write(f[0])

for line in f[1:]:
    line_list = line.rstrip("\n").split("\t")
    sample_wt_spec = line_list[4]
    sample_g17v_spec = line_list[7]
    # sample_2_spec = line_list[10]
    # sample_1_spec = line_list[13]
    sample_list = [sample_wt_spec, sample_g17v_spec]
    # print(sample_list)
    for i in range(len(sample_list)):
        if sample_list[i] == "":
            sample_list[i] = "0"

    Bool_wt_g17v = caompare_two_value(sample_list[0], sample_list[1])
    # Bool_4_1 = caompare_two_value(sample_list[0], sample_list[3])

    if Bool_wt_g17v is True:
        B.write(line)
B.close()
"""
f = open("Jingji_group_1.txt", 'r').readlines()
B = open("1_and_report.txt", "w")
B.write(f[0])
for line in f[1:]:
    line_list = line.rstrip("\n").split("\t")
    sample_1_spec = line_list[16]
    sample_5_spec = line_list[4]
    sample_4_spec = line_list[7]
    sample_3_spec = line_list[10]
    sample_2_spec = line_list[13]
    sample_list =[sample_5_spec, sample_4_spec, sample_3_spec, 
        sample_2_spec, sample_1_spec]
    print(sample_list)
    for i in range(len(sample_list)):
        if sample_list[i] == "":
            sample_list[i] = "0"
    Bool_5_4 = caompare_two_value(sample_list[0], sample_list[1])
    Bool_5_1 = caompare_two_value(sample_list[0], sample_list[4])
    Bool_3_2 = caompare_two_value(sample_list[2], sample_list[1])
    Bool_3_1 = caompare_two_value(sample_list[2], sample_list[3])

    bool_list_5 = [Bool_5_4, Bool_5_1]
    bool_list_3 = [Bool_3_2, Bool_3_1]
    if bool_list_5.count(True) == 2 and bool_list_3.count(True) == 2:
        B.write(line)
B.close()
"""
