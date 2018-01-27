import os
os.chdir(r"D:\E\Collabaration\Jingji_20170120")

def caompare_two_value(str1, str2):
    num1 = int(str1)
    num2 = int(str2)
    if num1 == 0:
        return False
    else:
        if num2 == 0:
            if num1 - num2 >9 :
                return True
            else:
                return False
        else:
            if num1/num2 >= 2.0:
                return True
            else:
                return False


f = open("Jingji_group3.txt", 'r').readlines()
B = open("3_report.txt", "w")
B.write(f[0])

for line in f[1:]:
    line_list = line.rstrip("\n").split("\t")
    sample_4_spec = line_list[4]
    sample_3_spec = line_list[7]
    sample_2_spec = line_list[10]
    sample_1_spec = line_list[13]
    sample_list =[sample_4_spec, sample_3_spec, sample_2_spec, sample_1_spec]
    print(sample_list)
    for i in range(len(sample_list)):
        if sample_list[i] == "":
            sample_list[i] = "0"
    
    Bool_4_3 = caompare_two_value(sample_list[0], sample_list[1])
    Bool_4_2 = caompare_two_value(sample_list[0], sample_list[2])
    Bool_4_1 = caompare_two_value(sample_list[0], sample_list[3])

    bool_list = [Bool_4_3, Bool_4_2, Bool_4_1]
    if bool_list.count(True) > 2:
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
    sample_list =[sample_5_spec, sample_4_spec, sample_3_spec, sample_2_spec, sample_1_spec]
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
