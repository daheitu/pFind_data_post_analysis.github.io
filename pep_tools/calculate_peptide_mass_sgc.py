pep = 'D(V/L/A)KPSTEHI(H/W)LK'

import re

a = re.findall("\((.*?)\)", pep)
b = re.split("\((.*?)\)", pep)

a_len = len(a):
if a_len == 0:
    print("wrong")
elif a_len == 1:
    vair_1_list = a[0].split('/')
    vair_idx = b.index(a[0])
    for char in vair_1_list:
        new_list = 

elif a_len == 2:
    vair_1_list = a[0].split('/')
    vair_2_list = a[1].split('/')
else:
    print("please contact Yong Cao")

vai_dic = {}
for wd in b:
    if '/' not in wd:
        wd_idex = b.index(wd)
        vai_dic[index] = wd.split('/')
    else:
        new_str += wd 
print(new_str)

[x*x x in ]