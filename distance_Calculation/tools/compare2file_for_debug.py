#比较两个的行是否一致，如果不一致，则输出该行
import os
os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2020.01.16.11.01.13_GST5_025UL\reports")
ori_file_path = r"./xlink_distance.txt"
new_file_path = r"./xlink_distance.txt1"

def get_fl_lines(fl_path):

    # for fl in os.listdir(fl_path):
    #     if fl.endswith(".csv"):
    #         fpath = os.path.join(fl_path, fl)
    return open(fl_path).readlines()


ori = get_fl_lines(ori_file_path)
new = get_fl_lines(new_file_path)

if len(ori) != len(new):
    print("文件长度不一致")
else:
    for i in range(len(ori)):
        if ori[i] != new[i]:
            print("the line %d is not equel" % (i+1))
            print(ori[i])
            print(new[i])
            break
    print("Right")
    