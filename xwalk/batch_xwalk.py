import os
import time
os.chdir(r"D:\Program Files (x86)\xwalk_v0.6\bin")

for root, dirs, files in os.walk(r"G:\merged_complex"):
    for file in files:
        if file[-4:] == ".pdb":
            print(os.path.join(root, file))
            output_path = r'E:\Script\xwalk_output\KR'
            output_name = file[:-4] + "_RK.txt"
            input_file = os.path.join(root, file)
            output_file = os.path.join(output_path, output_name)
            cmd_list = ["java -Xmx4096m Xwalk -infile ", input_file, " -out ", output_file, \
                " -aa1 ARG -aa2 LYS -a1 CA -a2 CA -max 60 -bb -homo -inter"]
            cmd = "".join(cmd_list)
            print(cmd)
            os.system(cmd)
            #time.sleep(5)

            