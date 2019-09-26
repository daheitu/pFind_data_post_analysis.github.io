# coding = utf-8
import os
path_pk = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.24.21.37.35_Integrin_PK\reports"
path_LT = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.21.15.42.25_INTEgrin_LT\reports"
path_LTA = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.21.21.22.17_Integrin_LTA\reports"
path_LTG = r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.24.21.18.43_Integrin_LTG\reports"
#
for path in [path_LT, path_LTA, path_LTG, path_pk]:
    enzyme = path.split("\\")[5].split("_")[-1]
    flpath = os.path.join(path, "report.csv")
    f = open(flpath, 'r').readlines()
    reportFl = "Integrin_aiib_" + enzyme + ".csv"
    b = open(reportFl, 'w')
    b.write(f[0])
    for line in f[1:]:
        if "sp|P08514|" in line:
            b.write(line)
    b.close()



#os.chdir(r"C:\Users\Yong Cao\Documents\pLink\pLink_task_2019.09.24.21.37.35_Integrin_PK\reports")
"""
f = open("report.csv", 'r').readlines()

b = open("Integrin_aiib.csv", 'w')
b.write(f[0])
for line in f[1:]:
    if "sp|P08514|" in line:
        b.write(line)

b.close()
"""