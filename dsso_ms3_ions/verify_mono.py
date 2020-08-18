# coding = utf-8


scans_old = set()
f = open('./mono1ppm.txt').readlines()
for line in f:
    scans_old.add(eval(line.split("\t")[0]))
f2 = open("./mono2ppm.txt").readlines()
for line in f2:
    scans_old.add(eval(line.split("\t")[0]))

print(len(scans_old))

scans_now = set()
f3 = open('./mono.txt').readlines()
for line in f3:
    scans_now.add(eval(line.split('\t')[0])[0])

print(len(scans_now))
ov = scans_now & scans_old
print(len(ov))
print(scans_old-scans_now)