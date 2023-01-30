import os
import sys
import webbrowser as web

#--------key parameters----------------
#pFind结果文件路径
pfind_cfg = r'E:\workspace\pFindTask114_YANZ_trans_id\param\pFind.cfg'
#目标蛋白质AC
target_pro = r'tr|A8J6H7|A8J6H7_CHLRE'
#输出HTML文件，程序运行成功后会自动打开
output_html_path = r'E:\workspace\pFindTask114_YANZ_trans_id\param\ppmap.html'

#--------advanced parameters-----------
qv_limit = 0.01 #FDR阈值
unique_pep = False #是否只考虑unique peptide，True/False，默认为False
segment_per_line = 5 #每行多少个氨基酸片段
aa_per_segment = 10 #每个氨基酸片段多少个氨基酸


with open(pfind_cfg, 'r') as f:
    lines = f.readlines()

fasta_path = r''
result_path = r''
for line in lines:
    if line.startswith('fastapath='):
        fasta_path = line.rstrip().split('=')[1]
    elif line.startswith('outputpath='):
        result_path = line.rstrip().split('=')[1] + os.sep + 'pFind.spectra'

if not (fasta_path and result_path):
    print('error!', fasta_path, result_path)
    exit(0)

proteins = []

#-----read FASTA file-----------
f = open(fasta_path, 'r')

while True:
    line = f.readline()
    if not line:
        break
    line = line.rstrip("\n")
    if line.startswith('>'):
        proteins.append(['', '', ''])
        if " " in line:
            proteins[-1][0], proteins[-1][1] = line[1:].split(maxsplit=1)
        else:
            proteins[-1][0] = line[1:].strip
    else:
        proteins[-1][2] += line

f.close()
print('proteins in fasta:', len(proteins))

ac, de, sq = '', '', ''
for pro in proteins:
    if pro[0] == target_pro:
        ac, de, sq = pro[0], pro[1], pro[2]

#------read result file and find valid PSMs
peptides = {}
f = open(result_path, 'r')
line = f.readline()
while True:
    line = f.readline()
    if not line:
        break
    line = line.rstrip()
    segs = line.split('\t')
    if float(segs[4]) > qv_limit:
        break
    pros = segs[12].split('/')[:-1]
    if target_pro not in pros:
        continue
    peptides[segs[5]] = 1
f.close()
print('peptides matched to %s'%(target_pro), len(peptides))

for i,pep in enumerate(peptides):
    print(i, pep)
    
v = ['0'] * len(sq)
sq_tmp = sq.replace('I', 'L')
for pep in peptides:
    pos = sq_tmp.find(pep.replace('I', 'L'))
    if pos == -1:
        print('error!', 'peptide:', pep, 'protein:', target_pro)
        exit(0)
    for i in range(pos, pos + len(pep)):
        v[i] = '1'

v = ''.join(v)


#------generate html---------
line_id_width = len(str(len(sq)))
out_str = '<html>\n<head>\n<title>peptide matched to %s</title>\n</head>\n<body>\n<b>\n<font face="Courier New" size="5">%s %s<br><br>'%(target_pro, ac, de)

start_pos = 0
start_seg = 0
while start_pos < len(sq):
    if start_seg == 0:
        for j in range(line_id_width-len(str(start_pos+1))):
            out_str += '&nbsp;'
        out_str += str(start_pos+1) + "\t"

    st = start_pos
    ed = start_pos + aa_per_segment
    while st < ed:
        pos = v.find('1', st)
        if pos == -1 or pos >= ed:
            out_str += sq[st:ed]
            break
        else:
            out_str += sq[st:pos]
            out_str += '<font color="#FF0000">'
            pos2 = v.find('0', pos)
            if pos2 == -1 or pos2 >= ed:
                out_str += sq[pos:ed]
            else:
                out_str += sq[pos:pos2]
            out_str += '</font>'
            st = pos2
    start_pos += aa_per_segment
    start_seg = (start_seg + 1)%segment_per_line
    if not start_seg:
        out_str += '\n<br>'
    else:
        out_str += '&nbsp;'


out_str += '</font></b>\n</body>\n</html>'

with open(output_html_path, 'w') as f:
    f.write(out_str)

web.open_new(output_html_path)
