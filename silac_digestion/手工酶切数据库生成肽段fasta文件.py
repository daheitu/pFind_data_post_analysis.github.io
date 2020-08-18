#-*-coding:utf8 -*-
import os
import sys
import time
import string
import platform


def main():
    #fastaPath = r"GST_sequence.fasta" #modify: fasta location
    fastaPath = r"./Integrin.fasta"

    file_fasta = open(fastaPath)

    output = open(fastaPath+"_enzyme_chymotryp_tryp", 'w')

    file_row = file_fasta.readlines()

    Miss_Cleave = 1 #modify: miss cleavage no.

    ProStart = 0
    ProEnd = 0
    for i in range(0, len(file_row)): 
        if file_row[i][0] == ">":
            ProEnd = i-1
            if ProEnd > ProStart:
                [sequence, position] = CutPro(ProStart, ProEnd, file_row)
                OutputFasta(sequence, position, output, Miss_Cleave)
                ProStart = i
        if i == len(file_row)-1 :
            ProEnd = i
            if ProEnd > ProStart:
                [sequence, position] = CutPro(ProStart, ProEnd, file_row)
                OutputFasta(sequence, position, output, Miss_Cleave)

    file_fasta.close()
    output.close()
    return

def CutPro(ProStart, ProEnd, file_row): 
    AC = file_row[ProStart].split()[0]
    sequence = []
    position = []
    s=''
    for i in range(ProStart+1, ProEnd+1): 
        s = s + file_row[i].rstrip()

    currtPos = 0
    for i in range(0, len(s)):
        if s[i] in ["F","Y", "W", "M", "L", "K", "R"]:#modify：add all the C-terminal cleave AA here
            if currtPos <= i:
                position.append(AC + "_" + str(currtPos+1) + "-" + str(i+1))#加一以从一开始记录
                sequence.append(s[currtPos:i+1])
                currtPos = i+1
        elif s[i] == 'U':#modify：add all the N-terminal cleave AA here
            if currtPos <= i-1:
                position.append(AC + "_" + str(currtPos+1) + "-" + str(i))
                sequence.append(s[currtPos:i])
                currtPos = i
        elif i == len(s)-1:
            position.append(AC + "_" + str(currtPos+1) + "-" + str(i+1))
            sequence.append(s[currtPos:i+1])

    return [sequence, position]

''' currtMiss = 0
    currtPos = 0
    for m in range(1, Miss_Cleave):
        for i in range(0, len(s)):
            if s[i] == 'K' or s[i] == 'R':
                if currtMiss < m:
                    currtMiss+=1
                else:
                    if currtPos <= i:
                        position.append(AC + "_" + str(currtPos+1) + "-" + str(i+1))#加一以从一开始记录
                        sequence.append(s[currtPos:i+1])
                        currtPos = i+1
                        currtMiss = 0
            elif s[i] == 'D' or s[i] == 'E':
                if currtMiss < m:
                    currtMiss+=1
                else:
                    if currtPos <= i-1:
                        position.append(AC + "_" + str(currtPos+1) + "-" + str(i))
                        sequence.append(s[currtPos:i])
                        currtPos = i
                        currtMiss = 0
            elif i == len(s)-1 and currtMiss == m:
                position.append(AC + "_" + str(currtPos+1) + "-" + str(i+1))
                sequence.append(s[currtPos:i+1])'''

def OutputFasta(sequence, position, output, Miss_Cleave):
    for i in range(0, len(sequence)):
        output.write("%s\n" % position[i]) 
        output.write("%s\n\n" % sequence[i])

    if Miss_Cleave < 1:
        return

    for m in range(1, Miss_Cleave+1):
        for i in range(0, len(sequence)-m):
            tempPos = joinpos(position[i], position[i+m])
            output.write("%s\n" % tempPos)
            for j in range(i, i+m+1):
                output.write("%s" % sequence[j])
            output.write("\n\n")

    return

def joinpos(pos1, pos2):
    temp1 = pos1[0:pos1.rfind('-')]
    temp2 = pos2[pos2.rfind('-')+1:]
    return temp1+'-'+temp2

if __name__ == '__main__':
	main()