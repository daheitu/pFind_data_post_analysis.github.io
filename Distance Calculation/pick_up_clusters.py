import os
# os.chdir(r"E:\Rosetta Learning\prepack\test_for_pickTop")
"""
os.chdir(r"E:\Rosetta Learning\prepack\test_for_pickTop")

cluInfo = open(r"RMSD-matrix-top-200-merge_sort_RMSD_filter_24_3.fasc_cluster_2.5", 'r').readlines()
scorefl = open("merge_sort_RMSD_filter_24_3.fasc", 'r').readlines()

def tgtConfInfo(cluInfo):
    repDic = {}
    for line in cluInfo:
        print(line)
        confName, clusNum = line.strip().split('\t')
        confName = confName[confName.find("_")+1:]
        repDic[confName] = [confName, clusNum]
    return repDic

def pickupConf(scorefl, )
print(tgtConfInfo(cluInfo))
#print(scorefl[2])
"""
def creatFlag(pdbName):
    flag = """-s ./%s
-partners AB_C
-ex1
-ex2aro
-dock_pert 3 8
-native ../../3u28_change.pdb
-norepack1
-docking:sc_min
-dock_rtmin
-mute core.util.prof ## dont show timing info
-out:overwrite
-constraints:cst_file ../../lowRes_cng_DSS
-out:file:scorefile %s_Local_refine.fasc
-out:file:fullatom #output in fullato scorefile
-mute core.io.database
#-out:prefix %s_
-out:nstruct 3
-out:path:pdb  ./
-out:path:score ./
""" % (pdbName, pdbName[:-4], pdbName[-4:])
    return flag

print(creatFlag("3u28"))


def writeFlag(pdbName, ppath):
    flag = creatFlag(pdbName)
    flagFlName = "flag_local_refine"
    flagPath = os.path.join(ppath, flagFlName)
    fg = open(flagPath, 'w')
    fg.write(flag)
    fg.close()


def mktargetPDB(confName, ppath):
    pdbName = confName+'.pdb'
    pdb_path = os.path.join(ppath, pdbName)
    b = open(pdb_path, 'w')
    commChainPath = '/nibs/home/ycao/Rosetta/common_chain/3u28_relaxed_AB.pdb'
    commChain = open(commChainPath, 'r').readlines()
    variChain = open('./top_2000_pose_filter/'+confName+".remove_chainABC.pdb")
    for line in commChain:
        b.write(line)
    for line in variChain:
        if line[:4] == "ATOM":
            b.write(line)
    b.write("TER"+'\n')
    b.close()


def localRefine(confName):
    path = os.path.join("local2refine" , confName)
    print(path)
    os.system("mkdir " + path)
    mktargetPDB(confName, path)
    writeFlag(confName+'.pdb', path)
    print("cd "+ path)
    #os.system("cd "+ path)
    os.chdir(path)
    dockCMD = "~/software/tools/Rosetta/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/docking_protocol.static.linuxgccrelease"
    print(dockCMD +  " @flag_local_refine")
    os.system("nohup " + dockCMD +  " @flag_local_refine")
    #os.system("cd ../../")
    os.chdir(r"/nibs/home/ycao/CNGP_medol/Rosetta_39/DSS")


def main():
    topClusNum = 3
    f = open('./matrix/RMSD-matrix-top-200-merge_sort_RMSD_filter_24_3.fasc_cluster_10A-report' ,'r').readlines()
    os.system("mkdir local2refine")
    for i in range(3):
        lineList = f[i].strip().split("\t")
        confmName = lineList[2][lineList[2].find("_")+1:]
        score = lineList[3]
        localRefine(confmName)

if __name__ == "__main__":
    main()

