from ftplib import FTP
import sys
import os
import re

def ftpconnet(ftpserver,port,username,password):
    ftp = FTP()
    try:
        ftp.connect(ftpserver,port)
    except:
        raise IOError,'FTP connect failed!'

    try:
        ftp.login(username,password)
    except:
        raise IOError,'FTP login failed!'
    else:
        return ftp

def ftpdownload(ftp,ori_path,dest_path):
    for each in ftp.nlst(ori_path):
        try:
            ftp.cwd(each)
        except:
            filename = re.search('\S+\/(\S+)',each).group(1)
            local_file = dest_path + '/' + filename
            if os.path.exists(local_file):
                lsize = os.stat(local_file).st_size
                rsize = ftp.size(each)
                if lsize > rsize:
                    sys.stderr.write('the local file %s is bigger than the remote!\n'%local_file)
                    return False
                elif lsize == rsize:
                    sys.stderr.write('the file %s has been completed!\n'%local_file)
                bufsize = 1024 * 1024
                fp = open(local_file,'ab')
                ftp.retrbinary('RETR '+each,fp.write,bufsize,rest=lsize)
            else:
                bufsize = 1024 * 1024
                fp = open(local_file,'wb')
                ftp.retrbinary('RETR '+each,fp.write,bufsize)
        else:
            dirname = re.search('\S+\/(\S+)',each).group(1)
            dirname = dest_path + '/' + dirname + '/'
            os.system('mkdir -p %s'%dirname)
            ftpdownload(ftp,each,dirname)

def ftpclose(ftp):
    ftp.quit()

if __name__ == '__main__':
    ftp = ftpconnet('192.168.54.211',22,'ycao','cy13471118')
    print(ftp)
    #ftpdownload(ftp,'/pub/10.5524/100001_101000/100145/','./')
    ftpclose(ftp)