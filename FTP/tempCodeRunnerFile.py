print(name)
        path = os.path.join(ftpIP, name)
        print(path)
        cmd = "wget " + path
        print(cmd)