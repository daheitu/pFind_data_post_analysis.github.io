这个脚本可以自动调用pParse和searcher程序，进行导出和搜索
# 4个重要的参数
1. raw_path raw文件所在的文件夹的路径
2. db_name 搜索使用的数据库名称，必须在db.ini中存在的，而且名称与之一致
3. plink_bin_path pLink的安装目录bin所在路径
4. plink_para_demo .plink参数文件的模板文件，根据它生成运行所需的参数问
# 使用方法
根据需要将4个参数填完成，保存脚本，使用python +脚本路径运行