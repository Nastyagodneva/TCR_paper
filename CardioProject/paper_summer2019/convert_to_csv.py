import pandas as pd
from os import listdir, makedirs
from os.path import isdir


############################################3
#parameters to change
dir_location = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Presentations and Manuscripts/CardioTCR paper/\
June2019/fig3_balanced_comparison/calc_fig3/'
dir_name = 'excel_files_new/'
#############################################



dir1=dir_location + dir_name

csv_dir_name = dir1.replace('excel_files_new/','csv_files_new/')
print (csv_dir_name)
if not isdir(csv_dir_name):
    makedirs(csv_dir_name)

print ('converting %s files from dir %s' %(len(listdir(dir1)),dir1))
for f in listdir(dir1):
    print f
    df = pd.read_excel(dir1+f)
    print df.head()
    new_f_name = f.replace('.xlsx', '.csv')
    print new_f_name
    print csv_dir_name+new_f_name
    df.to_csv(csv_dir_name+new_f_name)
    # except:
    #     print ('couldnt convert this file\n')
