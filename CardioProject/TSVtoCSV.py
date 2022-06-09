import re
from os import listdir
from os.path import isfile, join


## use the following for multiple files conversion:

# onlyfiles = [f for f in listdir("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/") if isfile(join("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/", f))]
# onlyfiles = [datafile for datafile in onlyfiles if datafile.startswith ('HIP') and datafile.endswith('.tsv')]
# print onlyfiles
# countfiles=0
# for datafile in onlyfiles[51::]:
#     tsv = open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s' % datafile, 'r')
#     fileContent =  tsv.read()
#     fileContent = re.sub("""(?ism)(,|"|')""", r"\\\1", fileContent) # escape all especial charaters (" ' ,) rfc4180
#     fileContent = re.sub("\t", ",", fileContent) # convert from tab to comma
#     csv_file = open("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles/%s.csv" % sample_name, "w")
#     csv_file.write(fileContent)
#     csv_file.close()
    


###use the following to convert one file:   
tsv = open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/HIP14060.tsv', 'r')
fileContent =  tsv.read()
fileContent = re.sub("""(?ism)(,|"|')""", r"\\\1", fileContent) # escape all especial charaters (" ' ,) rfc4180
fileContent = re.sub("\t", ",", fileContent) # convert from tab to comma
csv_file = open("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles/HIP14060.csv", "w")
csv_file.write(fileContent)
csv_file.close()
    