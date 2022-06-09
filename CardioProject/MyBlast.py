from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
f1 = '/net/mraid08/export/genie/LabData/Data/MBPipeline/Strains/tmp/UZP/UZP/strain012_v2_fullrun.fastq'
results_df  = pd.DataFrame()


blastp_results = []
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
record_iterator = SeqIO.parse(f1, "fastq")
E_VALUE_THRESH=0.0000001

for i, record in enumerate(record_iterator):
    print(i, record.seq)
    if i>100:
        break
    result_handle=NCBIWWW.qblast("blastn", "nt",record.seq)
    blast_record=NCBIXML.read(result_handle)
    res_read = pd.Series()
    number_of_alignments =  0
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <E_VALUE_THRESH:
                number_of_alignments = number_of_alignments+1
                if alignment.title in res_read.index:
                    res_read.loc[alignment.title] = res_read.loc[alignment.title]+1
                else:
                    res_read.loc[alignment.title] = 1

    res_read =  res_read/float(number_of_alignments)
    results_df[i] = res_read
results_df.to_csv('/net/mraid08/export/genie/LabData/Data/MBPipeline/Strains/strain_012.csv')
#'gi|2120027897|gb|CP085087.1| Streptococcus suis strain Ssuis_MA2 chromosome, complete genome'


#f1 = '/net/mraid08/export/genie/LabData/Data/MBPipeline/Strains/tmp/UZP/UZP/strain010_v2_fullrun.fastq'





#entry = str(">" + i.description + "\n" + i.seq) f1 = open("test.txt", "w") f1.write(entry) f1.close() f2 = open("test.txt", "r") blastp_cline = NcbiblastpCommandline(query = 'test.txt', db = 'nr -remote', evalue = 0.05, outfmt = '7 sseqid evalue qcovs pident') res = blastp_cline() blastp_results.append(res) f2.close()