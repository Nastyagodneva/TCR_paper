'''
This script is aimed to answer the final exam in the 'Python for Genomic Data Science'
course in coursera
'''

from Bio import SeqIO
import pandas as pd

'''
Q1:
find how many sequences are in the file: one method is to use the fasta parser from BioPython and check the
length of the sequences list. the other method is to use the string method 'count' to count the number of '>' in the file
'''

f1=open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Coursera/PythonIntro/dna2.fasta','r')
seq_pars=SeqIO.parse(f1, "fasta") ## using a module of biopython that knows to parse different file formats. 
                                  ## in case of fasta, it generates a list of entries, each contains the seq ID 
                                  ##and the sequence itslef
seq_list=list(seq_pars)
print('the number of sequences is %s' % len(seq_list)) #counts the number of entries after parsing
f1.close

## the following paragraph validates the result with another strategy. it reads the content from the fasta file
## as text, and counts the number of '>' signs, which should appear in the negining of each sequence header in the 
## fasta format

f2=open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Coursera/PythonIntro/dna2.fasta','r')
fasta_content=f2.read()
print('the number of times a ">" occures in the file: %s' % fasta_content.count('>'))
f2.close


'''
Q2:
What are the lengths of the sequences in the file? What is the longest sequence and what is the 
shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers?

Again0use the FASTA parser from biopython (SeqIO.parse), for each sequence it generates an entry that contains
the sequence ID (id) and the sequence itself (seq). therefore, the length of seq is the sequence length.

'''

f1=open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Coursera/PythonIntro/dna2.fasta','r')
ID=[]
sequence=[]
length=[]
for seq_record in SeqIO.parse(f1, "fasta"):
    ID.append(seq_record.id)
    sequence.append(repr(seq_record.seq))
    length.append(len(seq_record))
    
max1=max(length)
min1=min(length)
n_max=length.count(max1)
n_min=length.count(min1)

print('the maximal sequence length is %s bps' %max1)
print('the minimal sequence length is %s bps' %min1)
print('the number of sequences with the maximal length is %s' %n_max)
print('the number of sequences with the minimal length is %s' %n_min)

min1_id=[]
max1_id=[]
for i in range(len(ID)):
    if length[i]==max1:
        max1_id.append(ID[i])
    elif length[i]==min1:
        min1_id.append(ID[i])
print ('the longest sequences are: %s' %max1_id)
print ('the shortest sequences are: %s' %min1_id)

'''
Q3:
Given an input reading frame on the forward strand (1, 2, or 3) your program should be able to identify 
all ORFs present in each sequence of the FASTA file, and answer the following questions: what is the length 
of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? 
For a given sequence identifier, what is the longest ORF contained in the sequence represented by that 
identifier? What is the starting position of the longest ORF in the sequence that contains it? The position 
should indicate the character number in the sequence. 
***

the following code was adopted from the biopython tutorial: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc292
***there is some mistake in the code, the results are wrong! need to correct before using!***

'''



f1=open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Coursera/PythonIntro/dna2.fasta','r')
table = 1 #standard codon table. can be changed to accomodate the bacterial translation table (11) or any 
          # other codon table, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
          

# min_pro_len = 100 ## I din't use this parameter


## the following function returns all ORF translated, with their length, start and stop sites
def find_orfs_with_trans(seq, trans_table, frame_num):
    
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq)]: ##changed the code here to test only the forward strand
                                    ## if both strands are needed, go back to the original code.
        for frame in range(3):
            while frame==frame_num-1:
                trans = str(nuc[frame:].translate(trans_table))
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                    trans[aa_start:aa_end]))
                    aa_start = aa_end+1
                else:
                    break
    answer.sort()
    return answer

## now generate lists that contains all the IDs in the file, a list of ORF lengths for eac ID, a list of
## ORF start sites for each ID, and the length of the longest ORF in each ID"

ID_list=[]
length_list=[]
start_site_list=[]
max_length_list=[]
frame_num=1
for seq_record in SeqIO.parse(f1, "fasta"):   
#     print seq_record.id 
    orf_list = find_orfs_with_trans(seq_record.seq, table,frame_num)
    ## generate lists of ORF lengths and start sites within the specific ID:
    ORF_lengths=[]
    ORF_start_sites=[]
    
    for start, end, strand, pro in orf_list:
            ORF_lengths.append((len(pro)+1)*3)
            ORF_start_sites.append(start)
#         print("%s...%s - length %i, strand %i, %i:%i" \
#               % (pro[:30], pro[-3:], len(pro), strand, start, end))
    
    max_length_for_ID=max(ORF_lengths) #Calculate the maximal length of ORF withing a specific ID
    ##generate lists with all IDs in the file and their corresponding maximal ORF length, and lists of ORF lengths and start sites
    ID_list.append(seq_record.id) 
    max_length_list.append(max_length_for_ID)
    length_list.append(ORF_lengths)
    start_site_list.append(ORF_start_sites)

##generate a DF summarizing all the lists, to enable easy visualization of the data:
df1=pd.DataFrame({'ID': ID_list, 'max_length_ORF': max_length_list, 'Length list': length_list, 'start sites': start_site_list})
# print df1

max_length_ORF_in_all=max(max_length_list) #calculate who is the longest ORF among all IDs in the file

##the following loop identifies which ID contains the longest ORF in the file
for i in range(len(ID_list)):
    if max_length_list[i]==max_length_ORF_in_all:
        ID_with_max=ID_list[i]
print('the maximal length of any ORF in frame %d in the file is %s bps and it is contained within the sequence %s' % (frame_num, max_length_ORF_in_all, ID_with_max))


##the following paragraph identifies the longest ORF for a given ID, and its start site:
requested_ID='gi|142022655|gb|EQ086233.1|16' ##insert here the ID of the specific sequence to investigate: 
for i in range(len(ID_list)):
    if ID_list[i]==requested_ID:
        longest_ORF_for_id=max_length_list[i]
        for j in range(len(length_list[i])):
            if length_list[i][j]==max_length_list[i]:
                start_site_for_longest_in_id=start_site_list[i][j]
                
print('in sequence %s, the longest ORF start at position %s and its length is %s bps' %(requested_ID, start_site_for_longest_in_id, longest_ORF_for_id))

'''
Q4:
Given a length n, your program should be able to identify all repeats of length n in all sequences in the 
FASTA file. Your program should also determine how many times each repeat occurs in the file, and which is 
the most frequent repeat of a given length.
'''
f1=open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Coursera/PythonIntro/dna2.fasta','r')
n=7
specific_seq='AATGGCA' ##insert here specific repeat to check
seq_list=[]
for seq_record in SeqIO.parse(f1, "fasta"):
    seq_list.append(seq_record.seq)
seq_list=[str(i) for i in seq_list]
 
substring_list_in_all=[]
substring_index_in_all=[]
n_substrings_all=[]
for seq in seq_list:
    substring_list_in_record=[]
    substring_index_list=[]
    for i in range(len(seq)):
        substring=seq[i:i+n]
        if len(substring)>=n:
            substring_list_in_record.append(substring)
            substring_index_list.append(i+1)
    n_substrings=len(substring_list_in_record)
    for j in range(len(substring_list_in_record)):
        substring_list_in_all.append(substring_list_in_record[j])
        substring_index_in_all.append(substring_index_list[j])
    n_substrings_all.append(n_substrings)
substring_set_in_all=set(substring_list_in_all)
print n_substrings_all
print len(n_substrings_all)
print len(substring_list_in_all)
print len(substring_index_in_all)
print len(substring_set_in_all)
repeat_list=[]
for uniqe_substring in substring_set_in_all:
    repeats=substring_list_in_all.count(uniqe_substring)
    repeat_list.append(repeats)
    if uniqe_substring==specific_seq.upper():
        repeat_seq=substring_list_in_all.count(uniqe_substring)
        print('the sequence %s repeats %d times' %(specific_seq, repeat_seq))
max_repeats=max(repeat_list)
most_frequent_list=[]
for i in range(len(substring_set_in_all)):
    if repeat_list[i]==max_repeats:
        most_frequent_list.append(list(substring_set_in_all)[i])
      
print('the most frequent repeating sequences in length %d are %s and they repeat %d times each in the file' %(n, most_frequent_list, max_repeats))

    
