'''
this script is based on the jupyter notebook 'SampleFilesManipulations' which was used to organize lists
of blood DNA samples
'''


from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import time

cdate = str(time.strftime("%d%m%Y"))
print ('current date is %s' % cdate)

#--------------------------------------------------------------------------------
# # (1) defining file names:

# # a. file with all regNums that should be corrected using Daphna's function:
regNumListFile = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/File1.xlsx'
BloodSampleFile = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/File2.xlsx'
BloodDNASamplesFile = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/File3.xlsx'
UserAndRegNumbersFile = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/File4.xlsx'
SignedHelsinkiFile = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/File5.xlsx'
SampleNamesAndConcsTali = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/File6.xlsx'


#--------------------------------------------------------------------------------

# # (2)  finding the correct registration codes using Daphna's functions
 
from PNP_QC.getQCedMatchRegNum import getQCedBloodRegNumByRegNum
#from PNP_QC import getQCedBloodRegNumByRegNum



file_name = regNumListFile
regNumListFile = pd.read_excel(file_name)

print ('regNumListFile: ')
regNumListFile.head()

reg_list = list(regNumListFile['RegistrationCode'])

# # a loop to find the correct reg code for each original reg code:

count = 0
correct_reg_list = []
CorrectionStatus_list = []
only_corrected = []
for reg in reg_list:
    Correct = getQCedBloodRegNumByRegNum(reg)
    
    if Correct == reg:
        CorrectionStatus_list.append('same')
        correct_reg_list.append(Correct)
    else:
        count = count + 1
        if Correct == None:
            CorrectionStatus_list.append('correct reg unknown')
            correct_reg_list.append('unknown')
        else:
            CorrectionStatus_list.append('corrected')
            correct_reg_list.append(Correct)
            only_corrected.append(Correct)
        
        
        
CorrectRegDF = pd.DataFrame({'original registration code': reg_list,
                           'correct registration code':correct_reg_list,
                           'correction status (reg number)': CorrectionStatus_list })

print('CorrectRegDF: ')
CorrectRegDF.head()

grouping_count = CorrectRegDF.groupby('correction status (reg number)')['correct registration code'].count()
groupingCountsDF = pd.DataFrame({'counts':grouping_count})
groupingCountsDF['percentage'] = (groupingCountsDF['counts']) * 100 / groupingCountsDF['counts'].sum()
groupingCountsDF

# # saving the correct Reg table to pickles and excel:

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/CorrectRegDF_%s' %cdate, "wb") as f:
    pickle.dump(CorrectRegDF, f)
f.close()


writer = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/CorrectRegDF_%s.xlsx' % cdate
CorrectRegDF.to_excel(writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf', verbose=True)

#-------------------------------------------------------------------

# # generating comprehensive list of DNA samples extracted from blood

# # downloading the list of blood samples extracted from the DB using the following SQL script:

# #SELECT pnp.users.UserID, pnp.users.RegistrationCode, Lab.blood.UserID, 
# #Lab.blood.bid1, Lab.blood.bid2 

# #FROM pnp.users
# #LEFT JOIN Lab.blood ON pnp.users.UserID = Lab.blood.UserID

file_name = BloodSampleFile
BloodSamples_DBExtracted = pd.read_excel(file_name)
print ('BloodSamples_DBExtracted: ')
BloodSamples_DBExtracted.head()

# # removing samples with neither bid1 nor bid2, and dropping the UserID columns:

CleanBloodSamples = BloodSamples_DBExtracted[(BloodSamples_DBExtracted['bid1'].notnull()) | 
                                                  (BloodSamples_DBExtracted['bid2'].notnull())]
CleanBloodSamples = CleanBloodSamples.drop(['UserID', 'UserID.1'], axis=1)

print ('CleanBloodSamples: ')
CleanBloodSamples.head()

#-----------------------------------------------------------

# #downloading the list of blood DNA samples extracted from the DB using the following SQL script:

# #SELECT pnp.users.UserID, pnp.users.RegistrationCode, Lab.hostdna.UserID,
# #Lab.hostdna.DnaID

# #FROM pnp.users
# # LEFT JOIN Lab.hostdna ON pnp.users.UserID = Lab.hostdna.UserID

file_name = BloodDNASamplesFile
BloodDNASamples_DBExtracted = pd.read_excel(file_name)

print('BloodDNASamples_DBExtracted: ')
BloodDNASamples_DBExtracted.head()

# drop all rows that don't have BD num
BloodDNASamplesNoNan = BloodDNASamples_DBExtracted[~BloodDNASamples_DBExtracted['DnaID'].isnull()]
BloodDNASamplesNoNan

# merging the blood DNA table with the 
# correct reg numbers in order to correct this table:

CorrectBloodDNASamples = pd.merge(BloodDNASamplesNoNan, CorrectRegDF, how='left',
                                    left_on='RegistrationCode', right_on='original registration code')
print('CorrectBloodDNASamples: ')
CorrectBloodDNASamples.head()

# # make clean file with no wrong reg numbers and user IDs:
CleanCorrectBloodDNASamples = CorrectBloodDNASamples.drop(['RegistrationCode', 'original registration code',
                                                        'UserID', 'UserID.1'], axis=1)
print('CleanCorrectBloodDNASamples: ')
CleanCorrectBloodDNASamples.head()

CleanCorrectBloodDNASamples.groupby('correction status (reg number)')['DnaID'].count()

# # saving the CleanCorrectBloodDNASamples table to pickles and excel:

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/CleanCorrectBloodDNASamples_%s' % cdate, "wb") as f:
    pickle.dump(CleanCorrectBloodDNASamples, f)
f.close()

writer = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/CleanCorrectBloodDNASamples_%s.xlsx' % cdate
CleanCorrectBloodDNASamples.to_excel(writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf', verbose=True)


#-------------------------------------------------------------------
# # merge the clean blood DNA sample with the blood sample based on correct REG number:
BloodDNASamplesCorrectREG = pd.merge(CleanCorrectBloodDNASamples, CleanBloodSamples, how='left',
                                        left_on='correct registration code', right_on='RegistrationCode')
BloodDNASamplesCorrectREG = BloodDNASamplesCorrectREG.drop('RegistrationCode', axis=1)

print('BloodDNASamplesCorrectREG: ')
BloodDNASamplesCorrectREG.head()

#-------------------------------------------------------------------

# # downloading file with usersID and regCode:
file_name = UserAndRegNumbersFile
UserIDs_and_regNums = pd.read_excel(file_name)

print ('UserIDs_and_regNums: ')
UserIDs_and_regNums.head()

# add user IDs to the blood DNA file:
BloodDNASamplesCorrectREGwithUserID = pd.merge(BloodDNASamplesCorrectREG, UserIDs_and_regNums, how='left',
                                        left_on='correct registration code', right_on='RegistrationCode')

print('BloodDNASamplesCorrectREGwithUserID: ')
BloodDNASamplesCorrectREGwithUserID.head()

# remove redundent reg code column:
BloodDNASamplesCorrectREGwithUserID = BloodDNASamplesCorrectREGwithUserID.drop('RegistrationCode', axis=1)
print ('BloodDNASamplesCorrectREGwithUserID: ')
BloodDNASamplesCorrectREGwithUserID.head()

#-----------------------------------------------------------------------

file_name = SignedHelsinkiFile
signedHelsinki = pd.read_excel(file_name)
print('signedHelsinki: ')
signedHelsinki.head()

BloodDNASamples = pd.merge(BloodDNASamplesCorrectREGwithUserID, signedHelsinki, how='left',
                                    left_on='UserID', right_on='UserID')
print('BloodDNASamples: ')
BloodDNASamples.head()

# remove duplicates:

BloodDNASamples_2 = BloodDNASamples.drop_duplicates()

grouped_future_bids = BloodDNASamples_2.groupby('dnaid')['bid1'].apply(list)
grouped_future_bids = grouped_future_bids.reset_index()
grouped_future_bids[50:100] 
 

print ('BloodDNASamples list length is %s' %len(BloodDNASamples))
print ('BloodDNASamples_2 list length is %s' %len(BloodDNASamples_2))
print ('rouped_future_bids length is %s' %len(grouped_future_bids))




# MAKE THE BID1 column a summary of all bid1 samples per DNA sample:


BloodDNASamples_3 = pd.merge(BloodDNASamples_2, grouped_future_bids, how='left',
                                    left_on='DnaID', right_on='DnaID')
BloodDNASamples_3 = BloodDNASamples_3.rename(columns={'bid1_y':'all bid1s'}).drop('bid1_x', axis=1)

print('BloodDNASamples_3: ')
BloodDNASamples_3.head()

# remove duplicate lines for the same DNA sample:

BloodDNASamples_4 = BloodDNASamples_3.drop_duplicates(subset='DnaID')
print ('BloodDNASamples_4 list length is %s' %len(BloodDNASamples_4))
print ('BloodDNASamples_3 list length is %s' %len(BloodDNASamples_3))

print('BloodDNASamples_4: ')
BloodDNASamples_4.head()

#----------------------------------------------------------------
# # this is the final table for samples extracted from the database

# #this table can be merged with each updated file from Tali

# #Save the final table:

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/BloodDNASamples_4_%s' %cdate, "wb") as f:
    pickle.dump(BloodDNASamples_4, f)
f.close()


writer = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/BloodDNASamples_4_%s.xlsx' %cdate
BloodDNASamples_4.to_excel(writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf', verbose=True)


#-----------------------------------------------------------

# # merge the correct blood and DNA samples with Tali's table:
# downloading samples that Tali already cleaned with their updated concentrations

file_name = SampleNamesAndConcsTali
CleanedSamplesOnlyNamesConc = pd.read_excel(file_name)

print('CleanedSamplesOnlyNamesConc: ')
CleanedSamplesOnlyNamesConc.head()

# # merge the correct blood and DNA samples with Tali's table:

BloodDNASamples_CleaningInfo = pd.merge(BloodDNASamples_4, CleanedSamplesOnlyNamesConc, how='outer',
                                    left_on='DnaID', right_on='Sample name')
print ('BloodDNASamples_CleaningInfo: ')
BloodDNASamples_CleaningInfo.head()

print ('BloodDNASamples_CleaningInfo length is %s' %len(BloodDNASamples_CleaningInfo))

BloodDNASamples_CleaningInfo['enough for 0.5ug'] = np.where(BloodDNASamples_CleaningInfo[
        'conc. after cleanup (ng/ul)'] > 32.25, 1, 0)

BloodDNASamples_CleaningInfo['need to re-clean'] = np.where(BloodDNASamples_CleaningInfo[
        'conc. after cleanup (ng/ul)'] < 32.25, 1, 0)


print('BloodDNASamples_CleaningInfo: ')
BloodDNASamples_CleaningInfo.head()

# save the table with the cleaning info into excel and pickles:


with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/BloodDNASamples_CleaningInfo', "wb") as f:
    pickle.dump(BloodDNASamples_CleaningInfo, f)
f.close()


writer = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Old Computer Files/Sample files/BloodDNASamples_CleaningInfo.xlsx'
BloodDNASamples_CleaningInfo.to_excel(writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf', verbose=True)




