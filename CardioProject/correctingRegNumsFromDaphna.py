import pandas
from _mysql_exceptions import Error
from Utils import Load
resultsExcel='/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/PNPChip/tmpnonImputedChipFilesMAF_Single_StoolStoolVsChipAnalysis.xlsx'
matchingResultsExcel=pandas.read_excel(resultsExcel,sheetname='All')
matchingResultsExcel.index=matchingResultsExcel.index.astype('str')
matchingResultsExcel.correctStool_id=matchingResultsExcel.correctStool_id#.astype('str')#.apply(lambda x: str(x).split('.')[0])
matchingResultsExcel.batch_winners=matchingResultsExcel.batch_winners.apply(eval)
resolved=['1','3!','6!!']
multipleStoolForRegnum=Load('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/PNPChip/MultipleStoolForRegTableAllSamples.dat')
locationAndDate=Load('/net/mraid08/export/genie/Microbiome/Analyses/Daphna/PNPChip/ParticipantsDateLocationInFridgeWithDates.dat').set_index('regNum')
batch3regNums=pandas.read_csv('/net/mraid08/export/genie/Microbiome/Analyses/PNPChip/AdditionalExperiments/Batch3RegNumBDConversion.csv').set_index('RegNum')
import numpy as np
def resolvedType(conf):
    return conf[0] in resolved or\
    conf[:2] in resolved or\
    conf[:3] in resolved
    
def fillCorrectStoolID(df):
    df = df.reset_index()
    df.loc[df['correctStool_id'].notnull(),'correctStool_id']=df.loc[df['correctStool_id'].notnull(),'correctStool_id'].\
                                                                apply(lambda x: str(x).split('.')[0]).astype('str')
    df.loc[df['correctStool_id'].isnull(),'correctStool_id'] = df.loc[df['correctStool_id'].isnull(),'index'].\
                                                                apply(lambda x: str(x).split('@')[0]).astype('str')
    df=df.set_index("correctStool_id")
    return df

def participantHasNoStool_and_HisBlood_wasnt_assigned_to_someone_else(regNum,df):
    winners=sum(df[df['resolved']].batch_winners.tolist(),[])
    if regNum not in df.index and regNum not in winners:
        return True
    else:
        return False
def participantFromNewBatch(regNum,df):
    return (int(regNum) in batch3regNums.index)
def returnType(func):
    def wrapFunc(regNum):
        returnType=type(regNum)
        if returnType!=str:
            returnType=int
        regNum=str(regNum).split(".")[0]
        res=func(regNum)
        typeres=returnType(str(res).split(".")[0]) if res is not None and res!='nan' else None
#         if type(typeres)!=int:
#             raise Error
        return typeres
    return wrapFunc

@returnType
def getQCedBloodRegNumByRegNum(regNum):
    # tested in test_QCgetters 
    df = matchingResultsExcel
    df =fillCorrectStoolID(df)
    df['resolved']=df.confidenceBatch.apply(resolvedType)
    if participantHasNoStool_and_HisBlood_wasnt_assigned_to_someone_else(regNum,df):
        return regNum
    winners=df[df['resolved']].batch_winners
    if regNum in winners.index:
        match=winners.loc[regNum]
        if type(match)==list:
            summed_match=match
        else:
            summed_match = list(set(sum(match.values,[]))) 
        if len(summed_match)==1:
            return summed_match[0]
        else:
            raise Exception("Bad resolving, conflicting blood label for DB_ID %s"%regNum)
    if participantFromNewBatch(regNum,df):
        return regNum
    return None

# def getLocationInFridge(regnum):
#     correctRegNum=getQCedBloodRegNumByRegNum(regnum)
#     return locationAndDate.loc[correctRegNum]
# @returnType
# def getQCedStoolRegNumByRegNum(regNum):
#     df = matchingResultsExcel
#     df =fillCorrectStoolID(df)
#     if regNum in df.index.values:
#         if type(df.loc[regNum]['index'])==str:
#             return df.loc[regNum]['index']
#         else:
#             res= set([val.split('@')[0] for val in list(set(df.loc[regNum]['index']))])
# #             try:
# #                 assert len(res)==1
# #             except:
# #                 print regNum
#             res = list(res)[0]
#             return res
#     if regNum in df["index"].values:
#         return None
#     return regNum
# 
# def getFDbySample(sampleID):
#     return 'FD'+multipleStoolForRegnum[multipleStoolForRegnum['outputRegNum']==sampleID.replace('@','_')]\
#             .index.values[0].split('FD')[1].split('.')[0]
#             
# def getAllFDs():
#     return ['FD'+val.split('FD')[1].split('.')[0] for val in multipleStoolForRegnum.index.values if 'FD' in val]
# 
# def getQCedAllStoolRegNumByRegNum(regNum,ExtraFDs=None):
#     if type(regNum)==int:
#         regNum=str(regNum)
#     df = matchingResultsExcel
#     df =fillCorrectStoolID(df).reset_index().set_index('index')
#     
#     df_reindexed = df.reset_index()
#     stool_ids = list(df_reindexed[df_reindexed['correctStool_id']==regNum]['index'].values)
#     stool_ids = list(set([x + '@0' if '@' not in x else x for x in stool_ids]))
#     stool_fds = [getFDbySample(x) for x in stool_ids]
#     if ExtraFDs is not None:
#         newSamples=list(set(ExtraFDs)-set(getAllFDs()))
#         stool_fds+=newSamples
# #     print zip(stool_ids,stool_fds)
#     return stool_fds
# 
# 
# def getFDForBadStool():
#     signatureDF = pandas.read_excel('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/PNPChip/EladDaphaEran_Approved_and_Overriden.xlsx')
#     badStoolIds=signatureDF[signatureDF["correctStool_id"]=='X']['stool_id'].values
#     multipleStoolForRegnum=Load('/net/mraid08/export/jafar/Microbiome/Analyses/Daphna/PNPChip/MultipleStoolForRegTableAllSamples.dat')
#     fdnames=[]
#     for badStool in badStoolIds:
#         if '@' in str(badStool):
#             regnum=badStool.split('@')[0]
#             fdname='FD'+multipleStoolForRegnum[multipleStoolForRegnum['outputRegNum']==badStool.replace('@','_')]\
#             .index.values[0].split('FD')[1].split('.')[0]
#         else:
#             fdname='FD'+multipleStoolForRegnum[multipleStoolForRegnum['RegNum']==str(badStool)]\
#             .index.values[0].split('FD')[1].split('.')[0]
#         fdnames.append(fdname)
#     return fdnames
# 
# if __name__=="__main__":
#     getFDForBadStool()
#     regNums=pandas.read_csv('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/trial.csv')
#     regNums['registrationCode']=regNums['registrationCode'].apply(str) 
#     regNums['correctStoolRegistrationCode']=regNums.registrationCode.apply(getQCedStoolRegNumByRegNum)
#     regNums['correctBloodRegistrationCode']=regNums.registrationCode.apply(getQCedBloodRegNumByRegNum)
#     print getQCedBloodRegNumByRegNum('3')
    
# print getQCedBloodRegNumByRegNum(3)
    