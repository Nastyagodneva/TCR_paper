import matplotlib
matplotlib.use('Agg')

import sys
import numpy as np
from ShaniBA.PredictionPipeline import PredictionModuleShani
from SegalQueue.qp import qp, fakeqp
import os
from addloglevels import sethandlers
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists
import pandas as pd
from ShaniBA.myplots import *
from ShaniBA.MyFunctionsShani import *
import argparse
import json
import matplotlib.pyplot as plt

def run_predictions(Yname,pred_dir,model,n_permut_X,n_permut_Y, parameter_search_type,pathToModelParams,eval_metric1,eval_metric2):
    print Yname
    df=pd.DataFrame(Y_binary[Yname])
    f1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/%s.xlsx' %Yname
    f1=f1.replace(' ','')
    df.to_excel(f1)
    
    os.system('python /home/sbenari/workspace/Microbiome/ShaniBA/PredictionPipeline/permutation_predictions_for_cardio_phenotypes.py\
 /net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Cardio126_diseasePhenotypes/%(Yname)s/ -Yname %(Yname)s\
 -pred_dir %(pred_dir)s -model %(model)s -n_permut_X %(n_permut_X)s -n_permut_Y %(n_permut_Y)s -parameter_search_type %(parameter_search_type)s\
 -useShortVersion True -pathToModelParams %(pathToModelParams)s -bsY True -eval_metric1 %(eval_metric1)s -eval_metric2 %(eval_metric2)s'\
   %{'Yname':Yname,'pred_dir':pred_dir,'model':model,'n_permut_X':n_permut_X,'n_permut_Y':n_permut_Y,'parameter_search_type':parameter_search_type,
     'pathToModelParams':pathToModelParams,'eval_metric1':eval_metric1,'eval_metric2':eval_metric2})
    
    return 

pred_dir='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Cardio126_diseasePhenotypes/logRegStrat_gridSearch_byRepFeatPCA10percVDJ0999PredictedAgeGender_optByKappa/'
nPreSelectedFeatures=9999
n_permut_X=19
n_permut_Y=1
model= 'lin'
parameter_search_type='gridSearch'
pathToModelParams ='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Model_params_files/logReg_model_params_1.xlsx'
eval_metric1='kappa'
eval_metric2='Precision_Recall'

Y_file_binary='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/Cardio126phenAllInfo_Nov18_binary.xlsx'
Y_binary=pd.read_excel(Y_file_binary).set_index('BD')



tempDir = "/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir"
if not isdir(tempDir):
    makedirs(tempDir)
os.chdir(tempDir)
print ('setting handlers...')
sethandlers()
print 'start'
with qp(jobname='predWrap', q=['himem7.q'], mem_def="1G", trds_def=1, deleteCSHwithnoerr=True,
         tryrerun=False, max_u=350) as q:
    q.startpermanentrun()
    wait_on = []

for n,Yname in enumerate(Y_binary.columns):
    print ('now start working on %s (phenotype num %s)' %(Yname,n))
    wait_on.append(q.method(run_predictions, kwargs={'Yname':Yname,'pred_dir':pred_dir,'model':model,'n_permut_X':n_permut_X,'n_permut_Y':n_permut_Y,'parameter_search_type':parameter_search_type,
     'pathToModelParams':pathToModelParams,'eval_metric1':eval_metric1,'eval_metric2':eval_metric2}))                      
    try:
        res = q.waitforresults(wait_on)
    except:
        print 'some error occured'




