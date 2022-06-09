import sys
import numpy as np
from ShaniBA.PredictionPipeline import PredictionModuleShani
from queue.qp import qp, fakeqp
import os
from addloglevels import sethandlers
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists


# from joblib import Parallel, delayed
# import multiprocessing
# 
# 
# 
# from multiprocessing import Pool

def main(i):
        
    print i
    
    if i<10:
        Xpath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/X_withPredictedAgeGender_Cardio126_shuffling/shuf%s.dat' % i

        os.system('python /home/sbenari/workspace/Microbiome/ShaniBA/PredictionPipeline/PredictionModuleShani.py\
         /net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Cardio126_diseasePhenotypes/GRACEScore/\
shuff%s_XGBreg20LOO_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender_expVar/\
         -path_to_X %s\
         -path_to_Y /net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/GRACEScore.xlsx\
         -parameter_search_type randomSearch -n_random 25 -model XGB -useShortVersion True\
          -mem_def 1 -n_cols_per_job 1 -nPreSelectedFeatures 20 -k_folds 0\
          -predictor_params \'{"learning_rate":[0.1,0.05,0.01,0.005],"n_estimators":\
[200,600,1000,1400],"reg_alpha":[0,1,5,10],"gamma":[0,1,5,10],"subsample":[0.5,0.7],"max_depth":[1,3,5,7,9]}\'' % (i, Xpath))
    
    
    else:
        k=i-10
        print ('k= ', k)
        
        Xpath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/X_withPredictedAgeGender_Cardio126.dat' 
        Ypath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/GRACEScore_bsFiles/bs_%s.xlsx' %k
        
        
        os.system('python /home/sbenari/workspace/Microbiome/ShaniBA/PredictionPipeline/PredictionModuleShani.py\
         /net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Cardio126_diseasePhenotypes/GRACEScore/\
bs%s_XGBreg20LOO_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender_expVar/\
         -path_to_X %s\
         -path_to_Y %s\
         -parameter_search_type randomSearch -n_random 25 -model XGB -useShortVersion True\
          -mem_def 1 -n_cols_per_job 1 -nPreSelectedFeatures 20 -k_folds 0\
          -predictor_params \'{"learning_rate":[0.1,0.05,0.01,0.005],"n_estimators":\
[200,600,1000,1400],"reg_alpha":[0,1,5,10],"gamma":[0,1,5,10],"subsample":[0.5,0.7],"max_depth":[1,3,5,7,9]}\'' % (k, Xpath,Ypath))
    
        
        
    return





sethandlers()
outputDir = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Cardio126_diseasePhenotypes/GRACEScore/temps'
if not isdir(outputDir):
    makedirs(outputDir)
os.chdir(outputDir)
print 'start'
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp(jobname='predWrap', q=['himem7.q'], mem_def="1G", trds_def=4, deleteCSHwithnoerr=True,
         tryrerun=False, max_u=120) as q:
    q.startpermanentrun()
    wait_on = []
    
for i in range(20):
    print i
    wait_on.append(q.method(main, kwargs={'i':i}))
                       
q.wait(wait_on)
    
