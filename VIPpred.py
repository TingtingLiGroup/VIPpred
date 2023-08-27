import numpy as np
import pandas as pd
import os, re
import scipy
import csv
from collections import Counter
import pickle
from sklearn.metrics.cluster import *
from sklearn.metrics import hamming_loss
import xgboost as xgb
from xgboost import plot_importance
import subprocess
import sys



# Gets the param passed in
File = sys.argv[1]

print(File, ': read in')

FileName = File.split('.')[0]   #Example.txt

### Data input
df=pd.read_csv(File,sep='\t')
df['AAb']=df['AAchange'].apply(lambda x:x[0])
df['AAa']=df['AAchange'].apply(lambda x:x[-1])


### Import cache data
for pkl in [x for x in os.listdir('cachedata/') if x.split('.')[-1]=='pkl']:
    dictN=pkl.split('.')[0]
    with open('cachedata/%s'%pkl,'rb') as tf:
        globals()[dictN]=pickle.load(tf)


### Obtain functional score of nsSNVs based on dbNSFP
with open('dbNSFP/dfIn/cache.vepIn','w') as w:
    for i in df[['chr','chrl','NNchange','AAchange']].drop_duplicates().index:
        line=df[['chr','chrl','NNchange','AAchange']].loc[i]
        w.write(' '.join(line.astype(str)).replace('/',' ')+'\n')


result = subprocess.run(
    'java -Xmx8g search_dbNSFP43a -i dfIn/cache.vepIn -o dfIn/cache.vepOut -w 1-9,12-18,38-43,100-101,141-142,155-157,162-163',
    shell=True,
    capture_output=True,
    text=True, 
    cwd='dbNSFP')

if result.returncode == 0 and os.path.exists('dbNSFP/dfIn/cache.vepOut'):
    print(result.stdout)
    print('dbNSFP4.3 : finished')
else:
    print('dbNSFP4.3 : error')


## feature list
featureList = [
    'feature_paxdb_log10',
    'feature_length',
    'feature_distance_abs',
    'ASAquick_asa',
    'PSIPRED_ss3H', 'PSIPRED_ss3E', 'PSIPRED_ss3C',
    'feature_LCD', 'feature_IDR',
    'dist_3d',
    'feature_netpho_max_all', 'feature_netpho_max_kin',
    'feature_psitesnv_nbt21', 
    'feature_psitesnv_ptm21',
    'AAPolarity', 'AAVolume', 'AAHydrophobicity', 'AAGrantham',
    'site_coevolve',
    'feature_PWM_varpep', 'feature_PWM_refpep', 
    'phyloP17way_primate',
    'BayesDel_addAF_score',
    'integrated_fitCons_score',
    'GERP++_RS',
    'feature_sift_psite',
    'feature_pubmed_MS_LIT', 'feature_pubmed_LT_LIT',
    'feature_psite_ptmage0', 'feature_psite_ptmage3',
    'feature_psite_existed','feature_psite_ubi21'
]

import features_evolu
import features_public
import features_regu
import features_stru


def get_Feature(df):
    ### 'feature_PWM_refpep', 'feature_PWM_varpep', 'feature_PWM_maxdis',
    df=features_regu.get_PWM(df)
    
    ### 'feature_paxdb_log10'
    df=features_public.getPaxdb(df)

    # ### 'feature_distance_abs', 'feature_length',
    df=features_stru.getF_1d(df)

    ## 'PSIPRED_ss3H', 'PSIPRED_ss3E', 'PSIPRED_ss3C', 'ASAquick_asa', 'feature_psitesnv_LCD', 'feature_psitesnv_IDR',
    df=features_stru.getF_2d(df)

    ### 'feature_dist_3d',
    df=features_stru.getF_3d(df)

    ### 'feature_psite_ptmage0', 'feature_psite_ptmage3', 'feature_pubmed_n',
    df=features_evolu.get_ptmage_pubmed(df)

    ### 'feature_netpho_max_all', 'feature_netpho_max_kin', 'feature_netpho_max_nokin',
    df=features_regu.get_netpho(df)

    ### 'feature_psitesnv_ptm21', 'feature_psitesnv_nbt21', 'feature_psitesnv_nbt_score',
    df=features_regu.get_ptm21_nbt(df)
    df=features_regu.get_Ubi21(df)

    ### aaindex_SNVbox
    df=features_stru.get_AAindex(df)
    
    ### 'site_coevolve', 
    df=features_evolu.get_coevolve(df)

    ### 'feature_sift_psite', 'feature_provean_psite',
    df=features_evolu.get_sift(df)

    ###  Determine whether the phosphorylation site is a known phosphorylation site in PhosphositePLUS or as reported by Ochoa et al.
    df=features_public.get_psite_existed(df)
    return df


def get_predLabel(df,model,col1,col2,features=featureList):
    ## gain loss都预测一遍
    cachedf=df.fillna({
        'feature_pubmed_MS_LIT': 0,
        'feature_pubmed_LT_LIT': 0,
        'feature_psite_ptmage0': 0,
        'feature_psite_ptmage3': 0,
        'feature_psitesnv_compbias': 0,
        'feature_LCD': 0,
        'feature_IDR': 0,
        'feature_psitesnv_nbt21': 0,
        'feature_psitesnv_ptm21': 0,
    })
    
    Xtest = cachedf[features].values
    probas = model.predict_proba(Xtest)
    print('%s get'%col1)
    label = model.predict(Xtest)
    print('%s get'%col2)
    cachedf[col1]=probas[:, 1]
    cachedf[col2]=label
    
    return cachedf


def getLabel(row,col1,col2,col1s,col2s):
    if row[col1]==0 and row[col2]==0:
        return 'pairNoimpact'
    elif row[col1]==0 and row[col2]==1:
        return 'pairLoss'
    elif row[col1]==1 and row[col2]==0:
        return 'pairGain'
    elif row[col1s]>row[col2s]:
        return 'pairGain'
    elif row[col1s]<row[col2s]:
        return 'pairLoss'
    else:
        print(row.index)


## calculate feature values
df=features_public.get_vep(df, 'dbNSFP/dfIn/cache.vepOut')
df=get_Feature(df)

## Load VIPpred model
modelG=pickle.load(open('modelGain.pickle.dat','rb'))
modelL=pickle.load(open('modelLoss.pickle.dat','rb'))

df=get_predLabel(df,modelG,'gainScore','gainLabel')
df=get_predLabel(df,modelL,'lossScore','lossLabel')
df['VIP_Label']=df.apply(lambda row:getLabel(row,'gainLabel','lossLabel','gainScore','lossScore'),axis=1)

## Output
df.to_csv('%s.Output.txt'%(FileName),sep='\t')
