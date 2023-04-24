import numpy as np
import pandas as pd
import os, re
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import csv
from collections import Counter
import pickle
from sklearn.metrics.cluster import *
from sklearn.metrics import hamming_loss



for pkl in [x for x in os.listdir('cachedata/') if x.split('.')[-1]=='pkl']:
	dictN=pkl.split('.')[0]
	with open('cachedata/%s'%pkl,'rb') as tf:
		globals()[dictN]=pickle.load(tf)


## Obtain human consensus protein abundance deposited in PaxDb v4.2
def getPaxdb(df):
	"""
	Get protein expression from PAXdb v4.2
	"""
	# with open('cachedata/paxdbValue.pkl','rb') as tf:
	# 	paxdbValue=pickle.load(tf)

	for i in df.index:
		uni=df.loc[i,'uni']
		if paxdbValue.get(uni):
			df.loc[i,'feature_paxdb_log10']=paxdbValue[uni]
	return df



## Determine whether the phosphorylation site is a known phosphorylation site in PhosphositePLUS or as reported by Ochoa et al.
def get_psite_existed(df):
	nbt=pd.read_csv('cachedata/nbt.txt',sep='\t')
	ptmL=['Phosphorylation']
	ptmPic=pd.DataFrame()
	
	for mod in ptmL:
		cache=pd.read_csv('cachedata/psp/%s_site_dataset'%mod,sep='\t',skiprows=3,low_memory=False)
		cache=cache[cache.ORGANISM=='human']
		cache['Mod']=mod
		ptmPic=ptmPic.append(cache)
	ptmPic=ptmPic.reset_index(drop=True)
	ptmPic['siteLoc']=ptmPic['MOD_RSD'].str.split('-',expand=True)[0].apply(lambda x:int(x[1:]))
	
	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		loc=int(line.snvLoc)
		psite=line.psite
		ploc=int(line.ploc)
		
		tC=0
		if nbt[(nbt.uniprot==uni)&(nbt.position==ploc)].shape[0]>0:
			tC+=1
		if ptmPic[(ptmPic.ACC_ID==uni)&(ptmPic.siteLoc==ploc)].shape[0]>0:
			tC+=1
		df.loc[i,'feature_psite_existed']=tC
		
	return df
