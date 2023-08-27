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

## Calculate the sequential distance between between a phosphorylation site and nsSNV (non-synonymous single nucleotide variant)
def getF_1d(df):
	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		loc=int(line.snvLoc)
		ploc=int(line.ploc)
		length=len(uniprotSeq[uni])
		df.loc[i,'feature_length']=length
		df.loc[i,'feature_distance_abs']=abs(loc-ploc)
	return df



## Calculate the protein secondary structure, solvent accessibility, LCD and IDR
def getF_2d(df):
	for i in df.index:
		try:
			line=df.loc[i]
			uni=line.uni
			loc=int(line.snvLoc)-1
			ploc=int(line.ploc)-1

			PSIPREDseq=PSIPREDsecondaryStructure[uni]
			ASAquickseq=ASAquicksolventAccessibility[uni]
			LCDseq=LCD[uni]
			IDRseq=IDRdescrib[uni]

			df.loc[i,'ASAquick_asa']=int(ASAquickseq[ploc])+int(ASAquickseq[loc])
			df.loc[i,'feature_LCD']=int(LCDseq[ploc])+int(LCDseq[loc])
			df.loc[i,'feature_IDR']=int(IDRseq[ploc])+int(IDRseq[loc])
			
			H,E,C=0,0,0
			ss3p,ss3s=PSIPREDseq[ploc],PSIPREDseq[loc]
			df.loc[i,'PSIPRED_ss3p%s'%ss3p]=1
			df.loc[i,'PSIPRED_ss3s%s'%ss3s]=1
		except Exception as e:
#			 print(e)
			continue

	for ss3 in list('HEC'):
		if 'PSIPRED_ss3p%s'%ss3 not in list(df.columns):
			df['PSIPRED_ss3p%s'%ss3]=np.nan
		if 'PSIPRED_ss3s%s'%ss3 not in list(df.columns):
			df['PSIPRED_ss3s%s'%ss3]=np.nan
		df['PSIPRED_ss3%s'%ss3]=df['PSIPRED_ss3p%s'%ss3].fillna(0)+df['PSIPRED_ss3s%s'%ss3].fillna(0)
	return df



## Calculate the three-dimensional spatial distance between a phosphorylation site and an nsSNV based on the results of AlphaFold2
def dist_3d(pro_id, site1, site2, res1, res2):
	RV = np.zeros(site1.shape[0])
	RV[:] = None
	# annotation of residue name
	res_name1 = ['G','A','V','L','I','P','F','Y','W','S','T','C','M','N','Q','D','E','K','R','H']
	res_name2 = [x.upper() for x in ['Gly','Ala','Val','Leu','Ile','Pro','Phe','Tyr','Trp','Ser','Thr','Cys','Met','Asn','Gln','Asp','Glu','Lys','Arg','His']]
	data_dir = "cachedata/protein3d/"

	# check the protein data and load it
	if os.path.exists(data_dir + pro_id + '.txt') == False:
		print("No 3D structure data for %s in PDB!" %pro_id)
		return RV
	data3d = np.loadtxt(data_dir + pro_id + '.txt', delimiter='\t', skiprows=1, dtype="str")

	for i in range(site1.shape[0]):
		# check the transfer the residue name
		if res_name1.count(res1[i]) == 1:
			res1_tmp = res_name2[res_name1.index(res1[i])]
		else:
			print("No such residue of %s." %res1[i])
			continue
		if res_name1.count(res2[i]) == 1:
			res2_tmp = res_name2[res_name1.index(res2[i])]
		else:
			print("No such residue of %s." %res2[i])
			continue

		# find the same residue and the types
		idx1 = (data3d[:,2] == str(site1[i])) * (data3d[:,3] == res1_tmp)
		idx2 = (data3d[:,2] == str(site2[i])) * (data3d[:,3] == res2_tmp)
		data1, data2 = data3d[idx1,:], data3d[idx2,:]

		types1, types2 = [], []
		for j in range(data1.shape[0]):
			types1.append("_".join([data1[j,0],data1[j,1],data1[j,4]]))
		for j in range(data2.shape[0]):
			types2.append("_".join([data2[j,0],data2[j,1],data2[j,4]]))

		# calculate the average distance
		RV_tmp, cnt = 0, 0
		for j in range(len(types1)):
			if types2.count(types1[j]) == 1:
				_idx = types2.index(types1[j])
				dot1 = np.array(data1[j,5:8], "float")
				dot2 = np.array(data2[_idx,5:8], "float")
				dist_tmp = np.sqrt(np.sum((dot1-dot2)**2))
				RV_tmp += dist_tmp
				cnt += 1
		if cnt > 0:
			RV[i] = RV_tmp / (cnt+0.0)
	return RV


def getF_3d(df):
	dis_3dD=[]
	for uni in df.uni.drop_duplicates():
		cache=df[(df.uni==uni)]
		
		site1=np.array(cache.snvLoc)
		res1=np.array(cache.AAb)

		site2=np.array(cache.ploc)
		res2=np.array(cache.psite)
		
		dis_3d = dist_3d(uni, site1, site2, res1, res2)
		
		t=0
		for i in cache.index:
			dis_3dD.append([i,dis_3d[t]])
			t+=1
		
	dis_3dD=pd.DataFrame(dis_3dD,columns=['index','dist_3d'])
	dis_3dD.index=list(dis_3dD['index'])
	df=pd.merge(df,dis_3dD,left_index=True,right_index=True,how='left')
		
	return df.drop(columns=['index'])




## Calculate amino acid property score from snvbox database
def get_AAindex(df):

	"""The Grantham distance attempts to provide a proxy for the evolutionary distance between two amino acids
	 based on three key side chain chemical properties: composition, polarity and molecular volume.
	 In turn, evolutionary distance is used as a proxy for the impact of missense substitutions. 
	 The higher the distance, the more deleterious the substitution is expected to be."""

	aaindexdf=pd.read_csv('cachedata/AAindex_snvbox.txt',sep='\t',index_col=0)
	for i in df.index:
		line=df.loc[i]
		AAchange=line.AAchange
		cache=aaindexdf.loc[AAchange]
		df.loc[[i],['AAVolume','AAHydrophobicity','AAGrantham','AAPolarity','AAEx','AAPAM250','AAJM','AAVB']]=list(cache[['Volume','Hydrophobicity','Grantham','Polarity','Ex','PAM250','JM','VB']])
	return df
