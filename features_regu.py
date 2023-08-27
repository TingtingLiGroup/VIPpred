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

## Obtain phosphorylation regulation information based on the NetPho method
def get_netpho(df):
	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		psite=line.psite
		ploc=int(line.ploc)
		key=(uni,ploc,psite)
		if netPho.get(key):
			netphoL=netPho[key]
			df.loc[i,'feature_netpho_max_all']=netphoL[0]
			df.loc[i,'feature_netpho_max_kin']=netphoL[1]
			df.loc[i,'feature_netpho_max_nokin']=netphoL[2]
	return df



## Calculate the number of neighboring PTMs (+/-10 residues) included in PhosphositePLUS or as reported by Ochoa et al.
def get_ptm21_nbt(df):
	nbt=pd.read_csv('cachedata/nbt.txt',sep='\t')
	
	ptmL=['Acetylation','Methylation','Phosphorylation','Sumoylation','Ubiquitination','O-GalNAc','O-GlcNAc']
	ptmPic=pd.DataFrame()
	for mod in ptmL:
		cache=pd.read_csv('cachedata/PTM_21/%s_site_dataset'%mod,sep='\t',skiprows=3,low_memory=False)
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
		
		cache=nbt[(nbt.uniprot==uni)&(nbt.position==ploc)]['functional_score']
		if cache.shape[0]>0:
			df.loc[i,'feature_psite_nbtScore']=cache.iloc[0]
		
		
		locL21=set(range(loc-10,loc+11))|set(range(ploc-10,ploc+11))
		df.loc[i,'feature_psitesnv_nbt21']=nbt[(nbt.uniprot==uni)&(nbt.position.isin(locL21))].shape[0]
		df.loc[i,'feature_psitesnv_ptm21']=ptmPic[(ptmPic.ACC_ID==uni)&(ptmPic.siteLoc.isin(locL21))].shape[0]  
	return df



## Calculate the number of neighboring Ubiquitination (+/-10 residues) included in PhosphositePLUS or as reported by Ochoa et al.
def get_Ubi21(df):
	nbt=pd.read_csv('cachedata/nbt.txt',sep='\t')
	
	ptmL=['Ubiquitination']
	ptmPic=pd.DataFrame()
	for mod in ptmL:
		cache=pd.read_csv('cachedata/PTM_21/%s_site_dataset'%mod,sep='\t',skiprows=3,low_memory=False)
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
		locL21=set(range(loc-10,loc+11))|set(range(ploc-10,ploc+11))
		df.loc[i,'feature_psite_ubi21']=ptmPic[(ptmPic.ACC_ID==uni)&(ptmPic.siteLoc.isin(locL21))].shape[0]
	return df



## Calculate PWM similarity score
counts={'S': 942147, 'T': 602864, 'Y': 301037}

# counts={'S':0,'T':0,'Y':0}
# for uni in ASAquicksolventAccessibility:
# 	if uniprotSeq.get(uni):
# 		seq=uniprotSeq[uni]
# #	 else:
# #		 seq=ref[uni]
# 	for AA in counts:
# 		counts[AA]+=seq.count(AA)

def getseq_ref(site,siteloc,uni):

	if ref.get(uni):
		if siteloc-8>=0:
			subseqb=ref[uni][siteloc-8:siteloc-1]
		else:
			part_subseqb=ref[uni][0:siteloc-1]
			subseqb=''.join(['_']*(7-len(part_subseqb)))+part_subseqb

		subseqa=ref[uni][siteloc:siteloc+7]
		if len(subseqa)<7:
			subseqa=subseqa+''.join(['_']*(7-len(subseqa)))

		if ref[uni][siteloc-1]==site:
			if 'X' in subseqb+'_'+subseqa or 'U' in subseqb+'_'+subseqa:
				return 'error AA X or U'
			else:
				return site,subseqb+'_'+subseqa ## noSTY
		else:
			return 'psite wrong'

def getseq_var(site,siteloc,uni,snvloc,AAc):
	if siteloc!=snvloc:
		if ref.get(uni):

			seq=ref[uni]
			if seq[snvloc-1]==AAc[0]:

				seq=list(seq)
				seq[snvloc-1]=AAc[2]
				seq=''.join(seq)

				if siteloc-8>=0:
					subseqb=seq[siteloc-8:siteloc-1]
				else:
					part_subseqb=seq[0:siteloc-1]
					subseqb=''.join(['_']*(7-len(part_subseqb)))+part_subseqb

				subseqa=seq[siteloc:siteloc+7]
				if len(subseqa)<7:
					subseqa=subseqa+''.join(['_']*(7-len(subseqa)))

				if seq[siteloc-1]==site:
					
					if 'X' in subseqb+'_'+subseqa or 'U' in subseqb+'_'+subseqa:
						return 'error AA X or U'
					else:
						return site,subseqb+'_'+subseqa ## noSTY
				else:
					return 'Psite wrong'

			else:
				return 'Mutation wrong'
	else:
		return 'sameLoc mutLoc=psiteLoc'

def subseqNoSTY(df):
	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		site=line.psite
		siteloc=int(line.ploc)
	#	 AAchange=line.AAchange
	#	 site,siteloc=line['p-site'].split('_')[0][-1],int(line['p-site'].split('_')[0][:-1])
	#	 snvLoc=int(line.snvLoc)

		try:
			aa,subseq_ref=getseq_ref(site,siteloc,uni)

			df.loc[i,'subseq_ref']=subseq_ref
			df.loc[i,'subseq_ref_pA']=aa


		except Exception as e:
			print(e,site,siteloc,uni)


	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		site=line.psite
		siteloc=int(line.ploc)
		AAchange=line.AAchange
		snvLoc=int(line.snvLoc)

		try:
			aa,subseq_var=getseq_var(site,siteloc,uni,snvLoc,AAchange)
			if aa!='sameLoc':
				df.loc[i,'subseq_var']=subseq_var
				df.loc[i,'subseq_var_pA']=aa

		except Exception as e:
			print(e,site,siteloc,uni,snvLoc,AAchange)
	#		 print(line['p-site']+'_'+uni,AAchange,snvLoc,getseq_ref(site,uni,snvLoc),'---',e)
	
	return df

def MSSmatrix(df,kinL,Ksub,conser_df):
	## MSS score
	mss_ref=pd.DataFrame()
	mss_var=pd.DataFrame()

	for i in df.index:
		subseq_ref=df.loc[i,'subseq_ref']
		subseq_var=df.loc[i,'subseq_var']
		subseq_ref_pA=df.loc[i,'subseq_ref_pA']
		subseq_var_pA=df.loc[i,'subseq_var_pA']
	#	 break
		if pd.notna(subseq_ref):
			for kin in kinL:
				if subseq_ref_pA in Ksub[kin]:
					mins=globals()[kin].min()*conser_df.loc[kin]
					maxs=globals()[kin].max()*conser_df.loc[kin]

					mins=(mins).sum()-mins[0]
					maxs=(maxs).sum()-maxs[0]

					curr=[]
					locs=-7
					for aai in subseq_ref:
						if aai!='_':
							curr.append(globals()[kin].loc[aai,locs])
						else:
							curr.append(0)
						locs+=1
					currs=(curr*conser_df.loc[kin]).sum()
					mss_ref.loc[subseq_ref,kin]=(currs-mins)/(maxs-mins)


		if pd.notna(subseq_var):
			for kin in kinL:
				if subseq_var_pA in Ksub[kin]:
					mins=globals()[kin].min()*conser_df.loc[kin]
					maxs=globals()[kin].max()*conser_df.loc[kin]

					mins=(mins).sum()-mins[0]
					maxs=(maxs).sum()-maxs[0]

					curr=[]
					locs=-7
					for aai in subseq_var:
						if aai!='_':
							curr.append(globals()[kin].loc[aai,locs])
						else:
							curr.append(0)
						locs+=1
					currs=(curr*conser_df.loc[kin]).sum()
					mss_var.loc[subseq_var,kin]=(currs-mins)/(maxs-mins)
	
	return mss_ref,mss_var

def get_PWM(df):
	kin_sub=pd.read_csv('cachedata/Kinase_Substrate_Dataset',sep='\t',skiprows=3)
	kin_sub=kin_sub[(kin_sub.SUB_ORGANISM=='human')&(kin_sub.KIN_ORGANISM=='human')]
	kinase=kin_sub.KINASE.value_counts()
	kinL=kinase[kinase>=10].index
	
	# Obtain amino acid frequency distribution of whole human proteome
	aaf=pd.read_csv('cachedata/aa_count_percentage.csv',index_col=0)

	# background division
	Ksub = {}
	for kin in kinL:
		cache = pd.DataFrame(index=aaf.index)
		subdf = kin_sub[kin_sub.KINASE == kin]['SITE_+/-7_AA']
		subdf = pd.DataFrame([list(x.upper()) for x in subdf])
		subdf.columns = [x - 7 for x in subdf.columns]
		Ksub[kin] = set(subdf[0])

		for i in subdf.columns:
			# frequency
			cachee = pd.DataFrame(aaf.percentage)
			loc_counts = pd.merge(cachee,
								  subdf[i].value_counts(),
								  left_index=True,
								  right_index=True,
								  how='left')
			loc_counts[i] = loc_counts[i].fillna(0) + 1
			loc_counts[i] = (loc_counts[i] / (loc_counts[i].sum()))

			cache[i] = loc_counts[i] / loc_counts['percentage']
	#		 cache[i]=cache[i].apply(np.log2)

		globals()[kin] = cache

	# conservation
	conser_df=pd.DataFrame(columns=list(range(-7,8,1)))
	for kin in kinL:
		for i in conser_df.columns:
			conser_df.loc[kin,i]=sum([x*np.log(x) for x in globals()[kin][i]])
			
	global mss_ref
	global mss_var
	df=subseqNoSTY(df)
	mss_ref,mss_var=MSSmatrix(df,kinL,Ksub,conser_df)
	
	df = pd.merge(df, pd.DataFrame(mss_ref.max(axis=1), columns=['feature_PWM_refpep']),
					   left_on='subseq_ref',
					   right_index=True,
					   how='left')
#	 df = pd.merge(df, pd.DataFrame(mss_ref.idxmax(axis=1), columns=['PWM_refpep_name']),
#						left_on='subseq_ref',
#						right_index=True,
#						how='left')

	df = pd.merge(df, pd.DataFrame(mss_var.max(axis=1), columns=['feature_PWM_varpep']),
					   left_on='subseq_var',
					   right_index=True,
					   how='left')
#	 df = pd.merge(df, pd.DataFrame(mss_var.idxmax(axis=1), columns=['PWM_varpep_name']),
#						left_on='subseq_var',
#						right_index=True,
#						how='left')
	
	return df.drop(columns=['subseq_var_pA','subseq_ref_pA'])
