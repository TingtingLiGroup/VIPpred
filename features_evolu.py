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


## Obtain experimental evidence support for phosphorylation sites based on PhosphositePLUS
## Obtain ptmage score
def get_ptmage_pubmed(df):
	Phosphorylation=pd.read_csv('cachedata/psp/Phosphorylation_site_dataset',sep='\t',skiprows=3)
	Phosphorylation['siteLoc']=Phosphorylation['MOD_RSD'].str.split('-',expand=True)[0].apply(lambda x:int(x[1:]))
	age=pd.read_csv('cachedata/age_data.txt',sep=' ')
	
	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		loc=int(line.snvLoc)
		ploc=int(line.ploc)
		
		cache=Phosphorylation[(Phosphorylation.ACC_ID==uni)&(Phosphorylation.siteLoc==ploc)]
		if cache.shape[0]>0:
			if pd.notna(cache.MS_LIT.iloc[0]):
				df.loc[i,'feature_pubmed_MS_LIT']=int(cache.MS_LIT.iloc[0])
			else:
				df.loc[i,'feature_pubmed_MS_LIT']=cache.MS_LIT.iloc[0]
			
			if pd.notna(cache.LT_LIT.iloc[0]):
				df.loc[i,'feature_pubmed_LT_LIT']=int(cache.LT_LIT.iloc[0])
			else:
				df.loc[i,'feature_pubmed_LT_LIT']=cache.LT_LIT.iloc[0]
		else:
			df.loc[i,'feature_pubmed_MS_LIT']=np.nan
			df.loc[i,'feature_pubmed_LT_LIT']=np.nan
		
		cache=age[(age.acc==uni)&(age.position==ploc)]
		if cache.shape[0]>0:
			df.loc[i,'feature_psite_ptmage0']=cache.w0_mya.iloc[0]
			df.loc[i,'feature_psite_ptmage3']=cache.w3_mya.iloc[0]
		else:
			df.loc[i,'feature_psite_ptmage0']=np.nan
			df.loc[i,'feature_psite_ptmage3']=np.nan
	return df



## Calculate co-evolutionary information between the phosphorylation site and nsSNV
def readFastaEntry(fasta_file, rm_paralog=False, pro_id=""):
	# open file and read lines
	fid = open(fasta_file,"r")
	all_lines = fid.readlines()
	fid.close()

	# process all lines
	names, species, sequences = [], [], []
	seq = ""
	for i in range(len(all_lines)):
		line = all_lines[i].split()[0]
		if line.startswith( ">" ):
			names.append(line.split(">")[1])
			species.append(line.split(">")[1].split(".")[0])
			if i == 0:
				continue
			sequences.append(seq)
			seq = ""
		else:
			seq = seq + line
	sequences.append(seq)
	names, sequences = np.array(names), np.array(sequences)
	if rm_paralog == False:
		return names, sequences

	# remove paralog data
	human_idx = np.where(names == (pro_id))[0][0]
	human_seq = []
	human_seq[:0] = str(sequences[human_idx])
	species_uni = np.unique(species)
	kp_idx = np.zeros(species_uni.shape[0],"int")
	for i in range(species_uni.shape[0]):
		if species_uni[i] == "9606":
			kp_idx[i] = human_idx
		else :
			# choose the most similar paralogous protein via hamming distance
			_idx = np.where(np.array(species)==species_uni[i])[0]
			if _idx.shape[0] > 0:
				hm_loss = np.zeros(_idx.shape[0],"float")
				for j in range(_idx.shape[0]):
					seq_tmp = []
					seq_tmp[:0] = str(sequences[_idx[j]])
					hm_loss[j] = hamming_loss(human_seq, seq_tmp)
				_idx = _idx[np.argsort(hm_loss)]
			kp_idx[i] = _idx[0]
			# kp_idx[i] = species.index(species_uni[i])
	return names[kp_idx], sequences[kp_idx]


def get_align_site(pro_id, alignseq, site, res):
	'''obtain loc on NOG
	the pro_id is the uniport id here
	the site is the absolute site - 1'''
	
	RV = -1
	# load fasta with uniprot id
	
	if uniprotSeq.get(pro_id):
		seq=uniprotSeq[pro_id]
	else:
		print("There is no reference sequence with id as %s." %pro_id)
		return RV

	# check uniprot residue 
	if site>len(seq) or seq[site] != res:
		# print("The %d residue " %(site+1) + "on %s " %pro_id + "is not %s!" %res)
		return RV

	# map the residue to the aligned sequence
	seq_align = []
	seq_align[:] = alignseq
	seq_align = np.array(seq_align)
	idx_no_gap = np.where(seq_align != '-')[0]
	seq_align = seq_align[idx_no_gap]

	len_sur = 10
	left = min(len_sur, site)
	right = min(len(seq)-1-site, len_sur)

	for i in range(len(seq_align)):
		left_tmp = min(len_sur, i)
		right_tmp = min(len(seq_align)-1-i, len_sur)

		left_use = min(left, left_tmp)
		right_use = min(right, right_tmp)

		seq_sur_unip = []
		seq_sur_unip[:] = seq[site-left_use : site+right_use+1]

		seq_sur_alig = []
		seq_sur_alig[:] = seq_align[i-left_use : i+right_use+1]

		mapped_len = sum(np.array(seq_sur_unip) == np.array(seq_sur_alig))

		if seq[site] == seq_align[i] and mapped_len >= min(len(seq_sur_unip), len_sur+1):
			if RV != -1:
				print ("multiple matched surrouding sequence")
			RV = idx_no_gap[i]
	return RV

def get_msa_site(pro_id, seq, site, res):
	cnt = -1
	msa_site = -1
	for i in range(len(seq)):
		if seq[i] != "-" and seq[i] != "?":
			cnt += 1
		if cnt == site:
			msa_site = i
			break
	if msa_site == -1:
		print(str(site + 1) + " is out the length of %s." %pro_id)
	elif seq[msa_site] != res:
		#print("The %d site" %(site+1) + " of %s" %pro_id + " is %s." %seq[msa_site])
		msa_site = -1
#	 print(msa_site)
	return msa_site

def site_coevolve(unip_id, site1, site2, res1, res2, rm_paralog=True):
	"""
	standardize pro_id: Ensembl protein id; site1/2: int; res1/2: char
	output: save multiple sequence alignment for the two sites into a txt file
			return co-evolution score, float
			(optional)show the Phylogenetic Trees
	e.g. unip_id, site1, site2, res1, res2='P00533',np.array([1074,1074]),np.array([1070,1069]),np.array(['T','T']),np.array(['S','Y'])
	"""
	RV = np.zeros(site1.shape[0])
	RV[:] = None

	# load the id mapping file
	veNOG_Ensm_Unip_file = 'cachedata/eggNOG5/veNOG_Ensm_Unip_tri_40674.txt'
	veNOG_Ensm_Unip = np.loadtxt(veNOG_Ensm_Unip_file, delimiter='\t', skiprows=1, dtype="str")

	# id mapping and co-evolution
	idx_pro = np.where(veNOG_Ensm_Unip[:,2] == unip_id)[0]
	if idx_pro.shape[0] == 0:
		print("No MSA file in veNOG of protein %s!" %unip_id)
		return RV
	elif idx_pro.shape[0] > 1:
		if np.unique(veNOG_Ensm_Unip[idx_pro,0]).shape[0] == 1:
			pass
			# print("Multiple Ensembl ids in one veNOG file for %s!" %unip_id)
		else:
			print("Multiple veNOG files for %s!" %unip_id)
	# We will use the first one for simplisity
	NOG_id = veNOG_Ensm_Unip[idx_pro,0]
	Ensm_id = veNOG_Ensm_Unip[idx_pro,1]

	pro_id, NOG_id = Ensm_id[0], NOG_id[0]

	# load veNOG and Ensembl id map file
	align_dir = "cachedata/eggNOG5/40674/"

	# check the veNOG file
	if os.path.isfile(align_dir + NOG_id + ".trimmed_alg.faa") == False:
		print("No file of %s" %NOG_id + ".fa in the path %s" %align_dir)
		return RV

	# process the fasta file
	fasta_file = align_dir + NOG_id + ".trimmed_alg.faa"
	
	with open(fasta_file,'r') as f:
		if f.read().find('>9606')<0:
			print("No human MSA data in %s" %NOG_id)
			return RV
	
	names, seqs = readFastaEntry(fasta_file, rm_paralog, pro_id)
	
	# get the sequence of pro_id in use
	pro_idx = np.where(names == (pro_id))[0][0]
	human_seq = seqs[pro_idx]
	seq_len = len(human_seq) - human_seq.count("-")

	for i in range(site1.shape[0]):
		# check the site and residue in aligned sequences
#		 print(site1[i], type(site1[i]))
		msa_site1 = get_msa_site(pro_id, human_seq, site1[i]-1, res1[i])
		msa_site2 = get_msa_site(pro_id, human_seq, site2[i]-1, res2[i])

		if msa_site1 == -1 or msa_site2 == -1:
			msa_site1 = get_align_site(unip_id, human_seq, site1[i]-1, res1[i])
			msa_site2 = get_align_site(unip_id, human_seq, site2[i]-1, res2[i])

		# if msa_site1 == -1 or msa_site2 == -1:
		#	 print("The %dth or %dth site" %(site1[i]+1, site2[i]+1) + " of %s" %pro_id + " does not match with residue.")
		#	 continue

		unmatch = False
		if msa_site1 == -1:
			print("The %d residue %s of %s doesn't match in MSA file" %(site1[i], res1[i], pro_id))
			unmatch = True
		if msa_site2 == -1:
			print("The %d residue %s of %s doesn't match in MSA file" %(site2[i], res2[i], pro_id))
			unmatch = True
		if unmatch:
			continue

		# get the consevation sequences
		con_res1, con_res2 = [], []
		for j in range(len(names)):
			con_res1.append(seqs[j][msa_site1])
			con_res2.append(seqs[j][msa_site2])

		# get the site co-evolution score via normalized mutual information
		RV[i] = normalized_mutual_info_score(con_res1, con_res2)
#		 print(con_res1, con_res2)
	return RV

def site_coevolve(unip_id, site1, site2, res1, res2, rm_paralog=True):
	"""
	standardize pro_id: Ensembl protein id; site1/2: int; res1/2: char
	output: save multiple sequence alignment for the two sites into a txt file
			return co-evolution score, float
			(optional)show the Phylogenetic Trees
	e.g. unip_id, site1, site2, res1, res2='P00533',np.array([1074,1074]),np.array([1070,1069]),np.array(['T','T']),np.array(['S','Y'])
	"""
	RV = np.zeros(site1.shape[0])
	RV[:] = None

	# load the id mapping file
	veNOG_Ensm_Unip_file = 'cachedata/eggNOG5/veNOG_Ensm_Unip_tri_40674.txt'
	veNOG_Ensm_Unip = np.loadtxt(veNOG_Ensm_Unip_file, delimiter='\t', skiprows=1, dtype="str")

	# id mapping and co-evolution
	idx_pro = np.where(veNOG_Ensm_Unip[:,2] == unip_id)[0]
	if idx_pro.shape[0] == 0:
		print("No MSA file in veNOG of protein %s!" %unip_id)
		return RV
	elif idx_pro.shape[0] > 1:
		if np.unique(veNOG_Ensm_Unip[idx_pro,0]).shape[0] == 1:
			pass
			# print("Multiple Ensembl ids in one veNOG file for %s!" %unip_id)
		else:
			print("Multiple veNOG files for %s!" %unip_id)
	# We will use the first one for simplisity
	NOG_id = veNOG_Ensm_Unip[idx_pro,0]
	Ensm_id = veNOG_Ensm_Unip[idx_pro,1]

	pro_id, NOG_id = Ensm_id[0], NOG_id[0]

	# load veNOG and Ensembl id map file
	align_dir = "cachedata/eggNOG5/40674/"

	# check the veNOG file
	if os.path.isfile(align_dir + NOG_id + ".trimmed_alg.faa") == False:
		print("No file of %s" %NOG_id + ".fa in the path %s" %align_dir)
		return RV

	# process the fasta file
	fasta_file = align_dir + NOG_id + ".trimmed_alg.faa"
	
	with open(fasta_file,'r') as f:
		if f.read().find('>9606')<0:
			print("No human MSA data in %s" %NOG_id)
			return RV
	
	names, seqs = readFastaEntry(fasta_file, rm_paralog, pro_id)
	
	# get the sequence of pro_id in use
	pro_idx = np.where(names == (pro_id))[0][0]
	human_seq = seqs[pro_idx]
	seq_len = len(human_seq) - human_seq.count("-")

	for i in range(site1.shape[0]):
		# check the site and residue in aligned sequences
#		 print(site1[i], type(site1[i]))
		msa_site1 = get_msa_site(pro_id, human_seq, site1[i]-1, res1[i])
		msa_site2 = get_msa_site(pro_id, human_seq, site2[i]-1, res2[i])

		if msa_site1 == -1 or msa_site2 == -1:
			msa_site1 = get_align_site(unip_id, human_seq, site1[i]-1, res1[i])
			msa_site2 = get_align_site(unip_id, human_seq, site2[i]-1, res2[i])

		# if msa_site1 == -1 or msa_site2 == -1:
		#	 print("The %dth or %dth site" %(site1[i]+1, site2[i]+1) + " of %s" %pro_id + " does not match with residue.")
		#	 continue

		unmatch = False
		if msa_site1 == -1:
			print("The %d residue %s of %s doesn't match in MSA file" %(site1[i], res1[i], pro_id))
			unmatch = True
		if msa_site2 == -1:
			print("The %d residue %s of %s doesn't match in MSA file" %(site2[i], res2[i], pro_id))
			unmatch = True
		if unmatch:
			continue

		# get the consevation sequences
		con_res1, con_res2 = [], []
		for j in range(len(names)):
			con_res1.append(seqs[j][msa_site1])
			con_res2.append(seqs[j][msa_site2])

		# get the site co-evolution score via normalized mutual information
		RV[i] = normalized_mutual_info_score(con_res1, con_res2)
#		 print(con_res1, con_res2)
	return RV

def get_coevolve(df):
	coevD=[]
	for uni in df.uni.drop_duplicates():
		cache=df[(df.uni==uni)]

		site1=np.array(cache.snvLoc.astype('int'))
		res1=np.array(cache.AAb)

		site2=np.array(cache.ploc.astype('int'))
		res2=np.array(cache.psite)

		coevMI = site_coevolve(uni, site1, site2, res1, res2)

		t=0
		for i in cache.index:
			coevD.append([i,coevMI[t]])
			t+=1
			
	coevD=pd.DataFrame(coevD,columns=['index','site_coevolve'])
	coevD.index=list(coevD['index'])
	df=pd.merge(df,coevD,left_index=True,right_index=True,how='left')

	return df#.drop(columns=['index'])



## Obtain SIFT score
def get_sift(df):
	for i in df.index:
		line=df.loc[i]
		uni=line.uni
		AAp=line.psite
		ploc=int(line.ploc)
		siftIDp=(uni,ploc,AAp)
		if siftScore.get(siftIDp):
			score=siftScore[siftIDp]
			df.loc[i,'feature_sift_psite']=score[0]
			df.loc[i,'feature_provean_psite']=score[1]
	return df