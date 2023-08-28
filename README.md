# VIPpred: A novel model for predicting variant impact on phosphorylation events driving carcinogenesis
---

## Installation

Our model can be used to calculate features and predict VIP events. It is constructed sololy on Python 3. In order to run VIPpred on your terminal, it is imperative to have the following dependencies:
```
numpy, pandas, scipy, pickle, sklearn, xgboost
```
Should you lack any of the packages, install the corresponding package with
```
pip install XXX
```

In addition to direct download from GitHub, we have also uploaded VIPpred on Google Drive, and the link is:
```
https://drive.google.com/drive/folders/1FmEqDbwBL_hGjY89AZuyhx-e6xYierOT?usp=sharing
```
Due to the large size of databases such as AlphaFold2 and eggNOG5, which exceed the limitations of GitHub, we have uploaded the complete cache/ folder to Google Drive. Please download the cache folder and place it under the VIPpred directory.


The script uses dbNSFP (http://database.liulab.science/dbNSFP) to calculate functional scores of nsSNVs. Since the complete dbNSFP tool is over 30G, please download it from the link provided by dbNSFP, extract the files and place them under the “dbNSFP/” folder. The academic version currently used by VIPpred is dbNSFP v4.3a. Users should cite the tool properly as follows:
```
Liu X, Jian X, and Boerwinkle E. 2011. dbNSFP: a lightweight database of human non-synonymous SNPs and their functional predictions. Human Mutation. 32:894-899.
Liu X, Li C, Mou C, Dong Y, and Tu Y. 2020. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Medicine. 12:103.
```
If the user wish to use the tool for commerical use, it is imperitive to request a commerical license for dbNSFP.

After installing all dependencies and downloading & unzipping the code files of VIPpred, you should be ready to go!


## Running VIPpred on sample data

To start using VIPpred, make sure the input data file reside in VIPpred repositary folder. Run the following command in terminal:
```
python VIPpred.py Example.txt
```


## Sample file preview

```
uni snvLoc AAchange psite ploc chr chrl NNchange
Q14165 164 A/T S 162 12 120694899 G/A
P78364 708 P/L S 706 12 8934348 C/T
P78364 807 A/V Y 802 12 8936907 C/T
O60563 20 E/K Y 13 12 48716618 C/T
```
"uni" stands for UniProt ID.
"snvLoc" stands for nsSNV location on protein.
"AAchange" stands for amino acid substitution caused by nsSNV.
"psite" stands for amino acid type for phosphorylation site.
"ploc" stands for phosphorylation site location.
"chr" stands for the chromosome that nsSNV resides on.
"chrl" stands for the chromosome location that nsSNV resides on.
"NNchange" stands for nucleotide substitution of nsSNV.


## Contact

If you have any suggestions or questions, please visit http://bioinfolilab.phasep.pro/vippred for detailed information or open an issue on github.
