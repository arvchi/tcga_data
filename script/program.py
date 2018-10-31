import pandas as pd
import numpy as np
from lifelines import CoxPHFitter

# Read data
df = pd.read_csv("../data/data.txt",  sep='\t', index_col=0, header=None)
dft = df.transpose()

# make event boolean
d = {'dead': True, 'alive': False}
dft["Status"] = dft["Status"].map(d) 

# Remove irrelevant fields                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
dft = dft.iloc[:,4:]  # only transcriptomic data as indepenent variablea
dft = dft.drop("Age",1)

# Remove genes with low convergence (manually)
drop = ["ENSG00000100146",'ENSG00000100219', 'ENSG00000100249', 
        'ENSG00000105371', 'ENSG00000120242', 'ENSG00000122718',
        'ENSG00000129864', 'ENSG00000129873', 'ENSG00000130283',
        'ENSG00000132631', 'ENSG00000139574', 'ENSG00000151079',
        'ENSG00000163157', 'ENSG00000164708', 'ENSG00000166160',
        'ENSG00000167822', 'ENSG00000168787', 'ENSG00000169484',
        'ENSG00000169488', 'ENSG00000170803', 'ENSG00000171484',
        'ENSG00000172199', 'ENSG00000172288', 'ENSG00000172352',
        'ENSG00000172381', 'ENSG00000172769', 'ENSG00000172774',
        'ENSG00000173349']
dft = dft.drop(drop,1)

# Make numeric (except event, which is boolean)
dft = dft.apply(pd.to_numeric, errors='coerce')

cph = CoxPHFitter()

# Make cox regression with one transcript at the time
coxdata = pd.DataFrame()

#dft.columns.get_loc("ENSG00000174942")
# change back starting point to "2".
# changed range paraemetar from len(dft.columns) to 13476, to truncate DF
for i in range(2,13476):
    subset = [0,1,i] # create df with only status, LivingDays and one gene
    dftsubset = dft.iloc[:,subset]
    cph.fit(dftsubset, duration_col='LivingDays', event_col='Status', show_progress=True, step_size=0.1)
    coxdata = coxdata.append(cph.summary,ignore_index=False) # only append doesnt work
