import pandas as pd
from lifelines import CoxPHFitter

df = pd.read_csv("../data/data.txt",  sep='\t', index_col=0, header=None)
dft = df.transpose()

d = {'dead': True, 'alive': False}
dft["Status"] = dft["Status"].map(d) # make event boolean

dft = dft.iloc[:,4:]  # only transcriptomic data as indepenent variablea
dft = dft.drop("Age",1)

# Make numeric (except event, which is boolean)
dft = dft.apply(pd.to_numeric, errors='coerce')

cph = CoxPHFitter()
cph.fit(dft, duration_col='LivingDays', event_col='Status', show_progress=True)

