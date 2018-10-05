import pandas as pd
from lifelines import CoxPHFitter

df = pd.read_csv("../data/data.txt",  sep='\t', index_col=0, header=None)
dft = df.transpose()
dft = dft.iloc[:,4:]  # remove non-numeric columns

cph = CoxPHFitter()
cph.fit(df1, duration_col='LivingDays', event_col='Status', show_progress=True)
