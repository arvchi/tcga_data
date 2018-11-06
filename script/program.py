import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
import numpy.random as npr
from scipy import stats

# Qvals function
def qvalues(p_values,pi_0_average):
   m=len(p_values)

   # Sort p-values
   p_values=p_values.sort_values(ascending=True)

   # calculate q(p[m])
   q_p_m = pi_0_average*max(p_values)

   # Calculate q-values
   c=0
   q_values=[]

   for i in range(m,0,-1):
       if q_values:
           y=q_values[c]
           c=c+1
       else:
           y=q_p_m

       x=(pi_0_average*m*p_values[(i-1)])/(i)
       q_values.append(min(x,y))

   return q_values

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

# changed range parameter from len(dft.columns) to 13476, to truncate DF
for i in range(2,13476):
    subset = [0,1,i] # create df with only status, LivingDays and one gene at the time
    dftsubset = dft.iloc[:,subset]
    cph.fit(dftsubset, duration_col='LivingDays', event_col='Status', show_progress=True, step_size=0.1)
    coxdata = coxdata.append(cph.summary,ignore_index=False) # only append doesnt work

# Store coxdata for presentation

    
    
# Extract p-values

p_vals = coxdata.loc[:,"p"]    



# Evaluate pi0 for different lambdas
m=len(p_vals)
lam=np.arange(0.75,0.95,0.01)
pi0_estimates=[]
for n in lam:
   pi0_estimates.append(sum(i > n for i in p_vals)/(m*(1-n)))
pi_0_average = np.mean(pi0_estimates)

# Get qvals
q_values=qvalues(p_vals,pi_0_average)

# number of prognostic expression levels as a function of q-value.
count = 0
for i in q_values:
    if i < 0.05:
        count+=1
print("nr of significant features:", count)
