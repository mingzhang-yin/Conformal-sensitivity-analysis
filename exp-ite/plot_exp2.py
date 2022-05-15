import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plotnine as p9
import os
from plotnine import *
from glob import glob
from pathlib import Path

palette = ['#D55E00', '#E69F00', '#0072B2','#009E73', '#F0E442', '#CC79A7', '#56B4E9']


#%%
#### Coverage plot #######
alpha=0.2
#t="Interval_length"
t="Coverage"
dtype = "hete"

subdirs = "./out/{}/".format(dtype)

groups=['CSA-M','CSA-Q','CSA-B', 'ITE_NUS']
Ds = []

gmms=[]
for gamma in np.arange(1, 4.5, 0.5):
    gmms.append(gamma)
    if gamma%1 ==0:
        gamma = int(gamma)
    sub = subdirs+"gmm_{}/".format(gamma)
    if t=="Coverage":
        reader = pd.read_csv(sub+"coverage.csv")
    if t=="Interval_length":
        reader = pd.read_csv(sub+"len.csv")
    reader['gmm'] = gamma
    Ds.append(reader)

Data = pd.concat(Ds)
Data['group'] = Data['group'].map({'CSA-B':'Bonf.', 'CSA-M':'CSA-M', 'CSA-Q':'CSA-Q', 'ITE-NUS':'ITE-NUC'})
Data['gmm_cat'] = pd.Categorical(Data['gmm'], categories=list(sorted(gmms)))


p = ggplot(Data, aes(x='gmm_cat', y=t, fill='group')) +geom_boxplot( alpha=0.7, outlier_alpha=0.5, outlier_stroke=0.2, width = 0.5) +theme(figure_size=(10,6),      legend_position="top",      legend_entry_spacing_x=20,      #subplots_adjust={'right': 0.8},\
      text = element_text(size=18), \
      title=element_text(size=18),\
      legend_title = element_blank(),\
      legend_entry_spacing_y = 20, \
      legend_margin=-7
) +\
labs(x=r'Confounding strength ($\Gamma$)',y=t) +\
scale_fill_manual({'CSA-M':palette[0], 'CSA-Q':palette[2], 'ITE-NUC':palette[3], 'Bonf.':palette[1]})
if t=="Coverage":

    p += geom_hline(yintercept=1-alpha,linetype='dashed',size=1, color='red')

ggplot.draw(p)
plt.tight_layout()
#%%
#### Length plot #######

alpha=0.2
t="Interval_length"
dtype = "hete"

subdirs = "./out/{}/".format(dtype)

groups=['CSA-M','CSA-Q','CSA-B', 'ITE-NUS']
Ds = []
gmms=[]
for gamma in np.arange(1, 4.5, 0.5):
    gmms.append(gamma)
    if gamma%1 ==0:
        gamma = int(gamma)
    sub = subdirs+"gmm_{}/".format(gamma)
    reader = pd.read_csv(sub+"len.csv")
    reader['gmm'] = gamma
    for group in groups:
        ls = reader.loc[reader['group'] == group].iloc[:,0].values
        d = {"gmm": gamma,"center": np.median(ls), "min" : np.quantile(ls,0.), \
             "max" : np.quantile(ls,1), "group" : group }
        Ds.append(pd.DataFrame(data=d, index=[0]))


Data = pd.concat(Ds)
Data['group'] = Data['group'].map({'CSA-B':'Bonf.', 'CSA-M':'CSA-M', 'CSA-Q':'CSA-Q', 'ITE-NUS':'ITE-NUC'})
Data['gmm_cat'] = pd.Categorical(Data['group'])


p = ggplot(Data,aes(x='gmm',y='center',group='gmm_cat',colour='gmm_cat')) + \
        geom_line(aes(group='gmm_cat'),size=1.3) + \
        geom_point(size=4) + \
        geom_errorbar(aes(ymin='min', ymax='max'), width=.2, size=0.8) +\
        labs(x=r'Confounding strength ($\Gamma$)',y='Interval length') +\
        theme(figure_size=(9,7),\
            legend_position='right',\
              subplots_adjust={'right': 0.8},\
            text = element_text(size=18), \
            title=element_text(size=18),\
                legend_title = element_blank(),\
                    legend_entry_spacing_y = 20
            ) +\
        scale_color_manual({'CSA-M':palette[0], 'CSA-Q':palette[2], 'ITE-NUC':palette[3], 'Bonf.':palette[1]})

ggplot.draw(p)





#%%
#### Shrinkage plot #######

def get_fcts(subdirs, group):
    gmms=[]
    Median = []
    Max = []
    Min = []
 
    for gamma in np.arange(1, 4.5, 0.5):
        gmms.append(gamma)
        if gamma%1 ==0:
            gamma = int(gamma)
        sub = subdirs+"gmm_{}/".format(gamma)

        t = 'fct_val'

        reader =  pd.read_csv(sub+"fct.csv")
        df = reader[reader['group'] == group]
        vals  = df[t].values

        for i in range(df.shape[0]):
            if vals[i]> 1:
                df[t].iloc[i] = 1
        df[t] = 100 - 100*df[t]
        mean = np.mean(df[t])
        std = np.std(df[t])
        Median.append(mean)
        
        
        Max.append(mean + std)
        Min.append(max(0, mean - std))
    d = {"gmm":gmms, "min":Min, "max":Max, 'median':Median}
    return d


dtype = "hete"
groups=["CSA-M","CSA-Q"]
subdirs = "./out/loop_fct/{}/".format(dtype)
Ds = []

for group in groups:
    print(group)
    d = get_fcts(subdirs, group)
    d['group'] = group
    Ds.append(pd.DataFrame(data=d))
Data = pd.concat(Ds)



p = ggplot(Data,aes(x='gmm',y='median',group='group',colour='group')) + \
        geom_point(size=4) + \
        geom_line(aes(group='group'),size=1.3) + \
        labs(x='Confounding strength ($\Gamma$)',y='Shrinkage Factor (%)') +\
        geom_line(aes(group='group'),size=1.3) +\
geom_errorbar(aes(ymin='min', ymax='max'), width=0.2, size=0.8) +\
        theme(figure_size=(9,7),\
            legend_position='right',\
            legend_key = element_blank(),\
            legend_background = element_blank(),\
              subplots_adjust={'right': 0.8},\
            text = element_text(size=18), \
            title=element_text(size=18),\
                legend_title = element_blank(),\
                    legend_entry_spacing_y = 20
            ) +\
        scale_fill_brewer(type="qual", palette="Accent") +\
            scale_color_brewer(type="qual", palette="Accent")+\
        scale_color_manual({'CSA-M':palette[0], 'CSA-Q':palette[2]})


ggplot.draw(p)
plt.tight_layout()






