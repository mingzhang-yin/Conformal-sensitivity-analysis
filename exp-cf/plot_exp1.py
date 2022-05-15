import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9

from plotnine import *
from glob import glob
from pathlib import Path
import plotnine

palette = ['#D55E00', '#E69F00', '#0072B2',
                  '#009E73', '#F0E442', '#CC79A7', '#56B4E9']


#%% 
#### Coverage plot #######
alpha=0.2
t="Coverage"
#t="Interval_length"

subdirs = glob("./out/homo/*/")

groups=['CSA-M','CSA-Q','ITE-NUC']
name=Path(subdirs[0]).parent.name
Ds = []
for group in groups:
    ls=[]
    gmms=[]
    lens=[]
    for sub in subdirs:
        gamma = float(Path(sub).name)
        gmms.append(gamma)
        if t=="Coverage":
            reader = pd.read_csv(sub+"coverage.csv")
        if t=="Interval_length":
            reader = pd.read_csv(sub+"len.csv")
        reader['gmm'] = gamma
        Ds.append(reader)
        
Data = pd.concat(Ds)
Data['gmm_cat'] = pd.Categorical(Data['gmm'], categories=list(sorted(gmms)))


p = ggplot(Data, aes(x='gmm_cat', y=t, fill='group')) +\
    geom_boxplot( alpha=0.7, outlier_alpha=0.5, outlier_stroke=0.2, width = 0.5) +\
    theme(figure_size=(10,6),\
    legend_position="top",\
        legend_entry_spacing_x=20,\
    text = element_text(size=18), \
    title=element_text(size=18),\
        legend_title = element_blank(),\
            legend_entry_spacing_y = 20, \
        legend_margin=-7
    ) +\
    labs(x=r'Confounding strength ($\Gamma$)',y=t) +\
        ylim([0.4,1.0]) +\
    scale_fill_manual({'CSA-M':palette[0], 'CSA-Q':palette[2], 'ITE-NUC':palette[3]})
if t=="Coverage":
    p += geom_hline(yintercept=1-alpha,linetype='dashed',size=1, color='red')

ggplot.draw(p)
plt.tight_layout()
#%%
#### Length plot #######
alpha=0.2
t="Interval-length"
SAVE = False
subdirs = glob("./out/homo/*/")


groups=['CSA-M','CSA-Q','ITE-NUC']
name=Path(subdirs[0]).parent.name
Ds = []
for group in groups:
    ls=[]
    gmms=[]
    lens=[]
    for sub in subdirs:
        gmms.append(float(Path(sub).name))
        if t=="Coverage":
            reader = pd.read_csv(sub+"coverage.csv")
        if t=="Interval-length":
            reader = pd.read_csv(sub+"len.csv")
        df = reader[reader['group'] == group]
        tmp = df.iloc[:,0].values
        ls.append(tmp)
        
    gmms = np.array(gmms)
    ls = np.array(ls)
    idx = gmms.argsort() 
    gmms = gmms[idx] 
    ls = ls[idx]
    if t=="Interval-length":
        d = {"gmm": gmms,"center": np.median(ls,axis=1), "min" : np.quantile(ls,0.,axis=1), \
             "max" : np.quantile(ls,1,axis=1), "group" : group }
    if t=="Coverage":
        d = {"gmm": gmms,"center": np.median(ls,axis=1), "min" : np.quantile(ls,0.,axis=1), \
             "max" : np.quantile(ls,1,axis=1), "group" : group}
    Ds.append(pd.DataFrame(data=d))
Data = pd.concat(Ds)


p = ggplot(Data,aes(x='gmm',y='center',group='group',colour='group')) + \
        geom_line(aes(group='group'),size=1.3) + \
        geom_point(size=4) + \
        geom_errorbar(aes(ymin='min', ymax='max'), width=.2, size=0.8) +\
        labs(x=r'Confounding strength ($\Gamma$)',y=t) +\
        theme(figure_size=(9,7),\
            legend_position='right',\
              subplots_adjust={'right': 0.8},\
            text = element_text(size=18), \
            title=element_text(size=18),\
                legend_title = element_blank(),\
                    legend_entry_spacing_y = 20
            ) +\
        scale_color_manual({'CSA-M':palette[0], 'CSA-Q':palette[2], 'ITE-NUC':palette[3]})
    
if t=="Coverage":
    p += geom_hline(yintercept=1-alpha,linetype='dashed',size=1, color='red')

ggplot.draw(p)


#%% 
#### Shrinkage plot #######

def get_fcts(subdirs, group, exp='exp1'):
    gmms=[]
    Median = []
    Max = []
    Min = []
 
    for gamma in np.arange(1, 4.5, 0.5):
        gmms.append(gamma)
        if gamma%1 ==0:
            gamma = int(gamma)
        sub = subdirs+"gmm_{}/".format(gamma)
        if exp =='exp1':
            t = 'fct_min'
        else:
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
exp='exp1'
subdirs = "./out/loop_fct/{}/".format(dtype)
Ds = []

for group in groups:
    print(group)
    d = get_fcts(subdirs, group, exp)
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




#%% Individual level plot

import seaborn as sns
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
data = pd.read_csv('./out/individual.csv')
 
gmms = np.unique(data.gmm)
groups = np.unique(data.group)
groups = groups[[2,0,1]]
nplot = 70

rr = 5
xpos = []
labelpos = np.arange(len(groups))*rr*2
for r in labelpos:
    xpos.extend(list(np.linspace(r-(rr-1), r+(rr-1), num=nplot)))

data_combined = []

for gmm in [gmms[0],gmms[-3]]:

    subdata = []
    for group in groups:
        dta = data[(data['group']==group)&(data['gmm']==gmm)]
        subdata.append(dta.sample(nplot,replace=False))
    subdata = pd.concat(subdata, ignore_index=True)
    
    subdata['lower'] = subdata['lower']-subdata['ite']
    subdata['upper'] = subdata['upper']-subdata['ite']
    
    subdata['group_cat'] = pd.Categorical(subdata['group'], categories=list(groups))
    
    
    subdata = subdata.groupby(['group'], sort=False) \
        .apply(lambda x: x.sort_values(['cover'], ascending=False)) \
        .reset_index(drop=True)
    
    subdata['cover'] = pd.Categorical(subdata['cover'],categories=[1,0])
    
    subdata['xpos'] = np.array(xpos)
    
    subdata = subdata.round(1)
    
    data_combined.append(subdata)

data_combined = pd.concat(data_combined, ignore_index=True)

def col_func(s):
    return r'$\Gamma=$ ' + str(s) 

p = ggplot(data_combined,aes(x='xpos', color='group_cat',linetype='cover')) + \
            geom_errorbar(mapping=aes(ymin='lower', ymax='upper'), \
                           width=.2, size=0.9) + \
            labs(x='', y=r'$\hat{C}(X_i) - \tau_i$') +\
            theme(figure_size=(15,8),\
               subplots_adjust={'right': 0.8},\
               text = element_text(size=18), \
               title=element_text(size=18),\
               legend_position='none',\
               legend_title = element_blank(),\
               legend_entry_spacing_y = 20) +\
            scale_x_continuous(breaks=list(labelpos), labels=['CSA-M', 'CSA-Q', 'ITE-NUC']) +\
            scale_linetype_manual(values=['solid','dashed']) +\
            ylim([-5,5]) +\
            geom_hline(yintercept=0,linetype='dotted',size=1,color='red') +\
            facet_grid('gmm ~ .',labeller=labeller(rows=lambda s: r'$\Gamma=$ '+str(s) )) +\
            scale_colour_manual({'Mean':palette[0], 'CQR':palette[2], 'ITE':palette[3]})

ggplot.draw(p)








