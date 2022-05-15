import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import os
from plotnine import *
from glob import glob
from pathlib import Path

palette = ['#D55E00', '#E69F00', '#0072B2',
                  '#009E73', '#F0E442', '#CC79A7', '#56B4E9']

#%%

def load_files(gmm, alpha, method):
    record_ng = []
    record_ps = []
    for i in range(50):
        i= i+1
        reader = pd.read_csv(path + "alpha_{}/gmm_{}/ntrial_{}.csv".format(alpha, gmm, i))
        high = "{}_high".format(method)
        low = "{}_low".format(method)
        record_ng.append(np.mean(reader[high]<=0))
        record_ps.append(np.mean(reader[low] >= 0))
    return record_ps, record_ng



def get_alpha(gmm, method="cqr",pos = True):
    means = []
    variance = []
    for alpha in np.arange(0.05, 0.5, 0.05):
        alpha = int(alpha*100)/100
        record_ps, record_ng = load_files(gmm, alpha, method)
        if pos:
            means.append(np.mean(record_ps))
            variance.append(np.std(record_ps))
        else:
            means.append(np.mean(record_ng))
            variance.append(np.std(record_ng))   
    return means, variance



def get_gmm_alpha(alpha, gmm, method='mean', pos=True):

        record_ps, record_ng = load_files(gmm, alpha, method)
        if pos:
            return record_ps
           
        else:
            return record_ng


def get_df(method='mean', pos=True):
    percentage = []
    alpha_list = []
    gmm_list = []
    for gmm in np.arange(10, 32, 2):
        gmm = int(gmm)/10
        print(gmm)
        if gmm%1 ==0:
            gmm = int(gmm)
        for alpha in np.arange(0.1, 0.5, 0.05):
            alpha = int(alpha*100)/100
            percentage.append(np.mean(get_gmm_alpha(alpha, gmm, method, pos)))
            alpha_list.append(alpha)
            gmm_list.append(gmm)
    data = np.vstack([np.array(alpha_list),np.array(gmm_list), np.array(percentage)]).T
    df = pd.DataFrame(data, columns = ['alpha', 'gmm', 'percentage'])
    df['alpha'] = 1- df['alpha']
    return df


def make_plot(df, method, sign):
    g = (ggplot(df,p9.aes(x='alpha',y='gmm',fill='percentage' ))
         + geom_tile()
          + theme_minimal()
          + theme(figure_size=(6,6),\
                  text = element_text(size=18), \
                  title=element_text(size=18),\
                  axis_ticks_length=0,\
                  axis_ticks_length_major=0,\
                  axis_ticks_length_minor=0,
                  panel_grid_major = element_blank(), panel_grid_minor = element_blank()
                  )
          +labs(x=r'Coverage ($1-\alpha$)',y= r'Confounding strength ($\Gamma$)')
          + scale_x_continuous(expand = (0,0.01)) 
          + scale_y_continuous(expand = (0.02,0)) 
          + scale_fill_distiller(limits=(0,1),breaks=np.linspace(0,1,5), palette = "PuBu")
        ) 
    if sign=='positive':
        g = g + theme(legend_position="none")
    ggplot.draw(g)

#%%
path = "./out/"

for method in ["mean", "cqr"]:
    for pos in [False,True]:
        df = get_df(method,pos)
        if pos:
            make_plot(df, method,'positive' )
        else:
             make_plot(df,method,'negative' )
            

#%%
alpha = 0.2
path = "./out/alpha_" + str(alpha)
gmms = [1,2]
nplot = 70
labels = np.zeros([nplot,2])
reader = pd.read_csv("./out/alpha_" + str(alpha) + "/gmm_1/ntrial_1.csv")
n = np.shape(reader)[0]
labels[:,0] = np.sort(np.random.choice(n, nplot,replace=False))
labels[:,1] = np.random.choice(50, nplot,replace=True)

rr = 5
xpos = []
labelpos = np.arange(2)*rr*2
for r in labelpos:
    xpos.extend(list(np.linspace(r-(rr-1), r+(rr-1), num=nplot)))



record = []
for gmm in gmms:
    a = []
    for i in range(nplot):
        
        reader = pd.read_csv(path + "/gmm_{}/ntrial_{}.csv".format(gmm, int(labels[i,1]+1)))
        record.append([reader.loc[reader.index[int(labels[i,0])], 'mean_low'],
                       reader.loc[reader.index[int(labels[i,0])], 'mean_high'],
                       'CSA-M', gmm])
        a.append([reader.loc[reader.index[int(labels[i,0])], 'cqr_low'],
                       reader.loc[reader.index[int(labels[i,0])], 'cqr_high'],
                       'CSA-Q', gmm])
    record.extend(a)
data = pd.DataFrame(np.array(record))
data.columns =['lower', 'upper', 'group', 'gmm']

data['group_cat'] = pd.Categorical(data['group'])
data['lower'] = data['lower'].astype(float)
data['upper'] = data['upper'].astype(float)


data = data.sort_values(['gmm','group'], ascending=True) \
    .groupby(['gmm','group'], sort=False) \
    .apply(lambda x: x.sort_values(['lower'], ascending=True)) \
    .reset_index(drop=True)



data['xpos'] = np.array(xpos*2)

p = ggplot(data,aes(x='xpos', color='group')) + \
            geom_errorbar(mapping=aes(ymin='lower', ymax='upper'), \
                           width=.2, size=0.9) + \
            labs(x='', y=r'$\hat{C}(X_i)$') +\
            theme(figure_size=(15,8),\
               subplots_adjust={'right': 0.8},\
               text = element_text(size=18), \
               title=element_text(size=18),\
               legend_position='none',\
               legend_title = element_blank(),\
               legend_entry_spacing_y = 20) +\
            scale_x_continuous(breaks=list(labelpos), labels=['CSA-M', 'CSA-Q']) +\
            ylim([-2,6]) +\
            geom_hline(yintercept=0,linetype='dotted',size=1,color='red') +\
            facet_grid('gmm ~ .',labeller=labeller(rows=lambda s: r'$\Gamma=$ '+str(s) )) +\
            scale_colour_manual({'CSA-M':palette[0], 'CSA-Q':palette[2]})
ggplot.draw(p)





