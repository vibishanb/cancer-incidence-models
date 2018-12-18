
# coding: utf-8

# In[2]:


import numpy
import numpy.random as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm_notebook as tqdm
sns.set()


# In[3]:


sns.set_palette(sns.color_palette("hls", 20))


# In[4]:


parr = numpy.exp(numpy.arange(-20, -5))
narr = numpy.exp(numpy.arange(15, 35))


# In[67]:


fig, axarr = plt.subplots(1, 3, figsize=(16,5), sharey=True)

for ax, t in zip(axarr, [2, 5, 10]):
    for i, l in zip(parr, parr):
        ax.plot(numpy.log10(narr), 1-(1-(i**t))**narr, '+-', label=numpy.format_float_scientific(l, precision=2))
        ax.set_title('Oncogenic mutations = '+ str(t), color='k')

axarr[1].set_xlabel('log(Cell number)')
axarr[0].set_ylabel(r'$p_{can}$')
plt.legend(fontsize='medium', title='Mutation rate', bbox_to_anchor=(0.95, 0.95), edgecolor='k', framealpha=0.5)
plt.tight_layout()
plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/mutation-probability-31Oct.svg')


# In[5]:


p = numpy.median(parr)
n = numpy.median(narr)


# In[5]:


pcan_n1 = (1-(1-(p**2))**narr)
pcan_p1 = (1-(1-(parr**2))**n)

pcan_n2 = (1-(1-(p**5))**narr)
pcan_p2 = (1-(1-(parr**5))**n)


# In[9]:


p_age_n1 = numpy.array([1-(1-pcan_n1)**i for i in tqdm(range(100), desc='n')])
p_age_p1 = numpy.array([1-(1-pcan_p1)**i for i in tqdm(range(100), desc='p')])

p_age_n2 = numpy.array([1-(1-pcan_p2)**i for i in tqdm(range(100), desc='n')])
p_age_p2 = numpy.array([1-(1-pcan_p2)**i for i in tqdm(range(100), desc='p')])


# In[29]:


fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(10,8), sharey=True, sharex=True)

axarr[0,0].plot(p_age_n1)
axarr[0,0].set_title('Oncogenic mutations = 2')
axarr[0,1].plot(p_age_n2)
axarr[0,1].set_title('Oncogenic mutations = 5')
axarr[0,1].legend(numpy.array([numpy.format_float_scientific(i, precision=2) for i in narr]), fontsize='x-small', title='Cell number', bbox_to_anchor=(0.7, 0.1), framealpha=0.5)
axarr[1,0].plot(p_age_p1)
axarr[1,1].plot(p_age_p2)
axarr[1,1].legend(numpy.array([numpy.format_float_scientific(i, precision=2) for i in parr]), fontsize='x-small', title='Mutation rate', bbox_to_anchor=(0.7, 0.1), framealpha=0.5)

fig.text(0.025, 0.5, r'$p_{can}$', fontsize=20, rotation='vertical')
fig.text(0.5, 0.025, 'Time (years)', fontsize=15, rotation='horizontal', ha='center')
# plt.tight_layout()
plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/pcan-time-31Oct.svg')


# In[201]:


y = p_age_n1.T[-1]
x = numpy.arange(10000)
dy = numpy.zeros(y.shape, numpy.float)
dy[0:-1] = numpy.diff(y)/numpy.diff(x)
dy[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])


# In[17]:


plt.plot(numpy.diff(p_age_n1, n=1))


# In[252]:


fig, axarr = plt.subplots(nrows=1, ncols=4, figsize=(17, 4), sharey=True)

for p, ax in zip([parr[8], parr[10], parr[-2], parr[-1]], axarr):
    pmut = 1-(1-(p**5))**narr[4:]
    p_age = numpy.array([1-(1-i)**numpy.arange(time) for i in pmut])
    ax.plot(p_age.T)
    ax.set_title('p = ' + numpy.format_float_scientific(p, precision=2))
    
plt.tight_layout()
axarr[0].set_ylabel(r'$p_{can}$')
axarr[0].legend(numpy.array([numpy.format_float_scientific(i, precision=2) for i in narr[4:]]), title='Cell number', fontsize='small', framealpha=0.5)
plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/pcan-with-age-2Nov.svg')

