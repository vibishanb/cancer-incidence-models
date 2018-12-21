"""
V1: the primary formulation of the random-chance-model as given in the PNAS-MS, with pcan = 1-(1-p**k)**n, and p-age = 1-(1-pcan)**age

V2 [started on 18th December]:
1. Alternate formulation for pcan, with n being a function of time; the time series comes from logistic growth assuming g = 0.007, d = g/5, and one step of logistic growth = 1 day of human lifespan.
Under this formulation, 

pcan(t) = 1-(1-p**k)**m_delta(t), where m_delta(t) is the net change in cell number at time = t, as given by one step of the discrete logistic growth equation.
p-age(A) = 1-(1-p**k)**(sum(m_delta) from t=0 to t=A); in NumPy terms, this is given as 1-(1-p**k)**m_delta.cumsum()

2. Also added is the randomization of p and n separately

V2 corrected on 22 December; m_delta was counting cell divisions incorrectly. It is the positive part of the discrete logisitic equation that reflects cell growth in every step.

"""

import numpy
import numpy.random as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm_notebook as tqdm
sns.set()
sns.set_palette(sns.color_palette("viridis", 20))

def get_m(n):
	g = 0.007 #Standardized growth rate for one logistic step = 1 day
	d = g/5 #Constant death rate
	time = 365*100 #Lifespan
	m, m_delta = numpy.zeros(time), numpy.zeros(time)
	m[0] = 1

	for t in range(1, time):
		m[t] = m[t-1]
		m_delta[t]= m[t]*g*(n-m[t])/n
		m_inc = m_delta[t] - m[t]*d
		m[t] += m_inc

	return m, m_delta


parr = numpy.exp(numpy.arange(-20, -5))
narr = numpy.exp(numpy.arange(15, 30))
p = numpy.median(parr)
n = numpy.median(narr)


pcan_n1 = numpy.array([(1-(1-(p**2))**(get_m(i)[1].cumsum())) for i in narr])
pcan_n2 = numpy.array([(1-(1-(p**5))**(get_m(i)[1].cumsum())) for i in narr])
pcan_p1 = numpy.array([(1-(1-(i**2))**(get_m(n)[1].cumsum())) for i in parr])
pcan_p2 = numpy.array([(1-(1-(i**5))**(get_m(n)[1].cumsum())) for i in parr])

#p_age_n1 = pcan_n1.cumsum(axis=1)
#p_age_n2 = pcan_n2.cumsum(axis=1)
#p_age_p1 = pcan_p1.cumsum(axis=1)
#p_age_p2 = pcan_p2.cumsum(axis=1)


fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(10,8), sharex=True)

axarr[0,0].plot(pcan_n1.T)
axarr[0,0].set_title('k = 2')
axarr[0,1].plot(pcan_n2.T)
axarr[0,1].set_title('k = 5')
axarr[0,1].legend(numpy.array([numpy.format_float_scientific(i, precision=2) for i in narr]), fontsize='x-small', title='Cell number', loc='best', framealpha=0.5)
axarr[1,0].plot(pcan_p1.T)
axarr[1,1].plot(pcan_p2.T)
axarr[1,1].legend(numpy.array([numpy.format_float_scientific(i, precision=2) for i in parr]), fontsize='x-small', title='Mutation rate', loc='best', framealpha=0.5)

fig.text(0.001, 0.5, r'$p_{A}$', fontsize=20, rotation='vertical')
fig.text(0.5, 0.0025, 'Time (days)', fontsize=15, rotation='horizontal', ha='center')
plt.tight_layout()
plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/22Dec/p-age-alternate-formulation-22Dec.svg')


# Randomizing p and n
pdist = numpy.exp(np.uniform(-20, -5, size=1000))
ndist = numpy.exp(np.uniform(15, 25, size=1000))
pc_prand_k2 = 1-(1-(1-(1-pdist**2)**n))**36500
pc_nrand_k2 = 1-(1-(1-(1-p**2)**ndist))**36500
pc_prand_k5 = 1-(1-(1-(1-pdist**5)**n))**36500
pc_nrand_k5 = 1-(1-(1-(1-p**5)**ndist))**36500

fig, axarr = plt.subplots(nrows=2, ncols=2, figsize=(10,8), sharey=True)

axarr[0,0].scatter(numpy.log10(ndist), pc_nrand_k2)
axarr[0,0].set_title('k = 2')
axarr[0,0].set_xlabel('log(n)')
axarr[0,1].scatter(numpy.log10(ndist), pc_nrand_k5)
axarr[0,1].set_title('k = 5')
axarr[0,1].set_xlabel('log(n)')
axarr[1,0].scatter(numpy.log10(pdist), pc_prand_k2)
axarr[1,0].set_xlabel('log(p)')
axarr[1,1].scatter(numpy.log10(pdist), pc_prand_k5)
axarr[1,1].set_xlabel('log(p)')
fig.text(0.001, 0.5, r'$p_{can}$', fontsize=20, rotation='vertical')
plt.tight_layout()
plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/18Dec/randomizing-p-and-n-18Dec.svg')


"""

fig, axarr = plt.subplots(nrows=1, ncols=4, figsize=(17, 4), sharey=True)

for p, ax in zip([parr[8], parr[10], parr[-2], parr[-1]], axarr):
    pmut = 1-(1-(p**5))**narr[4:]
    p_age = numpy.array([1-(1-i)**numpy.arange(time) for i in pmut])
    ax.plot(p_age.T)
    ax.set_title('p = ' + numpy.format_float_scientific(p, precision=2))
    
plt.tight_layout()
axarr[0].set_ylabel(r'$p_{can}$')
axarr[0].legend(numpy.array([numpy.format_float_scientific(i, precision=2) for i in narr[4:]]), title='Cell number', fontsize='small', framealpha=0.5)
#plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/pcan-with-age-2Nov.svg')

"""
# In[67]:

"""
fig, axarr = plt.subplots(1, 3, figsize=(16,5), sharey=True)

for ax, t in zip(axarr, [2, 5, 10]):
    for i, l in zip(parr, parr):
        ax.plot(numpy.log10(narr), 1-(1-(i**t))**narr, '+-', label=numpy.format_float_scientific(l, precision=2))
        ax.set_title('Oncogenic mutations = '+ str(t), color='k')

axarr[1].set_xlabel('log(Cell number)')
axarr[0].set_ylabel(r'$p_{can}$')
plt.legend(fontsize='medium', title='Mutation rate', bbox_to_anchor=(0.95, 0.95), edgecolor='k', framealpha=0.5)
plt.tight_layout()
#plt.savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/random-model-figures/mutation-probability-31Oct.svg')
"""