
# coding: utf-8

# In[18]:


""" Import of modules """
import numpy.random as np
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import math
import pandas as pd
from tqdm import tqdm_notebook as tqdm
# from matplotlib import cm
import time
sns.set()

""" Initialization """
Npop=100000 #Population size
parr = numpy.exp(numpy.arange(-24, -14)) #Growth rate
narr = numpy.exp(numpy.arange(14, 25)) #Cell number
garr = numpy.array([0.03125, 0.0625, 0.125, 0.25, 0.5])
# rep = 100

ndiv=10 #Number of cell divisions per year
age=150 #Lifespan
time=ndiv*age #Duration of the simulation
wts=numpy.array([0.013818048555027355, 0.0553159434123515, 0.07253241028642805, 0.07303103455912367, 0.07216711636515384, 0.06647847243710951, 0.06452984736662379, 0.07104508339877749, 0.08076197744193335, 0.08185169462960405, 0.07211714069611326, 0.06271758577923968, 0.048454493422914295, 0.038794489715138394, 0.034263609991378986, 0.03177168658747205, 0.026997999757072733, 0.017842803104216928, 0.01550856249432107]) #Weights from the US 2000 standard population
n_class = len(wts)
crc, cmc, crr = numpy.zeros(len(parr)*age).reshape(len(parr),age), numpy.zeros(len(parr)*age).reshape(len(parr),age), numpy.zeros(len(parr)*age).reshape(len(parr),age)
# marr = numpy.zeros(len(parr)*time).reshape(len(parr),time)
# aa_rate = numpy.zeros(len(parr))

""" Main simulation """
# np_flag = 0
# for p in tqdm(garr, desc='Growth rate', leave=False):

flag = 0
g_inc = 0.5
n = narr[6]
# p = parr[6]
threshold = 5 #Mutation threshold for cancer

for p in tqdm(parr, desc='Mutation rate', leave=False):
    cancer_count = numpy.zeros(age) #Age-wise incidence of cancer
    num_surv = numpy.zeros(age) #Number of survivors in each age/generation
    cancer_fract = numpy.zeros(age) #Normalized incidence of cancer
    cumul_count = numpy.zeros(age) #Cumulative count of cancer
    crude_rate = numpy.zeros(age) #Calculated age-wise incidence per 100000
    g = numpy.array([0.1, 0.1+g_inc, 0.1+(g_inc*2), 0.1+(g_inc*3), 0.1+(g_inc*4), 0.1+(g_inc*5)]) #Growth rate progression
    d = g_inc/10 #Death rate
    #     cancer_time = numpy.zeros(Npop)
    #         g_initial = []
    # age_rate=numpy.zeros((19)) #Age-specific rates adjusted to the US 2000 standard population

    for j in tqdm(range(Npop), desc='Npop', leave=False):
        t=0 #Index to track time
        n_mut=numpy.zeros(time, dtype=int)
        m=numpy.zeros(threshold+1) #Initial cell populations
        m[0]=n
        p_mut=1-((1-p)**m[0]) #Initial probabiltiy of first mutation arising in the population

        for t in range(time):

            n_mut[t]=n_mut[t-1 or 0]

            if (p_mut > np.random_sample())*(n_mut[t] < 5): #New mutant population
                n_mut[t] += 1
                m[n_mut[t]] = 1.0
                p_mut = 1-((1-p)**m[n_mut[t]])
                m[n_mut[t]-1] -= 1.0

            elif n_mut[t] > 0: #Growth of existing mutant or normal population, as the case may be
                m += ((m*g*(n-m.sum())/n) - m*d)
    #             n += (n*1*(k-n-m.sum())/k) - n*d
                p_mut = 1-((1-p)**m[n_mut[t]])


            if n_mut[t] == 5: #Recording actual cancer cases
                cancer_count[int(t/ndiv)] += 1
    #                 cancer_time[j] = int(t/ndiv)
                break

    """ Calculations """
    cumul_count = cancer_count.cumsum()

    num_surv = numpy.array([Npop]*age, dtype=float)
    num_surv[1:] -= cumul_count[:-1]

    cancer_fract = cancer_count/(cancer_count+num_surv)
    crude_rate = cancer_fract*100000

    # cc_bin, crr_bin, cf_bin, age_rate = numpy.zeros(n_class), numpy.zeros(n_class), numpy.zeros(n_class), numpy.zeros(n_class)
    # age_rate[0]=cancer_count[0]
    # age_rate[1]=sum(cancer_count[1:4])
    # age_rate[-1]=sum(cancer_count[85:len(cancer_count)])
    # for i in range(2,18):
    #     age_rate[i]=sum(cancer_count[(5*(i-1)):(5*(i-1)+4)])

    # numsurv_bin = numpy.array([Npop]*n_class,dtype=float)
    # cc_bin = age_rate.cumsum()
    # numsurv_bin[1:] -= cc_bin[:-1]

    # cf_bin = age_rate/(age_rate+cc_bin)
    # crr_bin = cf_bin*100000
    
    crc[flag] = cancer_count[:]
    cmc[flag] = cumul_count[:]
    crr[flag] = crude_rate[:]
    flag += 1


# In[14]:


cdata = pd.DataFrame(cmc.T)


# In[11]:


cdata.index.names = ['Replicates']
cdata.columns.names = ['Age']


# In[59]:


x = numpy.arange(0,age)
y = cmc.mean(axis=0)
ystd = cmc.std(axis=0)
conf = 1.96*ystd/(rep**0.5)


# In[60]:


plt.plot(x, y)
plt.fill_between(x, y-conf, y+conf, alpha=0.5)
plt.xlabel('Age (years)')
plt.ylabel('Cumulative incidence')


# In[73]:


mdata = pd.DataFrame(marr.T)
mdata.columns = garr
mdata.columns.names = ['Growth rate']
mdata.index.names = ['m']


# In[78]:


ax = sns.boxplot(data=(mdata/n).where(mdata!=0))
ax = sns.stripplot(data=(mdata/n).where(mdata!=0), color='.2', alpha=0.5, marker='.')
plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/22Aug2018/linear_v1_mnratio_gsvg')


# In[74]:


# mdata.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/22Aug2018/linear_v1_m_g.xlsx')

df_crc = pd.DataFrame(crc, index = garr)
df_crc.index.names = ['Growth rate']
df_crc.columns.names = ['Age']
df_crc.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/22Aug2018/linear_v1_raw_counts_g.xlsx')

df_cc = pd.DataFrame(cmc, index = garr)
df_cc.index.names = ['Growth rate']
df_cc.columns.names = ['Age']
df_cc.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/22Aug2018/linear_v1_cumulative_g.xlsx')

df_crr = pd.DataFrame(cmc, index = garr)
df_crr.index.names = ['Growth rate']
df_crr.columns.names = ['Age']
df_crr.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/22Aug2018/linear_v1_crude_g.xlsx')

# df_aar = pd.DataFrame(aa_rate, index=garr)
# df_aar.index.names = ['Cell number']
# df_aar.columns.names = [r'$Age-adjusted\ rate$']
# df_aar.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/22Aug2018/linear_v1_aa_rate_neff.xlsx')


# In[31]:


m


# In[34]:


m*g*(n-m.sum())/n - m*d


# In[20]:


half_age_parr = numpy.zeros_like(parr)
for i in range(len(parr)):
    half_age_parr[i] = numpy.less_equal(cmc[i], cmc[i,-1]/2).sum()


# In[29]:


plt.plot(numpy.log10(parr), cmc[:,-1], 'o-')
plt.xlabel('log(p)')
plt.ylabel(r'$I_{max}$')
# plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/22Aug2018/linear_v1_maxcumulative_g.svg')


# In[30]:


plt.plot(numpy.log10(parr), half_age_parr, 'o-')
plt.xlabel('log(p)')
plt.ylabel(r'$Age\ at\ \frac{I_{max}}{2}$')
# plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/22Aug2018/linear_v1_halfmax_g.svg')


# In[27]:


for i,l  in zip(cmc, numpy.log10(parr).round(decimals=2)):
    plt.plot(i, label=l)
plt.legend()
plt.xlabel('Age')
plt.ylabel('Cumulative incidence')
# plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/22Aug2018/linear_v1_g.svg')


# ## Visualising data

# ### Age-adjusted rates

# In[4]:


fig, axarr = plt.subplots(1, 3, figsize=(12,4), sharey=True)
nl = [r'$n_{1}$', r'$n_{2}$', r'$n_{3}$', r'$n_{4}$']
pl = [r'$p_{1}$', r'$p_{2}$', r'$p_{3}$', r'$p_{4}$']
# ax.set_xticks(xs)
# ax.set_yticks(ys)
# ax.set_xticklabels(nl)
# ax.set_yticklabels(pl)
# aar_g = [aar_g1, aar_g2, aar_g3]
titles = [r'$g = 0.05$', r'$g = 0.1$', r'$g = 0.2$']
aar_g = aar_g.reshape(3,4,4)

for ax, p, title in zip(axarr, aar_g, titles):
    im = ax.imshow(p, interpolation=None, cmap=cm.coolwarm)
    ax.set_xticks(numpy.arange(4))
    ax.set_yticks(numpy.arange(4))
    ax.set_yticklabels(nl)
    ax.set_xticklabels(pl)
    ax.set_xlabel(r'$Mutation\ rate$')
    ax.set_title(title)
#     ax.set_xbound(lower=10**-8, upper=10**-6)
#     ax.set_ybound(lower=10**6, upper=10**8)
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', shrink=0.85, pad=0.2)
#     plt.setp(cbar.get_ticks(), size=8)
axarr[0].set_ylabel(r'$Cell\ number$')
plt.tight_layout()
plt.savefig('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/25Jun2018/linear_v1_all_effects_aa_rate.svg')
# fig.suptitle(r'$Effect\ of\ g$', horizontalalignment='center', y=1.05, fontsize=15)

# ax = fig.add_subplot(3, 1, 2)
# axarr[1].plot(xs, aar_g2.T, label=['n1', 'n2', 'n3', 'n4', 'n5'])

# ax = fig.add_subplot(3, 1, 3)
# axarr[2].plot(xs, aar_g3.T, label=['n1', 'n2', 'n3', 'n4', 'n5'])


# ### Crude rate per 100000

# In[7]:


fig, axarr = plt.subplots(3, 4, sharex=True, sharey=False, figsize=(15,9))
j=0
age_classes = [r'$0$', r'$1-4$', r'$5-9$', r'$10-14$', r'$15-19$', r'$20-24$', r'$25-29$', r'$30-34$', r'$35-39$', r'$40-44$', r'$45-49$', r'$50-54$', r'$55-59$', r'$60-64$', r'$65-69$', r'$70-74$', r'$75-79$', r'$80-84$', r'$85-89$']

fig.text(0.5, -0.03, r'$Age\ classes\ (yrs)$', ha='center', fontsize=20)
fig.text(-0.05, 0.75, r'$g=0.05$')
fig.text(-0.05, 0.5, r'$g=0.1$')
fig.text(-0.05, 0.25, r'$g=0.2$')
fig.text(-0.1, 0.5, r'$Crude\ rate\ per\ 100000$', va='center', rotation='vertical', fontsize=20)

for a in axarr[-1]:
#     a.set_xlabel('Age class')
    a.set_xticks(numpy.arange(19))
    a.set_xticklabels(age_classes)
    plt.setp(a.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=5)

for ax, garr in zip(axarr, crr_g):
    for i in range(len(ax)):
        ga = garr[i*4:(i+1)*4]
        for k in range(4):    
            ax[i].plot(numpy.arange(n_class), ga[k].T, label=pl[k])
            if j == 0:
                ax[i].set_title(nl[i])
    if j == 0:
        ax[-1].legend()
    j += 1

# fig.legend(p_labels, bbox_to_anchor=(0.6,0.02), ncol=5)
plt.tight_layout()
plt.savefig('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/25Jun2018/linear_v1_all_effects_crude_rate.svg')


# ### Cumulative incidence

# In[8]:


fig, axarr = plt.subplots(3, 4, sharex=True, sharey=False, figsize=(15,9))
j=0
age_classes = [r'$0$', r'$1-4$', r'$5-9$', r'$10-14$', r'$15-19$', r'$20-24$', r'$25-29$', r'$30-34$', r'$35-39$', r'$40-44$', r'$45-49$', r'$50-54$', r'$55-59$', r'$60-64$', r'$65-69$', r'$70-74$', r'$75-79$', r'$80-84$', r'$85-89$']

fig.text(0.5, -0.03, r'$Age\ classes\ (yrs)$', ha='center', fontsize=20)
fig.text(-0.05, 0.75, r'$g=0.05$')
fig.text(-0.05, 0.5, r'$g=0.1$')
fig.text(-0.05, 0.25, r'$g=0.2$')
fig.text(-0.1, 0.5, r'$Cumulative\ cancer\ incidence$', va='center', rotation='vertical', fontsize=20)

for a in axarr[-1]:
#     a.set_xlabel('Age class')
    a.set_xticks(numpy.arange(19))
    a.set_xticklabels(age_classes)
    plt.setp(a.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=5)

for ax, garr in zip(axarr, cc_g):
    for i in range(len(ax)):
        ga = garr[i*4:(i+1)*4]
        for k in range(4):    
            ax[i].plot(numpy.arange(n_class), ga[k].T, label=pl[k])
            if j == 0:
                ax[i].set_title(nl[i])
    if j == 0:
        ax[-1].legend()
    j += 1

# fig.legend(p_labels, bbox_to_anchor=(0.6,0.02), ncol=5)
plt.tight_layout()
plt.savefig('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/figures/25Jun2018/linear_v1_all_effects_cumul_rate.svg')


# ### Data export

# #### Crude rate

# In[9]:


df_crr_g1 = pd.DataFrame(crr_g[0], index = [numpy.ravel([[x,x,x,x] for x in garr]), numpy.array([garr]*4).ravel()])
df_crr_g2 = pd.DataFrame(crr_g[1], index = [numpy.ravel([[x,x,x,x] for x in garr]), numpy.array([garr]*4).ravel()])
df_crr_g3 = pd.DataFrame(crr_g[2], index = [numpy.ravel([[x,x,x,x] for x in garr]), numpy.array([garr]*4).ravel()])

df_crr_g1.index.names = ['Growth rate', 'Cell number']
df_crr_g1.columns.names = ['Age classes']
df_crr_g1.columns = age_classes
df_crr_g1.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_crude_rate_g_0.05.xlsx')

df_crr_g2.index.names = ['Growth rate', 'Cell number']
df_crr_g2.columns.names = ['Age classes']
df_crr_g2.columns = age_classes
df_crr_g2.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_crude_rate_g_0.1.xlsx')

df_crr_g3.index.names = ['Growth rate', 'Cell number']
df_crr_g3.columns.names = ['Age classes']
df_crr_g3.columns = age_classes
df_crr_g3.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_crude_rate_g_0.2.xlsx')


# #### Cumulative incidence

# In[10]:


df_cc_g1 = pd.DataFrame(cc_g[0], index = [numpy.ravel([[x,x,x,x] for x in garr]), numpy.array([garr]*4).ravel()])
df_cc_g2 = pd.DataFrame(cc_g[1], index = [numpy.ravel([[x,x,x,x] for x in garr]), numpy.array([garr]*4).ravel()])
df_cc_g3 = pd.DataFrame(cc_g[2], index = [numpy.ravel([[x,x,x,x] for x in garr]), numpy.array([garr]*4).ravel()])

df_cc_g1.index.names = ['Growth rate', 'Cell number']
df_cc_g1.columns.names = ['Age classes']
df_cc_g1.columns = age_classes
df_cc_g1.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_cumulative_rate_g_0.05.xlsx')

df_cc_g2.index.names = ['Growth rate', 'Cell number']
df_cc_g2.columns.names = ['Age classes']
df_cc_g2.columns = age_classes
df_cc_g2.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_cumulative_rate_g_0.1.xlsx')

df_cc_g3.index.names = ['Growth rate', 'Cell number']
df_cc_g3.columns.names = ['Age classes']
df_cc_g3.columns = age_classes
df_cc_g3.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_cumulative_rate_g_0.2.xlsx')


# #### Age-adjusted incidence

# In[11]:


# df_aar_g1 = pd.DataFrame(aar_g1, index = [numpy.ravel([[x,x,x,x,x] for x in p_arr]), numpy.array([n_arr]*5).ravel()])
# df_aar_g2 = pd.DataFrame(aar_g2, index = [numpy.ravel([[x,x,x,x,x] for x in p_arr]), numpy.array([n_arr]*5).ravel()])
# df_aar_g3 = pd.DataFrame(aar_g3, index = [numpy.ravel([[x,x,x,x,x] for x in p_arr]), numpy.array([n_arr]*5).ravel()])

# df_aar_g1 = df_aar_g1.unstack()
# df_aar_g2 = df_aar_g2.unstack()
# df_aar_g3 = df_aar_g3.unstack()
df_aar_g1 = pd.DataFrame(aar_g[0], index=garr, columns=garr)
df_aar_g2 = pd.DataFrame(aar_g[1], index=garr, columns=garr)
df_aar_g3 = pd.DataFrame(aar_g[2], index=garr, columns=garr)

df_aar_g1.index.names = ['Growth rate']
df_aar_g2.index.names = ['Growth rate']
df_aar_g3.index.names = ['Growth rate']

df_aar_g1.columns.names = ['Cell number']
df_aar_g2.columns.names = ['Cell number']
df_aar_g3.columns.names = ['Cell number']

df_aar_g1.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_aa_rate_g_0.05.xlsx')
df_aar_g2.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_aa_rate_g_0.1.xlsx')
df_aar_g3.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V1/data/25Jun2018/linear_v1_aa_rate_g_0.2.xlsx')

