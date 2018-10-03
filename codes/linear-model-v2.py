
# coding: utf-8

# In[10]:


""" Import of modules """
import numpy.random as np
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
from tqdm import tqdm_notebook as tqdm

""" Initialization """
Npop = 100000 #Population size
parr = numpy.exp(numpy.arange(-24, -14))
narr = numpy.exp(numpy.arange(14, 25))
gmean = numpy.array([-1, -0.5, 0, 0.5, 1, 1.5, 2])
gvar = numpy.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4])
# rep = 100

ndiv = 10 #Number of cell divisions per year
age = 150 #Lifespan
time = ndiv*age #Duration of the simulation
wts = numpy.array([0.013818048555027355, 0.0553159434123515, 0.07253241028642805, 0.07303103455912367, 0.07216711636515384, 0.06647847243710951, 0.06452984736662379, 0.07104508339877749, 0.08076197744193335, 0.08185169462960405, 0.07211714069611326, 0.06271758577923968, 0.048454493422914295, 0.038794489715138394, 0.034263609991378986, 0.03177168658747205, 0.026997999757072733, 0.017842803104216928, 0.01550856249432107]) #Weights from the US 2000 standard population
n_class = len(wts)
crc, cmc, crr, marr = numpy.zeros(len(parr)*age).reshape(len(parr),age), numpy.zeros(len(parr)*age).reshape(len(parr),age), numpy.zeros(len(parr)*age).reshape(len(parr),age), numpy.zeros(len(parr)*time).reshape(len(parr),time)
# aa_rate = numpy.zeros(len(parr))

""" Main simulation """
# np_flag = 0
# for p in tqdm(gvar, desc='Var(g)', leave=False):
n = narr[6]
# p = parr[6]
threshold = 5
flag = 0
for p in tqdm(parr, desc='Mutation rate', leave=False):
    growth_rate = np.normal(0, 0.5, Npop) #Growth rate distribution
    cancer_count = numpy.zeros(age) #Age-wise incidence of cancer
    num_surv = numpy.zeros(age) #Number of survivors in each age/generation
    cancer_fract = numpy.zeros(age) #Normalized incidence of cancer
    cumul_count = numpy.zeros(age) #Cumulative count of cancer
    crude_rate = numpy.zeros(age) #Calculated age-wise incidence per 100000
    #     cancer_time = numpy.zeros(Npop)
    #         g_initial = []
    # age_rate=numpy.zeros((19)) #Age-specific rates adjusted to the US 2000 standard population

    for j in tqdm(range(Npop), desc='Npop', leave=False):
        t=0 #Index to track time
        n_mut=numpy.zeros(time, dtype=int)
        m=numpy.zeros(threshold+1) #Initial cell populations
        m[0]=n
        p_mut=1-((1-p)**m[0]) #Initial probabiltiy of first mutation arising in the population

        g_inc = math.fabs(growth_rate[j])
        d = g_inc/10
        g = numpy.array([0.1, 0.1+g_inc, 0.1+(g_inc*2), 0.1+(g_inc*3), 0.1+(g_inc*4), 0.1+(g_inc*5)]) #Growth rate progression

        for t in range(time):

            n_mut[t]=n_mut[t-1 or 0]

            if (p_mut > np.random_sample())*(n_mut[t] < 5): #New mutant population
                n_mut[t] += 1
                m[n_mut[t]] = 1.0
                p_mut = 1-((1-p)**m[n_mut[t]])
                m[n_mut[t]-1] -= 1.0

            elif (g < 0).any(): #Negative growth is physiologically undefined and mathematically unbounded
                break

            elif n_mut[t] > 0: #Growth of existing mutant or normal population, as the case may be
                m += ((m*g*(n-m.sum())/n) - m*d)
    #             n += (n*1*(k-n-m.sum())/k) - n*d
                p_mut = 1-((1-p)**m[n_mut[t]])
    #             m_inc = math.copysign(1.0, g)*math.fabs(m*g*(1-(m/n)))
    #             m += 0.05*m_inc

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

    # num_surv[0]=Npop
    # for t in range(1, age):
    #     num_surv[t]=num_surv[t-1]-cancer_count[t-1] #Number of individuals surviving at given age
    #     cumul_count[t]=cumul_count[t-1] + cancer_count[t]

    # t=0
    # while (t <= (age-1) and num_surv[t] != 0 ):
    #     cancer_fract[t]=cancer_count[t] / (cancer_count[t]+num_surv[t]) #Fraction of surviving population getting cancer
    #     crude_rate[t]=cancer_fract[t]*100000
    #     t+=1

    # numsurv_bin, cc_bin, crr_bin, cf_bin, age_rate = numpy.zeros(n_class), numpy.zeros(n_class), numpy.zeros(n_class), numpy.zeros(n_class), numpy.zeros(n_class)
    # age_rate[0]=cancer_count[0]
    # age_rate[1]=sum(cancer_count[1:4])
    # age_rate[-1]=sum(cancer_count[85:len(cancer_count)])
    # for i in range(2,18):
    #     age_rate[i]=sum(cancer_count[(5*(i-1)):(5*(i-1)+4)])

    # numsurv_bin[0] = Npop
    # for t in range(1, n_class):
    #     numsurv_bin[t] = numsurv_bin[t-1] - age_rate[t-1]
    #     cc_bin[t] = cc_bin[t-1] + age_rate[t]

    # #     cf_bin = age_rate/(age_rate+numsurv_bin)
    # #     crr_bin = cf_bin*10000
    # t = 0
    # while (t <= (n_class-1) and numsurv_bin[t] != 0):
    #     cf_bin[t] = age_rate[t] / (age_rate[t]+numsurv_bin[t])
    #     crr_bin[t] = cf_bin[t]*100000
    #     t += 1

    crc[flag] = cancer_count[:]
    cmc[flag] = cumul_count[:]
    crr[flag] = crude_rate[:]
#     aa_rate[flag] = (crr_bin*wts).sum()
    flag += 1


# In[2]:


sns.set()


# In[48]:


m[n_mut[t]-1]*g[n_mut[t]-1]*((k-n-m.sum())/k)


# In[12]:


plt.plot(cmc.T)


# In[13]:


plt.plot(n_mut)


# In[14]:


g


# In[114]:


m*g*((n-m.sum())/n) - m*d


# In[15]:


m


# In[11]:


cdata = pd.DataFrame(cmc.T)


# In[12]:


ax = sns.lineplot(data=cdata[0], ci=95, estimator=numpy.median, n_boot=100, err_style="band")


# In[6]:


mdata = pd.DataFrame(marr.T)
mdata.columns = gvar
mdata.columns.names = [r'$\sigma_{g}$']
mdata.index.names = ['m']


# In[7]:


ax = sns.boxplot(data=(mdata/n).where(mdata!=0))
ax = sns.stripplot(data=(mdata/n).where(mdata!=0), color='.2', alpha=0.5, marker='.', jitter=True)
plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/22Aug2018/linear_v2_mnratio_gvar.svg')

ax = sns.boxplot(data=mdata.where(mdata!=0))
ax = sns.stripplot(data=mdata.where(mdata!=0), color='.2', alpha=0.5, marker='.', jitter=True)
plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/22Aug2018/linear_v2_m_gvar.svg')


# In[3]:


plt.boxplot([marr[i][marr[i].nonzero()] for i in range(len(marr))], showbox=False, whis=[5,95])


# In[8]:


mdata.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/22Aug2018/linear_v2_m_gvar.xlsx')

df_crc = pd.DataFrame(crc, index = gvar)
df_crc.index.names = ['Var(g)']
df_crc.columns.names = ['Age']
df_crc.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/22Aug2018/linear_v2_raw_counts_gvar.xlsx')

df_cc = pd.DataFrame(cmc, index = gvar)
df_cc.index.names = ['Var(g)']
df_cc.columns.names = ['Age']
df_cc.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/22Aug2018/linear_v2_cumulative_gvar.xlsx')

df_crr = pd.DataFrame(crr, index = gvar)
df_crr.index.names = ['Var(g)']
df_crr.columns.names = ['Age']
df_crr.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/7Aug2018/linear_v2_crude_gvar.xlsx')

# df_aar = pd.DataFrame(aa_rate, index = gvar)
# df_aar.index.names = [r'$Mutation\ rate$']
# df_aar.columns.names = [r'$Age-adjusted\ rate$']
# df_aar.to_excel('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/2Aug2018/linear_v2_aa_rate_gvar.xlsx')


# In[9]:


half_age_gvar = numpy.zeros(rep)
for i in range(rep):
    half_age_gvar[i] = numpy.less_equal(cmc[i], cmc[i,-1]/2).sum()


# In[11]:


plt.plot(gvar, half_age_gvar, 'o-')
plt.xlabel(r'$\sigma_{g}$')
plt.ylabel(r'$Age\ at\ \frac{I_{max}}{2}$')
plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/22Aug2018/linear_v2_halfmax_gvar.svg')


# In[12]:


plt.plot(gvar, cmc[:,-1], 'o-')
plt.xlabel(r'$\sigma_{g}$')
plt.ylabel('Maximum cumulative incidence')
plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/22Aug2018/linear_v2_incidence_cumulmax_gvar.svg')


# In[13]:


for i,l  in zip(cmc, gvar):
    plt.plot(i, label=l)
plt.legend()
plt.xlabel('Age')
plt.ylabel('Cumulative incidence')
plt.savefig('/home/iiser/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/22Aug2018/linear_v2_cumulative_incidence_gvar.svg')


# ## Plotting

# In[6]:


index = cancer_time.nonzero()


# In[7]:


plt.plot(growth_rate[index], cancer_time[index], 'b.')


# In[8]:


plt.plot(gvar[index], cancer_time[index], 'b.')


# In[9]:


plt.plot(gvar[index], cancer_time[index], 'b.')


# In[10]:


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(gvar[index], gvar[index], cancer_time[index], marker='.')


# In[11]:


ns, ni = stats.linregress(gvar[index], cancer_time[index])[:2]
ps, pi = stats.linregress(gvar[index], cancer_time[index])[:2]
gs, gi = stats.linregress(growth_rate[index], cancer_time[index])[:2]

Enc = cancer_time[index] - (ns*gvar[index] + ni)
Epc = cancer_time[index] - (ps*gvar[index] + pi)
Egc = cancer_time[index] - (gs*growth_rate[index] + gi)


# In[18]:


plt.plot(gvar[index], Egc, 'b.')


# In[60]:


p_slope, p_int = stats.linregress(gvar[index], cancer_time[index])[:2]
Ept = cancer_time[index] - (p_slope*gvar[index] + p_int)


# In[61]:


plt.plot(gvar[index], Ept, 'b.')


# In[62]:


g_slope, g_int = stats.linregress(growth_rate[index], cancer_time[index])[:2]
Egt = cancer_time[index] - (g_slope*growth_rate[index] + g_int)


# In[64]:


plt.plot(gvar[index], Egt, 'b.')


# In[35]:


stats.spearmanr(growth_rate[index], cancer_time[index])


# In[36]:


stats.spearmanr(gvar[index], cancer_time[index])


# In[37]:


stats.spearmanr(gvar[index], cancer_time[index])


# ### Age-adjusted rates

# In[75]:


fig, ax = plt.subplots()
nl = [r'$10^{6}$', r'$10^{7}$', r'$10^{8}$', r'$10^{9}$', r'$10^{10}$']
pl = [r'$10^{-9}$', r'$10^{-8}$', r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$']

# ax.set_xticks(xs)
# ax.set_yticks(ys)
# ax.set_xticklabels(nl)
# ax.set_yticklabels(pl)
# aar_g = [aar_g1, aar_g2, aar_g3]
# titles = [r'$\sigma = 0.5$', r'$\sigma = 1$', r'$\sigma = 2$', r'$\sigma = 4$']

# for ax, p, title in zip(axarr, aar_g, titles):
im = ax.imshow(aa_rate.reshape(5,5), interpolation=None, cmap='jet', aspect='auto')
ax.set_xticks(numpy.arange(5))
ax.set_yticks(numpy.arange(5))
ax.set_yticklabels(nl)
ax.set_xticklabels(pl)
ax.set_xlabel(r'$Mutation\ rate$')
#     ax.set_xbound(lower=10**-8, upper=10**-6)
#     ax.set_ybound(lower=10**6, upper=10**8)
plt.colorbar(im, ax=ax)
# fig.colorbar(im, ax=axarr.ravel(), orientation='horizontal', shrink=0.85, pad=0.2)
#     plt.setp(cbar.get_ticks(), size=8)
ax.set_ylabel(r'$Cell\ number$')
plt.tight_layout()
plt.savefig('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/28Jun2018/linear_v2_aa_rate_np.svg')


# ### Crude incidence per 100000

# In[85]:


fig, axarr = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(12.5,12.5))
age_classes = [r'$0$', r'$1-4$', r'$5-9$', r'$10-14$', r'$15-19$', r'$20-24$', r'$25-29$', r'$30-34$', r'$35-39$', r'$40-44$', r'$45-49$', r'$50-54$', r'$55-59$', r'$60-64$', r'$65-69$', r'$70-74$', r'$75-79$', r'$80-84$', r'$85-89$']

fig.text(0.5, -0.01, r'$Age\ classes\ (yrs)$', ha='center', fontsize=15)
# fig.text(-0.05, 0.15, r'$\sigma=4$')
# fig.text(-0.05, 0.4, r'$\sigma=2$')
# fig.text(-0.05, 0.6, r'$\sigma=1$')
# fig.text(-0.05, 0.85, r'$\sigma=0.5$')
# fig.text(-0.05, 0.9, r'$\mu=0$')
fig.text(-0.05, 0.5, r'$Crude\ rate\ per\ 100000$', va='center', rotation='vertical', fontsize=15)

for a in axarr[-1]:
#     a.set_xlabel('Age class')
    a.set_xticks(numpy.arange(19))
    a.set_xticklabels(age_classes)
    plt.setp(a.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=6)

j = 0
for a in range(len(axarr)):
    k = 0
    ax = axarr[a]
    cr = crr[::-1][a*5:(a+1)*5][::-1] #Everything's been reversed because the original gvar has been declared in descending order
    for i in range(len(ax)):
        ax[i].plot(numpy.arange(n_class), cr[i])
        if j == 0:
            ax[i].set_title(nl[i])
        if k == 0:
            ax[i].set_ylabel(pl[j], rotation=0, labelpad=25)
        k += 1
#     if j == 0:
#         ax[-1].legend()
    j += 1

# fig.legend(p_labels, bbox_to_anchor=(0.6,0.02), ncol=5)
plt.tight_layout()
plt.savefig('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/28Jun2018/linear_v2_crude_rate_np.svg')


# ### Cumulative incidence

# In[88]:


fig, axarr = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(12.5,12.5))
age_classes = [r'$0$', r'$1-4$', r'$5-9$', r'$10-14$', r'$15-19$', r'$20-24$', r'$25-29$', r'$30-34$', r'$35-39$', r'$40-44$', r'$45-49$', r'$50-54$', r'$55-59$', r'$60-64$', r'$65-69$', r'$70-74$', r'$75-79$', r'$80-84$', r'$85-89$']

fig.text(0.5, -0.01, r'$Age\ classes\ (yrs)$', ha='center', fontsize=15)
# fig.text(-0.05, 0.15, r'$\sigma=4$')
# fig.text(-0.05, 0.4, r'$\sigma=2$')
# fig.text(-0.05, 0.6, r'$\sigma=1$')
# fig.text(-0.05, 0.85, r'$\sigma=0.5$')
# fig.text(-0.05, 0.9, r'$\mu=0$')
fig.text(-0.05, 0.5, r'$Cumulative\ rate$', va='center', rotation='vertical', fontsize=15)

for a in axarr[-1]:
#     a.set_xlabel('Age class')
    a.set_xticks(numpy.arange(19))
    a.set_xticklabels(age_classes)
    plt.setp(a.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=6)

j = 0
for a in range(len(axarr)):
    k = 0
    ax = axarr[a]
    cr = cc[::-1][a*5:(a+1)*5][::-1] #Everything's been reversed because the original gvar has been declared in descending order
    for i in range(len(ax)):
        ax[i].plot(numpy.arange(n_class), cr[i])
        if j == 0:
            ax[i].set_title(nl[i])
        if k == 0:
            ax[i].set_ylabel(pl[j], rotation=0, labelpad=25)
        k += 1
#     if j == 0:
#         ax[-1].legend()
    j += 1

# fig.legend(p_labels, bbox_to_anchor=(0.6,0.02), ncol=5)
plt.tight_layout()
plt.savefig('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/figures/28Jun2018/linear_v2_cumulative_rate_np.svg')


# ## Data export

# In[101]:


df_crr = pd.DataFrame(crr, index = [numpy.ravel([[x,x,x,x,x] for x in gvar]), numpy.array([gvar]*5).ravel()])
df_crr.index.names = ['Var(g)', 'Var(g)']
df_crr.columns.names = ['Age classes']
df_crr.columns = age_classes
df_crr.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/28Jun2018/linear_v2_crude_rate_np.xlsx')

df_cc = pd.DataFrame(cc, index = [numpy.ravel([[x,x,x,x,x] for x in gvar]), numpy.array([gvar]*5).ravel()])
df_cc.index.names = ['Var(g)', 'Var(g)']
df_cc.columns.names = ['Age classes']
df_cc.columns = age_classes
df_crr.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/28Jun2018/linear_v2_cumulative_rate_np.xlsx')

df_aar = pd.DataFrame(aa_rate.reshape(5,5), index=gvar, columns=gvar)
df_aar.index.names = ['Var(g)']
df_aar.columns.names = ['Var(g)']
df_crr.to_excel('/home/vibishan/PhD/Research/cancer_project/cancer_incidence_model/linear_model/V2/data/28Jun2018/linear_v2_aa_rate_np.xlsx')

