# Redoing plots-23Dec
# Main plots from both models; plotted with viridis, and the same log scales.

pylab

df_cc=  pd.read_excel('/home/iiser/PhD/github-cancer-incidence-models/all-data/v1-data/22Oct2018/linear-v1-multipop-cumulative-count-n.xlsx', index_col=0)

df_crc = pd.read_excel('/home/iiser/PhD/github-cancer-incidence-models/all-data/v1-data/22Oct2018/linear-v1-multipop-cancer-count-n.xlsx', index_col=0)  

df_surv = pd.DataFrame(numpy.array([[Npop]*100]*11), index=df_cc.index)

df_crr = (df_crc/(df_surv-df_cc))*Npop

df_crr.T.plot(legend=False); legend(log10(narr).round(decimals=2), title='log(n)'); ylabel('Crude incidence'); xlabel('Time (years)'); show()

(df_cc*100/Npop).T.plot(legend=False); legend(log10(narr).round(decimals=2), title='log(n)'); ylabel('Cumulative incidence (%)'); xlabel('Time (years)'); show()

half_age_narr = array([less_equal(df_cc[i], df_cc[99]/2) for i in range(100)]).sum(axis=0)
df_cc['Half max age'] = half_age_narr
df_cc.reset_index().plot(x='Cell number', y='Half max age', logx=True, legend=False, marker='o'); ylabel(r'Age at $\frac{I_{max}}{2}$'); show()

df_cc.reset_index().plot(x='Cell number', y=99, logx=True, legend=False, marker='o'); ylabel(r'$I_{max}$'); show()


# Additional plots for sensitivity; instead of g vs p/n, g vs cancer-time for both distributions
df_nrand = pd.read_csv('/home/iiser/PhD/github-cancer-incidence-models/all-data/sensitivity-data/gdist/linear-v2-cancer-time-gdist-gumbel-nrand-16Nov.csv')
df_prand = pd.read_csv('/home/iiser/PhD/github-cancer-incidence-models/all-data/sensitivity-data/gdist/linear-v2-cancer-time-gdist-gumbel-prand-16Nov.csv')

sns.jointplot(x='g(k-1)', y='Time to cancer', data=df_nrand, kind='reg', stat_func=stats.pearsonr)
tight_layout()
savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/sensitivity-figures/gdist/linear-v2-g-vs-cancer-time-gumbel-nrand-23Dec.svg')
show()

sns.jointplot(x='g(k-1)', y='Time to cancer', data=df_prand, kind='reg', stat_func=stats.pearsonr)
tight_layout()
savefig('/home/iiser/PhD/github-cancer-incidence-models/all-figures/sensitivity-figures/gdist/linear-v2-g-vs-cancer-time-gumbel-prand-23Dec.svg')
show()

# Removing regression lines from S2-1 and S2-2
df_nrand=  pd.read_csv('/home/iiser/PhD/github-cancer-incidence-models/all-data/sensitivity-data/linear-v1-cancer-time-k&n-rand.csv', index_col=0)
df_prand = pd.read_csv('/home/iiser/PhD/github-cancer-incidence-models/all-data/sensitivity-data/linear-v1-cancer-time-k&p-rand.csv', index_col=0)

sns.jointplot(x='log(n)', y='Time to cancer', data=df_nrand)
sns.jointplot(x='log(p)', y='Time to cancer', data=df_prand)

df_nrand[r'$\Delta_{g}$'] = (df_nrand['g(k-1)'] - 0.007)/df_nrand['Threshold']
df_prand[r'$\Delta_{g}$'] = (df_prand['g(k-1)'] - 0.007)/df_prand['Threshold']

sns.jointplot(x=r'$\Delta_{g}$', y='Time to cancer', data=df_nrand.where(df_nrand['Threshold']==2), kind='reg', stat_func=stats.pearsonr)
show()
sns.jointplot(x=r'$\Delta_{g}$', y='Time to cancer', data=df_nrand.where(df_nrand['Threshold']==4), kind='reg', stat_func=stats.pearsonr)
show()
sns.jointplot(x=r'$\Delta_{g}$', y='Time to cancer', data=df_nrand.where(df_nrand['Threshold']==6), kind='reg', stat_func=stats.pearsonr)
show()

sns.jointplot(x=r'$\Delta_{g}$', y='Time to cancer', data=df_prand.where(df_prand['Threshold']==2), kind='reg', stat_func=stats.pearsonr)
show()
sns.jointplot(x=r'$\Delta_{g}$', y='Time to cancer', data=df_prand.where(df_prand['Threshold']==5), kind='reg', stat_func=stats.pearsonr)
show()
sns.jointplot(x=r'$\Delta_{g}$', y='Time to cancer', data=df_prand.where(df_prand['Threshold']==8), kind='reg', stat_func=stats.pearsonr)
show()

# Binned plots of crude and cumulative incidence in Figures 3A and E, and 5A and E
ncount1 = pd.read_excel('/home/iiser/PhD/github-cancer-incidence-models/all-data/v1-data/22Oct2018/linear-v1-multipop-cancer-count-n.xlsx', index_col=0)
pcount1 = pd.read_excel('/home/iiser/PhD/github-cancer-incidence-models/all-data/v1-data/22Oct2018/linear-v1-multipop-cancer-count-p.xlsx', index_col=0)
ncount2 = pd.read_excel('/home/iiser/PhD/github-cancer-incidence-models/all-data/v2-data/22Oct2018/linear-v2-multipop-cancer-count-n.xlsx', index_col=0)
pcount2 = pd.read_excel('/home/iiser/PhD/github-cancer-incidence-models/all-data/v2-data/22Oct2018/linear-v2-multipop-cancer-count-p.xlsx', index_col=0)

n_class = 20
n1_binned = pd.DataFrame(array([ncount1.iloc[:, 5*i:(5*i)+5].sum(axis=1) for i in range(n_class)]).T, index=ncount1.index)
p1_binned = pd.DataFrame(array([pcount1.iloc[:, 5*i:(5*i)+5].sum(axis=1) for i in range(n_class)]).T, index=pcount1.index)
n2_binned = pd.DataFrame(array([ncount2.iloc[:, 5*i:(5*i)+5].sum(axis=1) for i in range(n_class)]).T, index=ncount2.index)
p2_binned = pd.DataFrame(array([pcount2.iloc[:, 5*i:(5*i)+5].sum(axis=1) for i in range(n_class)]).T, index=pcount2.index)

n1_cumul = n1_binned.cumsum(axis=1)
p1_cumul = p1_binned.cumsum(axis=1)
n2_cumul = n2_binned.cumsum(axis=1)
p2_cumul = p2_binned.cumsum(axis=1)

n1_crude = n1_binned.copy()
p1_crude = p1_binned.copy()
n2_crude = n2_binned.copy()
p2_crude = p2_binned.copy()

"""The (k-1) in the sum term in the denominator is important; cancer cases up to previous time point go into calculating the surviving fraction in the current time point; sum up to the kth index would also include the cancer cases occurring in the current time point, which gives anomalous age-specific rates > 1."""
n1_crude.iloc[:, 1:] = pd.DataFrame(array([n1_binned.iloc[:,k]/(Npop-n1_binned.iloc[:, :(k-1)].sum(axis=1)) for k in range(1,n_class)]).T, index=n1_binned.index)
p1_crude.iloc[:, 1:] = pd.DataFrame(array([p1_binned.iloc[:,k]/(Npop-p1_binned.iloc[:, :(k-1)].sum(axis=1)) for k in range(1,n_class)]).T, index=p1_binned.index)
n2_crude.iloc[:, 1:] = pd.DataFrame(array([n2_binned.iloc[:,k]/(Npop-n2_binned.iloc[:, :(k-1)].sum(axis=1)) for k in range(1,n_class)]).T, index=n2_binned.index)
p2_crude.iloc[:, 1:] = pd.DataFrame(array([p2_binned.iloc[:,k]/(Npop-p2_binned.iloc[:, :(k-1)].sum(axis=1)) for k in range(1,n_class)]).T, index=p2_binned.index)

cols = array([str(5*i)+'-'+str((5*i)+4) for i in range(n_class)])
n1_crude.columns = cols
p1_crude.columns = cols
n2_crude.columns = cols
p2_crude.columns = cols

n1_cumul.columns = cols
p1_cumul.columns = cols
n2_cumul.columns = cols
p2_cumul.columns = cols

n1_cumul.columns.names = ['Age class']
p1_cumul.columns.names = ['Age class']
n2_cumul.columns.names = ['Age class']
p2_cumul.columns.names = ['Age class']

n1_crude.columns.names = ['Age class']
p1_crude.columns.names = ['Age class']
n2_crude.columns.names = ['Age class']
p2_crude.columns.names = ['Age class']

nl = array(log10(n1_crude.index)).round(decimals=2)
pl = array(log10(p1_crude.index)).round(decimals=2)

sns.set_palette(sns.color_palette("viridis", 13))
(n2_crude*100).T.plot(legend=False); legend(nl, title='log(n)'); xlabel('Age class', fontsize=12); ylabel('Age-specific incidence (%)', fontsize=15); xticks(range(n_class), cols, rotation=45, fontsize=10)
(p2_crude*100).T.plot(legend=False); legend(pl, title='log(p)'); xlabel('Age class', fontsize=12); ylabel('Age-specific incidence (%)', fontsize=15); xticks(range(n_class), cols, rotation=45, fontsize=10)
(n2_cumul*100/Npop).T.plot(legend=False); legend(nl, title='log(n)'); xlabel('Age class', fontsize=12); ylabel('Cumulative incidence (%)', fontsize=15); xticks(range(n_class), cols, rotation=45, fontsize=10)
(p2_cumul*100/Npop).T.plot(legend=False); legend(pl, title='log(p)'); xlabel('Age class', fontsize=12); ylabel('Cumulative incidence (%)', fontsize=15); xticks(range(n_class), cols, rotation=45, fontsize=10)