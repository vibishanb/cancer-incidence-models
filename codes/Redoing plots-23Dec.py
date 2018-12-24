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

# Removing lines from S2-1 and S2-2
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