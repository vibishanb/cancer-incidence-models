#########################################################################
######GETTING CUMULATIVE RATES BY CANCER SITE FROM SEER 18 DATA##########
#########################################################################
"""
10th December, 2018:

SEER 18 data sets give the probability of developing cancer from a given cancer-free age up to an upper age limit. The lifetime risk for a given cancer site is then given by the sum of the probability of cancer from cancer-free age of 0 years, up to the last age class, which in this case is 95+ years.
In the code below, the 'Risk' at cancer-free age, 0 is grouped by cancer site, and the lifetime risk of cancer is given by the value of the cancer probability at age 95+, for a cancer-free age of 0 years."""


df1 = pd.read_csv('SEER18-cancer-developing-probability-by-site-part-1.csv', skiprows=[0, 1, 2, 3, 4])
cr_series_1 = df1.where(df1['Cancer Free Age']=='Age 0')['Risk'].groupby([df1['Cancer Site'], df1['Probability of Developing Cancer (%) by Age']]).sum()
cr_summary_1 = cr_series_1.groupby(by='Cancer Site').max()

df2 = pd.read_csv('SEER18-cancer-developing-probability-by-site-part-2.csv', skiprows=[0, 1, 2, 3, 4])
cr_series_2 = df2.where(df2['Cancer Free Age']=='Age 0')['Risk'].groupby([df2['Cancer Site'], df2['Probability of Developing Cancer (%) by Age']]).sum()
cr_summary_2 = cr_series_2.groupby(by='Cancer Site').max()

df3 = pd.read_csv('SEER18-cancer-developing-probability-by-site-part-3-female-sites-only.csv', skiprows=[0, 1, 2, 3, 4])
cr_series_3 = df3.where(df3['Cancer Free Age']=='Age 0')['Risk'].groupby([df3['Cancer Site'], df3['Probability of Developing Cancer (%) by Age']]).sum()
cr_summary_3 = cr_series_3.groupby(by='Cancer Site').max()

df4 = pd.read_csv('SEER18-cancer-developing-probability-by-site-part-4-male-sites-only.csv', skiprows=[0, 1, 2, 3, 4])
cr_series_4 = df4.where(df4['Cancer Free Age']=='Age 0')['Risk'].groupby([df4['Cancer Site'], df4['Probability of Developing Cancer (%) by Age']]).sum()
cr_summary_4 = cr_series_4.groupby(by='Cancer Site').max()

devprob_final = pd.DataFrame([cr_summary_1, cr_summary_2, cr_summary_3, cr_summary_4]).T
devprob_final.to_csv('/home/iiser/PhD/github-cancer-incidence-models/all-data/random-chance-model-data/cancer-development-probability-by-site-seer18.csv')
