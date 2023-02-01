world['corr'] = [np.nan]*len(world)
for c in mapped.keys():
    if c == 'Democratic Republic of the Congo':
        cx = 'Dem. Rep. Congo'
    elif c == 'Central African Republic':
        cx = 'Central African Rep.'
    elif c == "South Sudan":
        cx = "S. Sudan"
    else:
        cx = c
    
#     if sum(world['name']==cx) == 0:
#         c_df = world.query('name =="{}"'.format(c))
#         c_df['corr'] = country_corr[c]
#         africa = pd.concat([africa, c_df])
#     else:
    idx = world.index[world['name']==cx]
    world.at[idx, 'corr'] = country_corr[c]


import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1)
world.plot(column='corr', ax=ax, legend=True, missing_kwds={'color': 'lightgrey'})

import pandas as pd
d = pd.read_csv('all_causes_country_timeline.csv')
d.index = pd.to_datetime(d.date)
d.info()

from sklearn.preprocessing import StandardScaler
mat_scaled = StandardScaler()
non_columns = ['Unnamed: 0', 'Unnamed: 0: news', 'Unnamed: 0: literature', 'date', 'country']
select_columns = []
for x in d.columns:
    if x not in non_columns:
        select_columns.append(x)
mat_sc = mat_scaled.fit_transform(d[select_columns])

fews = pd.read_csv('famine-country-province-district-years-CS.csv')

fews['CS'] = fews['CS'].replace(99.0, 0.0).replace(66.0, 0.0).replace(88.0, 0.0)
fews = fews.groupby(['country', 'year', 'month'])['CS'].max().reset_index()
fews['day'] = [1]*len(fews)
fews.index = pd.to_datetime(fews[['year', 'month', 'day']])
fews = fews.groupby('country').resample('MS').mean().interpolate().reset_index()
fews.index = fews['level_1']
# fews = fews.reindex(['2010-01-01:2017-10-01'])

import numpy as np
print (np.unique(fews['country'])) #no abyei, ilyemi triangle
print (np.unique(d['country']))
mapped = {'Afghanistan': 'AFGH', 'Angola': 'ANGOL', 'Burundi': 'BURUN', 'Cameroon': 'CAMER',
         'Central African Republic': 'CAFR', 'Chad': 'CHAD', 'El Salvador': 'ELSAL', 
         'Ethiopia': 'ETHPA', 'Mauritania': 'MAURTN', 'Guatemala': 'GUAT', 'Guinea': 'GUREP',
         'Haiti': 'HAIT', 'Honduras': 'HON', 'Kenya': 'KENYA', 'Liberia': 'LIBER', 'Madagascar': 'MALAG',
         'Malawi': 'MALAW', 'Mali': 'MALI', 'Mozambique': 'MOZAM', 'Niger': 'NIGER', 'Nigeria': 'NIGEA',
         'Rwanda': 'RWANDA', 'Senegal': 'SENEG', 'Sierra Leone': 'SILEN', 'Somalia': 'SOMAL', 'South Sudan': 'SOUSUD',
         'Sudan': 'SUDAN', 'Tajikistan': 'TADZK', 'Djibouti': 'TAI', 'Tanzania': 'TANZA', 'Uganda': 'UGANDA',
         'Burkina Faso': 'UPVOLA', 'Yemen': 'YEMAR', 'Democratic Republic of the Congo': 'ZAIRE', 'Zambia': 'ZAMBIA',
         'Zimbabwe': 'ZIMBAB'}

corr_df = pd.DataFrame()
k =0 
country_corr = {}
country_var = {}
for c in np.unique(fews['country']):
    if c in mapped:
        corr_max = 0
        corr_v = None
        for v in select_columns:
            corr = fews[fews['country']==c]['CS'].corr(d[d['country'] == mapped[c]][v])
            corr_df = corr_df.append(pd.Series([c, v, corr]), ignore_index=True)
            k += 1
            if corr_max < abs(corr):
                corr_max = abs(corr)
                corr_v = v
        country_corr[c] = corr_max
        country_var[c] = corr_v
        print (c, corr_v, corr_max)

corr_df = corr_df.rename(columns={0: 'country', 1: 'variable', 2: 'correlation'})
corr_df.to_csv('var_country_corr.csv')
