import pandas as pd
import numpy as np
from datetime import datetime

countries = ['AFGH', 'ANGOL', 'BURUN', 'CAFR', 'CAMER', 'CHAD', 'CONGO', 'ELSAL', 'ETHPA', 'GUAT', 'GUREP', 'HAIT', 'HON', 'KENYA', 'LIBER', 'MALAG', 'MALAW', 
    'MALI', 'MAURTN', 'MOZAM', 'NIGEA', 'NIGER', 'RWANDA', 'SENEG', 'SILEN', 'SOMAL', 'SOUSUD', 'SUDAN', 'TADZK', 'TAI', 'TANZA', 'UGANDA', 'UPVOLA', 'YEMAR',
    'ZAIRE', 'ZAMBIA', 'ZIMBAB']
full_index = {'date': [], 'country':[]}
for c in countries:
    dates = pd.date_range(start='1/1/2008', end='1/1/2018', freq='MS')
    full_index['date'].extend(dates)
    full_index['country'].extend([c]*len(dates))
all_vars = pd.DataFrame.from_dict(full_index)

phrases_file = '/scratch/ab7325/frame-models/logs/literature_causes.csv'
causes_df = pd.read_csv(phrases_file)
cause_phrases = np.unique(causes_df['cause'])
cause_words = set([])
for p in cause_phrases:
    cause_words.add(p)

for word in cause_words:
	wc_df = pd.read_csv('/scratch/ab7325/frame-models/logs/literature_causes/{}.csv'.format(word))
	wc_df['date'] = wc_df['date'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d'))
	all_vars = pd.merge(all_vars, wc_df[['date', 'country', word]], how='left', on=['date', 'country'])

all_vars.fillna(0, inplace=True)
all_vars.to_csv('/scratch/ab7325/frame-models/logs/literature_causes_word_country.csv')
