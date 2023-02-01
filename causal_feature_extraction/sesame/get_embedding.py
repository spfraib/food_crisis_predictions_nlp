# -*- coding: utf-8 -*-
import json
import os
import sys
import time
from optparse import OptionParser

from dynet import *
from evaluation import *
from raw_data import make_data_instance
from semafor_evaluation import convert_conll_to_frame_elements


optpr = OptionParser()
optpr.add_option("--output_dir", type='str', default='')
(options, args) = optpr.parse_args()

base_dir = "/scratch/ab7325/frame-models"
output_file_dir = "{}/logs/{}/".format(base_dir, options.output_dir)

USE_WV = True

sys.stderr.write("_____________________\n")
sys.stderr.write("COMMAND: {}\n".format(" ".join(sys.argv)))


post_train_lock_dicts()
lufrmmap, relatedlus = read_related_lus()
if USE_WV:
    pretrained_embeddings_map = get_wvec_map()
    PRETRAINED_DIM = len(pretrained_embeddings_map.values()[0])

lock_dicts()
UNKTOKEN = VOCDICT.getid(UNK)

from nltk.cluster import KMeansClusterer
import nltk

def add_embeddings(df, wvec_map, voc_dict):
    cause_embds = []
    effect_embds = []
    for r in df.iterrows():
    	if 'cause' in r[1] and 'effect' in r[1]:
	    cause = r[1]['cause'].split(" ")
	    effect = r[1]['effect'].split(" ")
	    cause_emb = [0.0]*PRETRAINED_DIM
	    for c in cause:
	    	try:
	    	    cause_emb = np.add(wvec_map[voc_dict.getid(c)], cause_emb)
	        except:
	    	    continue
	    cause_emb = np.divide(cause_emb, len(cause))
	    effect_emb = [0.0]*PRETRAINED_DIM
	    for e in effect:
	        try:
	    	    effect_emb = np.add(wvec_map[voc_dict.getid(c)], effect_emb)
	    	except:
	    	    continue
	    effect_emb = np.divide(effect_emb, len(effect))
	    cause_embds.append(cause_emb)
	    effect_embds.append(effect_emb)

    NUM_CLUSTERS=10
    kclusterer = KMeansClusterer(NUM_CLUSTERS, distance=nltk.cluster.util.cosine_distance, repeats=25, avoid_empty_clusters=True)
    assigned_cause_clusters = kclusterer.cluster(cause_embds, assign_clusters=True)
    assigned_effect_clusters = kclusterer.cluster(effect_embds, assign_clusters=True)
    return df.assign(cause_emb=cause_embds, effect_emb=effect_embds, cause_cluster=assigned_cause_clusters, effect_cluster=assigned_effect_clusters)


import pandas as pd
C = pd.read_csv('{}args-tuples.csv'.format(output_file_dir), delimiter='|', names=['cause','word','effect'], header=None)
D = C.dropna()
corrected = pd.DataFrame()
for r in D.iterrows():
    if 'cause' in r[1] and 'effect' in r[1] and 'word' in r[1]:
	if r[1]['effect'].startswith('by'):
            corrected = corrected.append({'cause': r[1]['effect'], 'effect': r[1]['cause'], 'word': r[1]['word']}, ignore_index=True)
    	else:
            corrected = corrected.append(r[1], ignore_index=True)
        
corrected.to_csv('pruned_args_tuples.csv')

X = add_embeddings(corrected, pretrained_embeddings_map, VOCDICT)

cause_F = pd.DataFrame()
effect_F = pd.DataFrame()
for r in X.iterrows():
    if 'famine' in r[1]['effect']:
        cause_F = cause_F.append(r[1])
    elif 'famine' in r[1]['cause']:
        effect_F = effect_F.append(r[1])
cause_F.to_csv('{}causes_of_famine.csv'.format(output_file_dir))
effect_F.to_csv('{}effects_of_famine.csv'.format(output_file_dir))

