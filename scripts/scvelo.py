
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import scvelo as scv
import os
import pickle

## Variables.

clusters = {
  'KY_Mononuclear' : 'integrated_snn_res.0.4',
  'tdT_Parental' : 'integrated_snn_res.0.4',
  'Pax7_DTA_4day' : 'integrated_snn_res.0.5',
  'Pax7_tdT_4day' : 'integrated_snn_res.0.5'
}

h5ad_files = {
  'KY_Mononuclear' : 'results/py_objects/KY_Mononuclear.h5ad',
  'tdT_Parental' : 'results/py_objects/tdT_Parental.h5ad',
  'Pax7_DTA_4day' : 'results/py_objects/Pax7_DTA_4day.h5ad',
  'Pax7_tdT_4day' : 'results/py_objects/Pax7_tdT_4day.h5ad'
}

####################################
## RNA Velocity and PAGA Analysis ##
####################################

## Load the samples.

samples = {x:scv.read(y) for x,y in h5ad_files.items()}

## Change the metadata to categorical.

for key in samples.keys():
    samples[key].obs = samples[key].obs.astype('category')

## Preprocess the data.

for key in samples.keys():
    scv.pp.filter_and_normalize(samples[key], min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(samples[key], n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

for key in samples.keys():
    scv.tl.velocity(samples[key], mode_neighbors = 'connectivities')
    scv.tl.velocity_graph(samples[key])

## Plot RNA velocity streams.

outdir = 'results/trajectory/velocity/velocity_streams'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding_stream(
      value, basis = 'umap', color = clusters[key],
      save = '{}.png'.format(key), title = key, show = False,
      figsize = (10, 10), size = 50, dpi = 300
    )

## Plot RNA velocity arrows.

outdir = 'results/trajectory/velocity/velocity_arrows'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_embedding(
      value, arrow_length = 3, arrow_size = 2, dpi = 300,
      basis ='umap', color = clusters[key],
      figsize = (10, 10), size = 50, show = False,
      save = '{}.png'.format(key), title = key
    )

## Plot velocity speed and coherence.

outdir = 'results/trajectory/velocity/velocity_speed_coherence'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

metrics = ['velocity_length', 'velocity_confidence']
for key,value in samples.items():
    scv.tl.velocity_confidence(value)
    scv.pl.scatter(
      value, c = metrics, cmap = 'coolwarm', perc= [5, 95],
      size = 50, show = False, dpi = 300, figsize = (10, 10),
      save = '{}.png'.format(key)
    )

## Plot cell connections.

outdir = 'results/trajectory/velocity/velocity_connections'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.velocity_graph(
      value, threshold = .2, size = 50, show = False, dpi = 300,
      figsize = (10, 10), color = clusters[key],
      save = '{}.png'.format(key), title = key
    )

## Plot velocity pseudotime.

outdir = 'results/trajectory/velocity/velocity_pseudotime'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.tl.velocity_pseudotime(value)
    scv.pl.scatter(
      value, color = 'velocity_pseudotime', cmap = 'gnuplot', dpi = 300,
      show = False, figsize = (10, 10), title = key, size = 50,
      save = '{}.png'.format(key)
    )

## PAGA.

outdir = 'results/trajectory/velocity/velocity_paga'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key in samples.keys():
    samples[key].uns['neighbors']['distances'] = samples[key].obsp['distances']
    samples[key].uns['neighbors']['connectivities'] = samples[key].obsp['connectivities']

for key,value in samples.items():
    scv.tl.paga(value, groups = clusters[key])
    scv.pl.paga(
      value, basis = 'umap', color = clusters[key],
      dpi = 300, show = False, figsize = (10, 10), title = key, size = 50,
      save = '{}.png'.format(key)
    )

## Export UMAP plots.

outdir = 'results/trajectory/velocity/velocity_umap'
if not os.path.exists(outdir):
    os.makedirs(outdir)

scv.settings.figdir = outdir

for key,value in samples.items():
    scv.pl.umap(
      value, color = clusters[key], show = False, figsize = (10, 10), title = key,
      size = 50, save = '{}.png'.format(key)
    )
    scv.pl.umap(
      value, show = False, figsize = (10, 10), title = key,
      size = 50, save = '{}_nocolor.png'.format(key)
    )

## Get important genes.

outdir = 'results/trajectory/velocity/velocity_genes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

for key,value in samples.items():
    scv.tl.rank_velocity_genes(value, groupby = clusters[key], min_corr = .3)

for key,value in samples.items():
    df = scv.DataFrame(value.uns['rank_velocity_genes']['names'])
    df.to_csv("{}/{}.tsv".format(outdir, key), sep = '\t', header = True, index = False)

for key,value in samples.items():
    

## Save the velocities.

with open('results/py_objects/velocities.pickle', 'wb') as handle:
    pickle.dump(samples, handle)

## Load the velocities if required.

with open('results/py_objects/velocities.pickle', 'rb') as handle:
    samples = pickle.load(handle)

