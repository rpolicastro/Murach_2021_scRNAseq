
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import scvelo as scv
import cellrank as cr

## Velocity
## ----------

## Load the data.

seu_spliced = Path("results/py_objects/Pax7_DTA_4day.h5ad")
scv_data = scv.read(seu_spliced)

scv_data.obs = scv_data.obs.astype('category')

## Preprocess the data.

scv.pp.filter_and_normalize(scv_data, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(scv_data, n_pcs=30, n_neighbors=30)

## Calculate RNA velocities.

scv.tl.recover_dynamics(scv_data)
scv.tl.velocity(scv_data, mode='dynamical')
scv.tl.velocity_graph(scv_data)

## Plot RNA velocities.

scv.pl.velocity_embedding_stream(scv_data, basis='umap', color='integrated_snn_res.0.4')

scv.pl.velocity_embedding(
  scv_data, arrow_length=3, arrow_size=2, dpi=120,
  basis = 'umap'
)

## PAGA.

scv_data.uns['neighbors']['distances'] = scv_data.obsp['distances']
scv_data.uns['neighbors']['connectivities'] = scv_data.obsp['connectivities']

scv.tl.paga(scv_data, groups = 'integrated_snn_res.0.4')
scv.pl.paga(scv_data, color = 'integrated_snn_res.0.4', basis = 'umap')

## Latent time.

scv.tl.latent_time(scv_data)
scv.pl.scatter(scv_data, color='latent_time', color_map='gnuplot', size=80)

## Run cellrank.

scv.tl.velocity_graph(scv_data, mode_neighbors='connectivities', compute_uncertainties=True)

cr.tl.final_states(
  scv_data, cluster_key='integrated_snn_res.0.4', weight_connectivities=0.2,
  use_velocity_uncertainty=True
)

cr.tl.root_states(scv_data, cluster_key='integrated_snn_res.0.6')
