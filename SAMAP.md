```
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder, transfer_annotations,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import scanpy
import scanpy.external as sce
import scanpy as sc

fn1 = '/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_7/data_filt_clean.h5ad'
fn2='/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_6/s_fca_biohub_ovary_10x_clean.h5ad'
sm = SAMAP(
        filenames,
        f_maps = '/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/8.integrate_with_drosophila/5.integrate/maps_v4/',
        save_processed=True, #if False, do not save the processed results to `*_pr.h5ad`
        resolutions = {'AR':0.7,'DM':0.7}
    )
sm.run(pairwise=True)
keys = {'AR':'seurat_clusters','DM':'annotation'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
sankey_plot(MappingTable, align_thr=0, species_order = ['AR','DM'])
chord_plot(MappingTable, align_thr=0.2)
fig,ax = plt.subplots();
fig.set_size_inches((10,10))
sm.scatter(axes=ax)
gpf = GenePairFinder(sm,keys=keys)
gene_pairs = gpf.find_all(align_thr=0.05)
table=transfer_annotations(sm,reference_id='DM',keys='annotation')

import matplotlib.pyplot as plt
#ax.set_figure_params(figuresize=(8,8))
fig,ax = plt.subplots();
fig.set_size_inches((15,15)) 
sm.samap.scatter(c='annotation_transfer',axes=ax, cmap=make_random_cmap(24)) #
#plt.savefig('/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_6/annotation_transfer.png',bbox_inches='tight',facecolor='white', edgecolor='white')
