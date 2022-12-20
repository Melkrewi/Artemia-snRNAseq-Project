```
dnucl_2=scanpy.read_h5ad("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_6/s_fca_biohub_ovary_10x.h5ad", chunk_size=6000)
features=pd.DataFrame({'features': dnucl_2.raw.var_names})
dnucl_4=ad.AnnData(dnucl_2.raw.X,obs=dnucl_2.obs[['n_counts','n_genes','annotation']],var=features)
dnucl_4.write_h5ad("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_6/s_fca_biohub_ovary_10x_clean.h5ad", as_dense=())

artemia=scanpy.read_h5ad("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_7/data_filt.h5ad", chunk_size=6000)
features_AR=pd.DataFrame({'features': artemia.raw.var_names})
artemia_2=ad.AnnData(artemia.raw.X,obs=artemia.obs,var=features_AR)
artemia_2.write_h5ad("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_7/data_filt_clean.h5ad", as_dense=())
```
