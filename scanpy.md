```
artemia=scanpy.read_h5ad("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_7/data_filt.h5ad", chunk_size=6000)
features_AR=pd.DataFrame({'features': artemia.raw.var_names})
artemia_2=ad.AnnData(artemia.raw.X,obs=artemia.obs,var=features_AR)
artemia_2.write_h5ad("/nfs/scistore18/vicosgrp/melkrewi/Project_snRNA_ovaries/22.cellblender/analysis_7/data_filt_clean.h5ad", as_dense=())
```
