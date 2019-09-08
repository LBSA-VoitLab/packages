Custom package to calculate enrichment of various cell populations based on gene expression data. The backend database is based on single cell markers db created by Xinxin Zhang et al. Some parts of this package have been borrowed from iaconogi/enrichCM package on github.

Sample usage
```
results = binf.cellMarker::marker_enrich(gene_list,filter_category_1,filter_category_1_list,filter_category_2,filter_category_2_list)
results = binf.cellMarker::marker_enrich(gene_list,"tissueType",c("Peripheral blood","Blood"),"cellType",c("Normal cell"))
```
Enrichment results could be viewed by
```
enrichment_results = results$enrichments
enrichment_genes = results$genes
```

In progress
1. This package currently works for human cell markers and will be updated with other species
