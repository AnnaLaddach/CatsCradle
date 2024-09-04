# CatsCradle 

<img src="vignettes/CatsCradleLogo.png" alt="" width="200"/>

CatsCradle is an R package which provides tools for exploring duality
present in single cell datasets. It interfaces with analyses performed
with the Seurat Package, although many functions take simple lists,
dataframes or matrices as input.

In a typical scRNAseq dataset duality exists between cells and the
genes they express. Cells are typically assigned cell types based on
gene expression; however, programmes of gene expression can cut across
cell types. Furthermore, information on how genes are coexpressed can
be leveraged to infer functionality. Please see vignette
[CatsCradle](vignettes/CatsCradle.Rmd) to find out how you can apply
CatsCradle to explore this duality.

Even more relationships exist in high resolution spatial
transcriptomics datasets, between genes, cells and their position in
space. Here we use graph-based representations of cell position, and
provide functionality for extending these graphs to account for
longer-range interactions.  We also provide functions for identifying
and exploring neighbourhoods; in terms of gene expression and cell
type composition.  Finally, we provide tools for analysing
ligand-receptor interactions, taking spatial information into
account. Please see vignette
[CatsCradleSpatial](vignettes/CatsCradleSpatial.Rmd) for details on
how to use CatsCradle to analyse spatial data.

To dive straight into CatsCradle please see our quickstart vignette
[CatsCradleQuickStart](vignettes/CatsCradleQuickStart.Rmd)!

## Installation

```
# Install devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("AnnaLaddach/CatsCradle")
```


#### Developed by Anna Laddach and Michael Shapiro in the Pachnis lab at the Francis Crick Institute.
     

