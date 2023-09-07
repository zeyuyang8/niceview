# Request

## Requirements for whole-slide image visualization

- Image
  - Original image [fix zoom range]
  - Original image overlayed with cell masks
    - Cell masks are colored by gene expression level
    - Cell masks are colored by cell type
  - Original image overlayed with spot given centroid coordinates and diameters. The spot size is fixed.
    - Spots are colored by gene expression level
    - Spots are transparent (only show the boundary)
    - Spots are colored by spot type (maybe)
  - Heatmap
    - Heatmap is colored by gene expression level (pathway enrichment score)
- Gene expression
- Callbacks
  - Draw a region and save it to a file. There should be a function to convert the GIS coordinate to the WSI pixels.

[high performance](https://numba.pydata.org/)
[leaflet demo for map analysis](https://www.youtube.com/watch?v=6-KTrjtfS9U)
