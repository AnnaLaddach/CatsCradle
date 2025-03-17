# CatsCradle 

<img src="CatsCradleLogoTransparent.png" alt="" width="200"/>

CHANGES IN VERSION 0.99.10
--------------------------

* We have now broken out our code into files pertaining to different
  areas of functionality.
* We have improved functionality for
  handling SingleCellExperiment and SpatialExperiment objects.
* CatsCradle is currently under review at Bioconductor.  Watch this
  space.

CHANGES IN VERSION 1.1.1
--------------------------

* Default behaviour of nbhdsAsEdgesToNbhdsAsList changed to exclude 
  self. This fixes a bug leading to inflated values of MoransI.
  
CHANGES IN VERSION 1.1.2
--------------------------

* Function getSubsetComponents added. This is designed to dectect the 
  components of a gene subset in the case where median complement distance
  detects clustering.
* Function edgeLengthsAndCellTypePairs rewritten for faster runtime.
