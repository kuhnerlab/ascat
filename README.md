Allele-Specific Copy number Analysis of Tumors 
======

Description
--------
This contains a modified version of v2.1 of the ASCAT code.  ASCAT is described in [this](http://www.ncbi.nlm.nih.gov/pubmed/20837533) publication.  The main modifications are described in [this](https://www.nature.com/articles/s41467-017-02621-x) publication.

The modifications include:
* Replacement of ASCAT's segmentation algorithm with one from 'copynumber' (described in [this](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-591) publication, and available from https://git.bioconductor.org/packages/copynumber and https://bioconductor.org/packages/release/bioc/html/copynumber.html)
* Allow purities of 1.0, if no other maximum was found, according to a user setting.

Further changes were made and will be described in a forthcoming paper from Smith, Yamato, and Kuhner.  They include:
* Allow a constrained search of 'diploid-like' or 'tetraploid-like' ploidies (psi values), according to a user setting.

