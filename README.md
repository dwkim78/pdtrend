# PDT

<div align="center"><img src="./pdt/datasets/images/PDT_logo.png"></div>

PDT (<b>P</b>hotometric <b>D</b>e<b>T</b>rending Algorithm  aims to remove systematic trends in the light curves. For details about the algorithm, see [Kim et al. 2009](http://adsabs.harvard.edu/abs/2009MNRAS.397..558K).


The latest PDT uses [Birch](https://en.wikipedia.org/wiki/BIRCH) to find highly-correlated light curves rather than [Hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering) that [Kim et al. 2009](http://adsabs.harvard.edu/abs/2009MNRAS.397..558K) originally used. This is mainly because 1) Birch is scalable (i.e. applicable to large dataset), and 2) Birch does not need to set the number of clusters.


Note that the input light curves <b>must</b> have the same number of data points. PDT then assumes that each light curve is synced in time. PDT is not designed to deal with missing data points or desynced data points. Also note that the light curves must be cleaned beforehand (e.g. highly-fluctuated data points, etc). Nevertheless, such functionality (e.g. dealing with desynced data, pre-processing, etc.) might be implemented in the future.  

## Index
1. [Dependency](#1-dependency)
2. [Installation](#2-installation)
3. [Test the Installation](#3-test)
4. [How to Use PDT](#5-how-to-use-pdt)

- [ChangeLog](#changelog)
- [Citation](#citation)


## Dependency

[Python 2.7+](https://www.python.org/) 

 * Not tested with Python 3.0+

[Numpy 1.10+](http://www.numpy.org/)

[Scipy 0.17+](http://www.scipy.org/)
 
[Scikit-learn 0.17+](http://scikit-learn.org/)

[Matplotlib 1.5+](http://matplotlib.org/)


These libraries will be automatically installed if your machine does not have them installed. If you encounter errors during the installation of these dependencies, try to install them individually. Your machine may not have other required libraries by these dependencies.


## Installation


## Test


## How to Use PDT


## ChangeLog

### v.0.1
- pre-alpha version

### v.0.0.0
- create the GitHub repository 

## Citation

If you use PDT in publication, we would appreciate citations to the paper, 
[Kim et al. 2009](http://adsabs.harvard.edu/abs/2009MNRAS.397..558K).


## Contact
Dae-Won Kim, email: dwkim78 at gmail.com

Webpage: https://sites.google.com/site/dwkim78/

#### Keywords

astronomy - light curves - trend removal - detrend - machine learning

