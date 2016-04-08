# PDT

<div align="center"><img src="./pdtrend/datasets/images/PDT_logo.png"></div>

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

The easiest way to install the PDT package is:

```python
pip install pdtrend
```

Or,

```python
pip install git+https://github.com/dwkim78/pdtrend
```

If you do not want to install/upgrade the dependencies, execute the above commend with the ```--no-deps``` option. PDT possibly works with older version of Python and other libraries. 


Alternatively, you can download the PDT package from the Git repository as:

```python
git clone https://github.com/dwkim78/pdtrend

cd pdtrend
python setup.py install
```

You can edit ```setup.py```, if you do not want to update your own Python libraries (i.e. edit the ```install_requires``` variable).


## Test

To check if PDT is correctly installed, type following commands in the Python console.

```python
from pdtrend import test

test()
```

The command will print messages like:
```
yyyy-mm-dd hh:mm:ss,sss INFO - Loading the light curve set.
yyyy-mm-dd hh:mm:ss,sss INFO - 	The number of light curves is 57.
yyyy-mm-dd hh:mm:ss,sss INFO - Initializing pdtrend.
yyyy-mm-dd hh:mm:ss,sss INFO - Calculating the distance matrix.
yyyy-mm-dd hh:mm:ss,sss INFO - Searching for clusters using Birch.
yyyy-mm-dd hh:mm:ss,sss INFO - Filtering the clusters.
yyyy-mm-dd hh:mm:ss,sss INFO - Building master trends.
yyyy-mm-dd hh:mm:ss,sss INFO - Detrending one light curves using the master trends.
yyyy-mm-dd hh:mm:ss,sss INFO - Ploting results.
yyyy-mm-dd hh:mm:ss,sss INFO - Done.
```

This command reads the sample data set cosisting of 57 Pan-STARRS light curves (Python pickled and bzipped), run the clustering algorithm (i.e. Birch) to find clusters, construct master trends using the clusters, and detrend one sample light curve. In addition, it generates two images under the "./output" folder.

<div align="center"><img src="./pdtrend/datasets/images/master_trends.png" width="100%"œ><br/>[ Master Trends ]</div>

The above image shows the master trend constructed by the clustering algorithm. In this example data set, PDT found one master trend. For details about what is a master trend, see [Kim et al. 2009](http://adsabs.harvard.edu/abs/2009MNRAS.397..558K).

The following image is an example light curve before and after the detrending. Note that when PDT detrends a light curves, it minimized RMS of residuals while constraining weights for each master trend to be positive. The positive constraint is mandatory to avoid undesirable RMS minimization. For instance, if the weights are negative while the master trends are monotonically <b>increasing</b>, RMS minimization can reduce monotonically <b>decreasing</b> signals in light curves, which is unwanted. 

<div align="center"><img src="./pdtrend/datasets/images/detrended.png" width="100%"œ>[ An example of a detrended light curve ]</div>


## How to Use PDT

Using PDT is relatively simple because PDT assumes that light curves are synced. Nevertheless, note that PDT requires enough number of light curves to find clusters and master trends. We recommend to use PDT with at least 50 light curves.

The following pseudo code shows how to use PDT.

```
# Import PDT.
from pdtrend import PDTrend

# Read light curves.
lcs = ...

# Create PDT instance.
pdt = PDTrend(lcs)

# Run PDT to find clusters and then master trends.
pdt.run()

# Detrend each light curve.
for lc in lcs:
    detrended = pdt.detrend(lc)
```

In order to use PDT, light curve set must be read beforehand (e.g. the line  ```lcs = ...```). The ```lcs``` must consists of N rows and each row must contain M columns. N is the number of light curves and M is the number of data points. ```lcs``` could be either Python list or numpy.ndarry. For example:

```
lcs = [
        [1, 2, 3, 4, 5]
        [5, 4, 3, 2, 1]
        [3, 3, 3, 3, 3]
      ]
```

is the light curve set consisting of three light curves, each of which contains 5 data points.

When creating the PDT instance, you can set additional two options as:

| Option | Description |
|---:|:---|
| n_min_member | The minimum number of members in each cluster. If a cluster has fewer members, PDT discards the cluster. Default is 5. If you have a lot of light curves (e.g. several hundreds), you can increase this number to 10, 20, 30 or so. |
| dist_cut | The distance matrix that PDT uses is (1 - correlation matrix) / 2. Thus, if a cluster found by Birch consists of light curves of random Gaussian noise (i.e. no clear variability), it is likely that the median distance between the light curves is close to 0.5. Thus we can remove clusters whose median distance is larger than 0.5. Nevertheless, the default value is set to 0.45 in order to discard less-correlated clusters as well. If you increase this value (e.g. to 0.6 or so), PDT will construct master trends consisting of non-varying light curves. |


## ChangeLog

### v.0.2
- release alpha version

### v.0.1
- release pre-alpha version

### v.0.0.0
- create the GitHub repository 

## Citation

If you use PDT in publication, we would appreciate citations to the paper, 
[Kim et al. 2009](http://adsabs.harvard.edu/abs/2009MNRAS.397..558K) and this GitHub repository as well.


## Contact
Dae-Won Kim, email: dwkim78 at gmail.com

Webpage: https://sites.google.com/site/dwkim78/

#### Keywords

astronomy - light curves - trend removal - detrend - machine learning

