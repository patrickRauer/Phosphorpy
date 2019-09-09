# Phosphorpy
**Phosphorpy** is python package to mine large photometric sky surveys. 
It is designed to allow to do common and regular task, which are done if
large photometric data sets are used, in just few lines. 
The aim is to provide a simple interface to make such data sets and 
speially data set combination more accessable to the community.

These task are for example, cross-match different catalogs (photometric and 
astrometric), make basic selection, plotting results, simple image and/or
light curves downloads.

All accessable data are publicly available and for the references, see the
[survey config file](https://gitlab.sron.nl/asg/jonker/Phosphorpy/blob/master/Phosphorpy/local/survey.conf).
This file includes the references to the corresponding papers.

Details are on the [wiki page](https://gitlab.sron.nl/patrickr/Phosphorpy/wikis/home)

## Installation
**Phosphorpy** is available via [Pypi](https://pypi.org/) and can be installt with pip
```
pip install Phosphorpy
```

## Requirements
**Phosphorpy** uses planty of other packages to perform the queries, downloads
and do the calculations.

* astropy
* astroquery
* matplotlib
* numpy
* numba
* pandas
* seaborn
* scikit-learn
* scipy
* urllib
* armapy
* extinction

## Example
See [example](https://gitlab.sron.nl/patrickr/Phosphorpy/wikis/Example) for
a almost complete example.
