# Phosphorpy
Phosphorpy is python package to mine large photometric sky surveys. It (will) provide a simplified interface to reduce the standard work, which has to be done if you work with photometric surveys, to a minimal amount.

## Example

### load data
At the moment a set of coordinates is needed to start. In future releases we want to implement system to perform initial queries without an initial coordinate file.
However, you can load the coordinates with load_coordinates diretly.
```python
from phosphorpy.data.data import DataSet

# load a csv file with two columns (RA and Dec)
ds = DataSet.load_coordinates('coordinates.csv', format='csv')

```
### download data from surveys
To do one of the main task, the mining of photometric surveys, we can use the method 'load_from_vizier'. Currently, we have implement the largest current surveys like 2MASS, UKIDSS, VIKING, KiDS, SDSS, PAN-STARRS, GALEX and WISE. Just put the name of the required survey as a string as parameter into the method and all the data are loaded automatically.
```python
# load the magnitudes from 2MASS and PAN-STARRS of the sources at that coordinates
ds.load_from_vizier('2mass')
ds.load_from_vizier('panstarrs')

```
### flux
From time to time, it is easier to work with the fluxes than with magnitudes. Therefore we have also a simple way to get the fluxes from the magnitudes. The converting uses the bandpass information from SVO.
```python
# get the flux
flux = ds.flux

# plot the SED of the first source with log-log axis
flux.plot.sed(1, x_log=True, y_log=True)
```
