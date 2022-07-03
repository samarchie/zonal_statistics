<br />
<p align="center">
  <h1 align="center">Zonal Statistics</h3>
</p>
<br />

## About The Project

Collect statistics (minimum, maximum, mean, median, standard deviation, sum and count) on pixel values (raster) which lie within polygons (vector). If allow_touching_pixels is set to True, then any pixel that touches the exterior of the vector geometry is included in the statistics calculation. The returned dictionary has the vector FID as the key, where the values is a dictionary mapping the statistic's name to value.

This script is better than previously used packages [pygeoprocessing](https://github.com/natcap/pygeoprocessing) and [rasterstats](https://github.com/perrygeo/python-rasterstats) because:
1. It is faster than both packages using the same test case
2. The option to allow touching pixels exists, where these two packages did not include this

## Getting Started

This module requires a version of Python 3.8 or above (tested on 3.9.1) and an environment with the packages listed in ```requirements.txt``` (which can be installed by running ```pip install -r requirements.txt```). 

## Acknowledgements

Inspiration is taken from https://opensourceoptions.com/blog/zonal-statistics-algorithm-with-python-in-4-steps/, https://github.com/natcap/pygeoprocessing/blob/c364ba8e2429c841625c95a34cc871883d3a9cfd/src/pygeoprocessing/geoprocessing.py#L1144 and https://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides4.pdf
