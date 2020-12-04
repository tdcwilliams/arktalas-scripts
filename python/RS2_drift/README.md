# Install sea ice drift with conda
```
conda create -q --yes -n py3drift -c conda-forge python=3.6 numpy scipy matplotlib cartopy \
      netcdf4 cartopy gdal opencv nose ipython jupyter \
&& conda activate py3drift
```

```
export PYTHONPATH=$HOME/Github-Repos/nansat:$HOME/Github-Repos/sea_ice_drift:$PYTHONPATH
cd $HOME/Github-Repos/nansat
python setup.py install
```
