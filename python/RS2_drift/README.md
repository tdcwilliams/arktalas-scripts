# Download RS2 data
They are on johansen at `/Data/sat/downloads/Radarsat2/beaufort_sea/2013/`.

# Install sea ice drift with conda
```
conda create -q --yes -n py3drift -c conda-forge python=3.6 numpy scipy matplotlib cartopy \
      netcdf4 cartopy gdal opencv nose ipython jupyter mock \
&& conda activate py3drift
```

## Install nansat
```
export PYTHONPATH=$HOME/Github-Repos/nansat:$PYTHONPATH
cd $HOME/Github-Repos/nansat
python setup.py build_ext --inplace
```
Download MODIS landmask from ftp://ftp.nersc.no/nansat/MOD44W.tgz - untar to directory `$MOD44WPATH`.

Until gdal fix is merged checkout this branch:
```
git checkout hotfix493-import-gdal
```

## Install sea_ice_drift
```
export PYTHONPATH=$HOME/Github-Repos/sea_ice_drift:$PYTHONPATH
```
Until gdal fix is merged checkout this branch:
```
git checkout issue24-gdal
```

# Match scenes that are close enough in time
`./01_match_scenes.py`
Saves results to `out/RS2_pairs.csv`.
