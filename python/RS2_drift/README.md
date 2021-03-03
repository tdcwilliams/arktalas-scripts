# Download RS2 data
They are on johansen at `/Data/sat/downloads/Radarsat2/beaufort_sea/2013/`.
Put the files in `$RS2_dir`.

# Install sea ice drift with conda
```
conda create -q --yes -n py3drift -c conda-forge python=3.6 numpy scipy matplotlib \
      netcdf4 cartopy gdal opencv \
      pandas pyyaml requests xdg pyproj \
      nose ipython jupyter mock \
&& conda activate py3drift
```

## Install nansat
```
git clone https://github.com/nansencenter/nansat.git
export PYTHONPATH=$HOME/Github-Repos/nansat:$PYTHONPATH
cd $HOME/Github-Repos/nansat
python setup.py build_ext --inplace
```
Download MODIS landmask from ftp://ftp.nersc.no/nansat/MOD44W.tgz - untar to directory `$MOD44WPATH`.

Until gdal fix is merged checkout this branch:
```
git checkout hotfix493-import-gdal
```
Also need `pythesint`: if `pip install https://github.com/nansencenter/py-thesaurus-interface/archive/master.tar.gz
` doesn't work, install manually with:
```
git clone https://github.com/nansencenter/py-thesaurus-interface.git
export PYTHONPATH=$HOME/Github-Repos/py-thesaurus-interface:$PYTHONPATH
```
Check installation with
```
nosetests nansat
```

## Install sea_ice_drift
```
git clone https://github.com/nansencenter/sea_ice_drift.git
export PYTHONPATH=$HOME/Github-Repos/sea_ice_drift:$PYTHONPATH
```
Until gdal fix is merged checkout this branch:
```
git checkout issue24-gdal
```
Download the following example files from colhub.met.no and put them in `S1B_dir`:
```
S1B_EW_GRDM_1SDH_20200123T120618_20200123T120718_019944_025BA1_D4A2
S1B_EW_GRDM_1SDH_20200125T114955_20200125T115055_019973_025C81_EC1A
```

# Match scenes that are close enough in time
`./01_match_scenes.py`
Saves results to `out/RS2_pairs.csv`.

# Calculate drift and make plots
`./02_get_drift.py out`
Saves results to `out/npz_files` and saves some other plots to `out` directory.
