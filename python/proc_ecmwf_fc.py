#!/usr/bin/env python
'''
INPUT FILES:

analysis
hours: 0,6,12,18
od.ans.20130201-20130531.sfc.ci.nc      conc
od.ans.20130201-20130531.sfc.d2m.nc     2D
od.ans.20130201-20130531.sfc.lsm.nc     land-sea mask
od.ans.20130201-20130531.sfc.sp.nc      MSL = pressure
od.ans.20130201-20130531.sfc.t2m.nc     2T
od.ans.20130201-20130531.sfc.u10.nc     10U
od.ans.20130201-20130531.sfc.v10.nc     10V

00:00 forecast
hours: 3,9,15
od.for00.20130201-20130531.sfc.sf.nc    snowfall - DEACCUM!!
od.for00.20130201-20130531.sfc.ssrd.nc  SSRD = short wave down - DEACCUM!!
od.for00.20130201-20130531.sfc.strd.nc  STRD = long wave down - DEACCUM!!
od.for00.20130201-20130531.sfc.tp.nc    TP = total precip - DEACCUM!!

12:00 forecast
hours: 15,21,3
od.for12.20130201-20130531.sfc.sf.nc    snowfall - DEACCUM!!
od.for12.20130201-20130531.sfc.ssrd.nc  SSRD = short wave down - DEACCUM!!
od.for12.20130201-20130531.sfc.strd.nc  STRD = long wave down - DEACCUM!!
od.for12.20130201-20130531.sfc.tp.nc    TP = total precip - DEACCUM!!

1. instantaneous variables are fine
2. accumulated variables:
    0=2*for00:3, 6=for00:9-for00:3, 12=2*for12:15, 18=for12:21-for12:15
'''

import os
import argparse
import numpy as np
import datetime as dt
from string import Template
from netCDF4 import Dataset

# Celsius to Kelvin conversion
_KELVIN = 273.15 # [C]

# filenames
# - output
DST_FILEMASK = os.path.join(
    '/cluster/projects/nn2993k/sim/data/WRF/ECMWF_forecast_arctic',
    'ec2_start%Y%m%d.nc')
DST_REFDATE = dt.datetime(1950, 1, 1) #ref date hard-coded into neXtSIM
# - input
SRC_FILE_TEMPLATE = os.path.join(
    '/cluster/projects/nn2993k/sim/data/WRF/ECMWF_forecast_arctic_clemens',
    'od.${TYPE}${CYCLE}.20130201-20130531.sfc.${VARNAME}.nc')
SRC_REFDATE = dt.datetime(1900, 1, 1) #ref date in downloaded file

# Destination variables
DST_DIMS = {
        'time' : 'time',
        'lat': 'latitude',
        'lon': 'longitude',
        }
DST_VARS = {
    '10U'   : dict(VARNAME='u10',  TYPE='ans', CYCLE=''),
    '10V'   : dict(VARNAME='v10',  TYPE='ans', CYCLE=''),
    '2T'    : dict(VARNAME='t2m',  TYPE='ans', CYCLE=''),
    '2D'    : dict(VARNAME='d2m',  TYPE='ans', CYCLE=''),
    'MSL'   : dict(VARNAME='sp',   TYPE='ans', CYCLE=''),
    'SSRD'  : dict(VARNAME='ssrd', TYPE='for'),
    'STRD'  : dict(VARNAME='strd', TYPE='for'),
    'TP'    : dict(VARNAME='tp',   TYPE='for'),
    'SF'    : dict(VARNAME='sf',   TYPE='for'),
    }
GRIDFILE = Template(SRC_FILE_TEMPLATE).safe_substitute(DST_VARS['10U'])
if 0:
    # test on smaller subset of variables
    #dst_var = 'TP'
    dst_var = '2T'
    DST_VARS = {dst_var: DST_VARS[dst_var]}


VALID_DATE = lambda x : dt.datetime.strptime(x, '%Y%m%d')
KW_COMPRESSION = dict(zlib=True)

def parse_args():
    ''' parse command line arguments '''
    parser = argparse.ArgumentParser(
            description="""
            Split downloaded ECMWF file into daily file
            and deaccumulate the required variables""")
    parser.add_argument('date', type=VALID_DATE,
            help='input date (YYYYMMDD)')
    return parser.parse_args()

def get_time_slice(src_ds, date, ftype):
    '''
    Parameters:
    -----------
    src_ds : netCDF4.Dataset
    date : datetime.datetime
    ftype : str
        ans, for00 or for12

    Returns:
    --------
    t_slice : slice
    '''
    src_time_raw = src_ds.variables['time'][:].astype(float)
    src_time = [SRC_REFDATE + dt.timedelta(hours=h) for h in src_time_raw]
    shifts = dict(ans=0, for00=3, for12=15)
    recs = dict(ans=4, for00=2, for12=2)
    i0 = src_time.index(date + dt.timedelta(hours=shifts[ftype]))
    return slice(i0, i0+recs[ftype])

def get_instantaneous_var(var_info, date):
    '''
    Parameters:
    -----------
    var_info : dict
        info from 1 member of DST_VARS
    date : dt.datetime

    Returns:
    --------
    var : np.ndarray(float)
    atts: dict
        variable attributes
    '''
    fname = Template(SRC_FILE_TEMPLATE).safe_substitute(var_info)
    with Dataset(fname, 'r') as src_ds:
        t_slice = get_time_slice(src_ds, date, 'ans')
        src_var = src_ds.variables[var_info['VARNAME']]
        return src_var[t_slice], vars(src_var)

def get_accumulated_var(var_info, date):
    '''
    Parameters:
    -----------
    var_info : dict
        info from 1 member of DST_VARS
    date : dt.datetime

    Returns:
    --------
    var : np.ndarray(float)
    atts: dict
        variable attributes
    '''
    v = []
    atts = None
    for cycle in ['00', '12']:
        fname = Template(SRC_FILE_TEMPLATE).safe_substitute(dict(**var_info, CYCLE=cycle))
        with Dataset(fname, 'r') as src_ds:
            t_slice = get_time_slice(src_ds, date, f'for{cycle}')
            src_var = src_ds.variables[var_info['VARNAME']]
            v_ = src_var[t_slice]
            v += [2*v_[0], v_[1] - v_[0]] # convert to rate*6h for model
            if atts is None:
                atts = vars(src_var)
    return np.array(v), atts

def get_var(var_name, date):
    '''
    Parameters:
    -----------
    var_name : str
        variable name in output file
    date : dt.datetime

    Returns:
    --------
    var : np.ndarray(float)
    atts: dict
        variable attributes
    '''
    var_info = DST_VARS[var_name]
    if var_info['TYPE'] == 'ans':
        return get_instantaneous_var(var_info, date)
    return get_accumulated_var(var_info, date)

def get_destination_coordinates(date):
    """ Load dimensions

    Parameters
    ----------
    date: dt.datetime

    Returns
    -------
    dst_vec : dict
        three vectors with destination coordinates: time, lat, lon
    """
    # coordinates on destination grid
    time = [date + dt.timedelta(hours=i*6) for i in range(4)]
    time = [(t - DST_REFDATE).days*24 for t in time]
    with Dataset(GRIDFILE, 'r') as ds:
        return {
            'time': np.array(time),
            'lat': ds.variables['latitude'][:],
            'lon': ds.variables['longitude'][:],
        }

def export(outfile, dst_dims, dst_data):
    """
    Export output netcdf file

    Parameters
    ----------
    outfile : str
        netcdf output filename
    dst_dims : dict
        three vectors with destination coordinates (time, lat, lon)
    dst_data : dict
    """
    # Create dataset for output
    skip_var_attr = ['_FillValue', 'grid_mapping']
    # create dataset
    print(f'Exporting {outfile}')
    with Dataset(GRIDFILE, 'r') as src_ds, Dataset(outfile, 'w') as dst_ds:
        # add dimensions
        for dim_name, dim_vec in dst_dims.items():
            dlen = {'time': None}.get(dim_name, len(dim_vec)) #time should be unlimited
            dtype = {'time': 'f8'}.get(dim_name, 'f4') #time should be double
            dst_dim = dst_ds.createDimension(dim_name, dlen)
            dst_var = dst_ds.createVariable(
                    dim_name, dtype, (dim_name,), **KW_COMPRESSION)
            src_var = src_ds.variables[DST_DIMS[dim_name]]
            for ncattr in src_var.ncattrs():
                if [dim_name, ncattr] == ['time', 'units']:
                    units = DST_REFDATE.strftime(
                            'hours since %Y-%m-%d 00:00:00.0') #need to change ref time
                    dst_var.setncattr(ncattr, units)
                else:
                    dst_var.setncattr(ncattr, src_var.getncattr(ncattr))
            dst_var[:] = dim_vec

        # add processed variables
        for dst_var_name, (data, atts) in dst_data.items():
            dst_var = dst_ds.createVariable(dst_var_name, 'f4',
                    ('time', 'lat', 'lon'),
                    **KW_COMPRESSION)
            for att, val in atts.items():
                if att in skip_var_attr:
                    continue
                dst_var.setncattr(att, val)
            dst_var[:] = data

def run(args):
    '''
    make the file

    Parameters:
    -----------
    args : argparse.Namespace
    '''
    # load data
    dst_dims = get_destination_coordinates(args.date)
    dst_data = {v: get_var(v, args.date) for v in DST_VARS}
    # export
    outfile = args.date.strftime(DST_FILEMASK)
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    export(outfile, dst_dims, dst_data)

if __name__ == '__main__':
    run(parse_args())
