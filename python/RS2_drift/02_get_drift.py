#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.tri import Triangulation

from shapely.geometry import Polygon
import pyproj
import datetime as dt
import pandas as pd
from osgeo import gdal
from string import Template

from nansat import Nansat, Domain, NSR, Figure

from sea_ice_drift import get_n
from sea_ice_drift.lib import get_spatial_mean
from sea_ice_drift.ftlib import feature_tracking
from sea_ice_drift.pmlib import pattern_matching

# output drift in neXtSIM projection
NS_SRS = f'+proj=stere +lat_0={90} +lat_ts={60} +lon_0={-45} +x_0=0 +y_0=0 +a=6378273 +b=6356889.44891059 +units=m +no_defs'
NS_PROJ = pyproj.Proj(srs)
NS_GLOBE = ccrs.Globe(semimajor_axis=6378273, semiminor_axis=6356889.44891059)
NS_CRS = ccrs.PolarStereographic(central_latitude=90, central_longitude=-45,
                         true_scale_latitude=60, globe=globe)
# hi-res landmask
LAND_50M = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])


def get_nansat(f):
    """
    Returns:
    --------
    nft : nansat.Nansat
        sigma0_HH with spatial mean removed (helps feature tracking and plotting)
    npm: nansat.Nansat
        sigma0_HH without spatial mean removed (so pattern matching doesn't crash)
    """
    if 0:
        n = get_n(f, bandName='sigma0_HH', remove_spatial_mean=True)
        return n, n
    n = Nansat(f)
    start = n.time_coverage_start
    ds = gdal.Open(n.filename+'/imagery_HH.tif')
    hh = ds.ReadAsArray()
    npm = Nansat.from_domain(n, array=hh)
    hh -= get_spatial_mean(hh).astype('uint8') #remove mean
    nft = Nansat.from_domain(n, array=hh)
    for n in [nft, npm]:
        n.set_metadata(dict(time_coverage_start=start))
    return nft, npm

def get_bbox_approx():
    lon, lat = np.array([
        (-109.4, 74.7), # TR
        (-118.8, 67.5), # BR
        (-161.3, 70.5), # BL: Barrow
        (-153.9, 75.8), # TL
        ]).T
    return lon, lat


def get_bbox():
    lon, lat = get_bbox_approx()
    x, y = NS_PROJ(lon, lat)
    xav, yav = np.mean(x), np.mean(y)
    dx = x.max() - x.min()
    dy = y.max() - y.min()
    factor = 1.4
    return xav - factor*dx/2, xav + factor*dx/2, yav - factor*dy/2, yav + factor*dy/2

def get_filename(t, n1, n2):
    fmt = '%Y%m%dT%H%M%SZ'
    return t.substitute(dict(
        dto1 = n1.time_coverage_start.strftime(fmt),
        dto2 = n2.time_coverage_start.strftime(fmt),
        ))

def save_fig(fig, t, n1, n2):
    figname = get_filename(t, n1, n2)
    os.makedirs(os.path.basename(figname), exist_ok=True)
    fig.savefig(figname, bbox_inches='tight')
    plt.close()

def save_npz(t, n1, n2, **kwargs):
    fname = get_filename(t, n1, n2)
    os.makedirs(os.path.basename(fname), exist_ok=True)
    np.savez(fname, **kwargs)

# image in stereographic projection
def get_projected_hh(n, d, **kwargs):
    n.reproject(d, **kwargs)
    n1_hh_pro = n[1]
    n.undo()
    return n1_hh_pro.astype(float)

def plot_ft(n1ft, n2ft, c1, r1, c2, r2):
    lon1b, lat1b = n1ft.get_border()
    lon2b, lat2b = n2ft.get_border()
    x1b, y1b = NS_PROJ(lon1b, lat1b)
    x2b, y2b = PROJ_NS(lon2b, lat2b)
    # convert row/column coordinates of matched features to lon/lat
    lon1ft, lat1ft = n1ft.transform_points(c1, r1)
    lon2ft, lat2ft = n2ft.transform_points(c2, r2)
    # plot keypoints
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70))
    ax.add_feature(LAND_50M, zorder=0, edgecolor='black')
    ax.plot(lon1ft, lat1ft, '.', label='keypoints_1', transform=ccrs.PlateCarree())
    ax.plot(lon2ft, lat2ft, '.', label='keypoints_2', transform=ccrs.PlateCarree())
    ax.plot(lon1b, lat1b, '.-', label='border_1', transform=ccrs.PlateCarree())
    ax.plot(lon2b, lat2b, '.-', label='border_2', transform=ccrs.PlateCarree())
    ax.legend()
    t = Template('out/ft_keypoints/ft_keypoints_${dto1}-${dto2}.png')
    save_fig(fig, t, n1ft, n2ft)

    # Plot ice drift on top of image_1
    # - end points in image_1 coordinate system
    n1c2, n1r2 = n1ft.transform_points(lon2ft, lat2ft, DstToSrc=1)
    # - ice drift components in image_1 coordinate system
    dr = n1r2 - r1
    dc = n1c2 - c1
    # - border of image_2 in image_1 coordinate system
    n1lon2b, n1lat2b = n1ft.transform_points(lon2b, lat2b, DstToSrc=1)
    # plot of ice drift. Arrows are 5 times longer than actual drift
    fig = plt.figure(figsize=(10,10))
    hh = n1ft[1]
    (vmin,), (vmax,) = Figure(hh).clim_from_histogram(ratio=.9)
    plt.imshow(hh, cmap='gray', vmin=vmin, vmax=vmax)
    plt.quiver(c1, r1, dc, dr, color='r', angles='xy', scale_units='xy', scale=0.2)
    plt.plot(n1lon2b, n1lat2b, 'k.-')
    t = Template('out/ft_quiver/ft_quiver_${dto1}-${dto2}.png')
    save_fig(fig, t, n1ft, n2ft)

def plot_pm(n1ft, n2ft, lon1pm, latpm, upm, vpm, apm, rpm, hpm):
    # compute ice drift speed [m/s]
    delta_t = (n2ft.time_coverage_start - n1ft.time_coverage_start).total_seconds()
    u = upm / delta_t
    v = vpm / delta_t
    # start points in stereoprojection
    x1pm, y1pm = proj(lon1pm, lat1pm)

    # quality (correlation * peakiness)
    fig = plt.figure((10,10))
    ax = fig.add_subplot(121)
    q = rpm*hpm
    im = ax.imshow(q)
    ax.contour(q, levels=[2.5, 3, 3.5, 4], inline=True)
    fig.colorbar(im)
    ax.set_title('rpm*hpm')
    #
    ax = fig.add_subplot(122)
    im = ax.imshow(u)
    ax.contour(q, levels=[2.5, 3, 3.5, 4], inline=True)
    ax.title('x velocity, m/s')
    fig.colorbar(im)
    t = Template('out/pm_quality/pm_quality_${dto1}-${dto2}.png')
    save_fig(fig, t, n1ft, n2ft)

    # plot final drift on image 1
    res = 2000 #m
    xmin, xmax, ymin, ymax = x1pm.min(), x1pm.max(), y1pm.min(), y1pm.max()
    d = Domain(NS_SRS, '-te %f %f %f %f -tr %f %f' % (
        xmin, ymin, xmax, ymax, res, res))
    d_extent = [xmin, xmax, ymin, ymax]
    hh = n1ft[1]
    (vmin,), (vmax,) = Figure(hh).clim_from_histogram(ratio=.9)
    # n1pro = get_projected_hh(n1ft, d)# nearest neighbour
    hh = get_projected_hh(n1ft, d, resample_alg=1)#1:bilinear - better than NN, tps=True

    # plot valid vectors in stereographic projection
    gpi = (q > 3.5) 
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=NS_CRS)
    im = ax.imshow(hh, vmin=vmin, vmax=vmax, cmap='gray',
            origin='upper', zorder=0, extent=d_extent)
    # plt.colorbar(im, shrink=0.5)
    quiv = ax.quiver(x1pm[gpi], y1pm[gpi], u[gpi], v[gpi], rpm[gpi], scale=2)
    fig.colorbar(quiv, shrink=0.5)
    #plt.quiverkey(quiv, 110000, -770000, 0.05, '0.05 m/s', coordinates='data')
    ax.set_title('Ice drift speed [m/s]')
    ax.add_feature(LAND_50M, edgecolor='black')

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    ax.gridlines(linewidth=2, color='m', 
        draw_labels=True, alpha=0.5, linestyle=':')
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    t = Template('out/pm_quiver/pm_quiver_${dto1}-${dto2}.png')
    save_fig(fig, t, n1ft, n2ft)

def process_1pair(f1, f2):
    print(f'1st image start: {n1ft.time_coverage_start}')
    print(f'2nd image start: {n2ft.time_coverage_start}')

    # create Nansat objects with one band only. 
    n1ft, n1pm = get_nansat(f1)
    n2ft, n2pm = get_nansat(f2)

    # Run Feature Tracking
    # - get start/end coordinates in the image coordinate system (colums/rows)
    # - works best with spatial mean removed (use `n1ft`, `n2ft`)
    c1, r1, c2, r2 = feature_tracking(n1ft, n2ft, nFeatures=100000,
            ratio_test=0.6, domainMargin=0)
    t = Template('out/npz_files/ft_${dto1}-${dto2}.npz')
    save_npz(t, n1ft, n2ft, c1=c1, r1=r1, c2=c2, r2=r2)
    plot_ft(n1ft, n2ft, c1, r1, c2, r2)

    # Pattern matching
    # lon/lat grids for image_1
    lon1pm, lat1pm = n1ft.get_geolocation_grids(200)
    # Run Pattern Matching for each element in lon1pm/lat1pm matrix
    # ice displacement upm and vpm are returned in meters in Stereographic projection
    upm, vpm, apm, rpm, hpm, lon2pm, lat2pm = pattern_matching(
        lon1pm, lat1pm, n1pm, c1, r1, n2pm, c2, r2,
        min_border=300, max_border=300,
        img_size=50, srs=NS_SRS)
    t = Template('out/npz_files/pm_${dto1}-${dto2}.npz')
    save_npz(t, n1ft, n2ft, upm=upm, vpm=vpm, apm=apm, rpm=rpm, hpm=hpm)
    plot_pm(n1ft, n2ft, upm, vpm, apm, rpm, hpm)

def run():
    df = pd.read_csv('out/RS2_pairs.csv',
            names=["File1","File2","Interval","Overlap"], skiprows=1,
            dtype={"File1":str, "File2":str, "Interval":float, "Overlap":float})
    for index, (f1,  f2, interval, overlap) in df.iterrows():
        print(f'Processing pair {index} out of {len(df)}')
        print(f1, f2, sep='\n')
        print(f'Overlap = {overlap}')
        print(f'Time interval = {interval}')
        process_1pair(f1, f2)

if __name__ == "__main__":
    run()
