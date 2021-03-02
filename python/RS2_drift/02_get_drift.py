#!/usr/bin/env python
# coding: utf-8

import os
from argparse import ArgumentParser
import numpy as np
import pyproj
import datetime as dt
import pandas as pd

from scipy.ndimage import median_filter
from scipy.ndimage.morphology import distance_transform_edt
from matplotlib.tri import Triangulation

import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from shapely.geometry import Polygon

from osgeo import gdal
from string import Template

from nansat import Nansat, Domain, NSR, Figure
from sea_ice_drift import get_n
from sea_ice_drift.lib import get_spatial_mean
from sea_ice_drift.ftlib import feature_tracking
from sea_ice_drift.pmlib import pattern_matching

# output drift in neXtSIM projection
NS_SRS = f'+proj=stere +lat_0={90} +lat_ts={60} +lon_0={-45} +x_0=0 +y_0=0 +a=6378273 +b=6356889.44891059 +units=m +no_defs'
NS_PROJ = pyproj.Proj(NS_SRS)
NS_GLOBE = ccrs.Globe(semimajor_axis=6378273, semiminor_axis=6356889.44891059)
NS_CRS = ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=60,
        globe=NS_GLOBE)
# hi-res landmask
LAND_50M = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
def parse_args():
    parser = ArgumentParser('Extract sea ice drift from pairs of Radarsat 2 images')
    parser.add_argument('outdir', type=str,
            help='Where to save the figures')
    parser.add_argument('-t', '--test', action='store_true',
            help='Just run for one day for testing')
    parser.add_argument('-urf', '--use-real-features', action='store_true',
            help='Use real features for feature tracking')
    parser.add_argument('-qt', '--quality-threshold', type=float,
            default=2, help='Threshold on quality, require rpm*hpm > this value')
    return parser.parse_args()

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

def fake_feature_tracking(n1, n2):
    print('Creating fake features')
    # create fake feature tracking points
    n1rows, n1cols = n1.shape()
    n2rows, n2cols = n2.shape()

    stp = 500
    c1, r1 = [a.flatten() for a in np.meshgrid(
        np.arange(100, n1cols-100, stp),
        np.arange(100, n1rows-100, stp),
        )]
    c1lon, r1lat = n1.transform_points(c1, r1)
    c2, r2 = n2.transform_points(c1lon, r1lat, DstToSrc=1)
    gpi = (c2 > 0) * (c2 < n2cols) * (r2 > 0) * (r2 < n2rows)
    return c1[gpi], r1[gpi], c2[gpi], r2[gpi]

def fill_nan_gaps(array, distance=15):
    """ Fill gaps in input raster

    Parameters
    ----------
    array : 2D numpy.array
        Input ratser with nan values
    distance : int
        Minimum size of gap to fill

    Returns
    -------
    array : 2D numpy.array
        Ratser with nan gaps field with nearest neigbour values

    """
    array = np.array(array)
    dist, indi = distance_transform_edt(
        np.isnan(array),
        return_distances=True,
        return_indices=True)
    gpi = dist <= distance
    r,c = indi[:,gpi]
    array[gpi] = array[r,c]
    return array

def clean_velo_field(args, a, rpm, hpm, fill_size=1, med_filt_size=3):
    """ Replace gaps with median filtered values """
    a2 = np.array(a)
    a2[(rpm * hpm) < args.quality_threshold] = np.nan
    a3 = fill_nan_gaps(a2, fill_size)
    a4 = median_filter(a3, med_filt_size)
    gpi = np.isnan(a2) * np.isfinite(a4) 
    a2[gpi] = a4[gpi]
    return a2

def get_deformation_elems(x, y, u, v, a):
    """ Compute deformation for given elements.
    Input X, Y, U, V are organized in three columns: for each node of M elements.
    To convert deformation rates from 1/s to %/day outputs should be multiplied by 8640000.
    Parameters
    ----------
    x : 3xM ndarray
        X-coordinates of nodes, m
    y : 3xM ndarray
        Y-coordinates of nodes, m
    u : 3xM ndarray
        U-component of nodes, m/s
    v : 3xM ndarray
        V-component of nodes, m/s
    a : Mx1 ndarray
        area of elements, m2
    Returns
    -------
    e1 : Mx1 array
        Divergence, 1/s
    e2 : Mx1 array
        Shear, 1/s
    e3 : Mx1 array
        Total deformation, 1/s
    """
    # contour integrals of u and v [m/s * m ==> m2/s]
    ux = uy = vx = vy = 0
    for i0, i1 in zip([1, 2, 0], [0, 1, 2]):
        ux += (u[i0] + u[i1]) * (y[i0] - y[i1])
        uy -= (u[i0] + u[i1]) * (x[i0] - x[i1])
        vx += (v[i0] + v[i1]) * (y[i0] - y[i1])
        vy -= (v[i0] + v[i1]) * (x[i0] - x[i1])
    # divide integral by double area [m2/s / m2 ==> 1/day]
    ux, uy, vx, vy =  [i / (2 * a) for i in (ux, uy, vx, vy)]

    # deformation components
    e1 = ux + vy
    e2 = ((ux - vy) ** 2 + (uy + vx) ** 2) ** 0.5
    e3 = np.hypot(e1, e2)
    return e1, e2, e3

def get_deformation_nodes(x, y, u, v):
    """ Compute deformation for given nodes.
    Input X, Y, U, V are given for individual N nodes. Nodes coordinates are triangulated and
    area, perimeter and deformation is computed for M elements.
    Parameters
    ----------
    x : Nx1 ndarray
        X-coordinates of nodes, m
    y : Nx1 ndarray
        Y-coordinates of nodes, m
    u : Nx1 ndarray
        U-component of nodes, m/s
    v : Nx1 ndarray
        V-component of nodes, m/s
    Returns
    -------
    e1 : Mx1 array
        Divergence, 1/s
    e2 : Mx1 array
        Shear, 1/s
    e3 : Mx1 array
        Total deformation, 1/s
    a : Mx1 array
        Area, m2
    p : Mx1 array
        Perimeter, m
    t : 3xM array
        Triangulation (indices of input nodes for each element)
    """
    tri = Triangulation(x, y)

    # coordinates and speeds of corners of each element
    xt, yt, ut, vt = [i[tri.triangles].T for i in (x, y, u, v)]

    # side lengths (X,Y,tot)
    tri_x = np.diff(np.vstack([xt, xt[0]]), axis=0)
    tri_y = np.diff(np.vstack([yt, yt[0]]), axis=0)
    tri_s = np.hypot(tri_x, tri_y)
    # perimeter
    tri_p = np.sum(tri_s, axis=0)
    s = tri_p/2
    # area
    tri_a = np.sqrt(s * (s - tri_s[0]) * (s - tri_s[1]) * (s - tri_s[2]))

    # deformation components
    e1, e2, e3 = get_deformation_elems(xt, yt, ut, vt, tri_a)

    return e1, e2, e3, tri_a, tri_p, tri.triangles

def get_filename(args, t, index, n1, n2):
    fmt = '%Y%m%dT%H%M%SZ'
    return os.path.join(
            args.outdir,
            t.substitute(dict(
                index = index,
                dto1 = n1.time_coverage_start.strftime(fmt),
                dto2 = n2.time_coverage_start.strftime(fmt),
                )),
            )

def save_fig(args, fig, t, index, n1, n2):
    figname = get_filename(args, t, index, n1, n2)
    os.makedirs(os.path.dirname(figname), exist_ok=True)
    print(f'Saving {figname}')
    fig.savefig(figname, bbox_inches='tight')
    plt.close()

def save_npz(args, t, index, n1, n2, **kwargs):
    fname = get_filename(args, t, index, n1, n2)
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    print(f'Saving {fname}')
    np.savez(fname, **kwargs)

def load_npz(args, t, index, n1, n2):
    fname = get_filename(args, t, index, n1, n2)
    if os.path.exists(fname):
        print(f'Loading {fname}')
        return dict(np.load(fname))

# image in stereographic projection
def get_projected_hh(n, d, **kwargs):
    n.reproject(d, **kwargs)
    n1_hh_pro = n[1]
    n.undo()
    return n1_hh_pro.astype(float)

def plot_ft(args, index, n1ft, n2ft, c1, r1, c2, r2):
    lon1b, lat1b = n1ft.get_border()
    lon2b, lat2b = n2ft.get_border()
    x1b, y1b = NS_PROJ(lon1b, lat1b)
    x2b, y2b = NS_PROJ(lon2b, lat2b)
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
    t = Template('ft_keypoints/ft_keypoints_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1ft, n2ft)

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
    plt.quiver(c1, r1, dc, dr, color='r', angles='xy', scale_units='xy')#, scale=0.2)
    plt.plot(n1lon2b, n1lat2b, 'k.-')
    t = Template('ft_quiver/ft_quiver_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1ft, n2ft)

def plot_pm(args, index, n1ft, n2ft, lon1pm, lat1pm, upm, vpm, apm, rpm, hpm):
    # compute ice drift speed [m/s]
    delta_t = (n2ft.time_coverage_start - n1ft.time_coverage_start).total_seconds()
    u = upm / delta_t
    v = vpm / delta_t
    # start points in stereoprojection
    x1pm, y1pm = NS_PROJ(lon1pm, lat1pm)

    # quality (correlation * peakiness)
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(121)
    quality = rpm*hpm
    im = ax.imshow(quality)
    ax.contour(quality, levels=[args.quality_threshold])
    fig.colorbar(im, shrink=.4)
    ax.set_title('rpm*hpm')
    #
    ax = fig.add_subplot(122)
    im = ax.imshow(u)
    ax.contour(quality, levels=[args.quality_threshold])
    ax.set_title('x velocity, m/s')
    fig.colorbar(im, shrink=.4)
    t = Template('pm_quality/pm_quality_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1ft, n2ft)

    # plot final drift on image 1
    res = 2000 #m
    xav = x1pm.mean()
    yav = y1pm.mean()
    xmin, xmax, ymin, ymax = x1pm.min(), x1pm.max(), y1pm.min(), y1pm.max()
    dx = .2*(xmax-xav)
    dy = .2*(ymax-yav)
    xmin -= dx
    xmax += dx
    ymin -= dy
    ymax += dy
    d = Domain(NS_SRS, '-te %f %f %f %f -tr %f %f' % (
        xmin, ymin, xmax, ymax, res, res))
    d_extent = [xmin, xmax, ymin, ymax]
    hh = n1ft[1]
    (vmin,), (vmax,) = Figure(hh).clim_from_histogram(ratio=.9)
    # n1pro = get_projected_hh(n1ft, d)# nearest neighbour
    hh = get_projected_hh(n1ft, d, resample_alg=1)#1:bilinear - better than NN, tps=True
    hh[hh==0] = np.nan

    # plot valid vectors in stereographic projection
    gpi = (quality > args.quality_threshold)
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=NS_CRS)
    im = ax.imshow(hh, vmin=vmin, vmax=vmax, cmap='gray',
            origin='upper', zorder=0, extent=d_extent)
    # plt.colorbar(im, shrink=0.5)
    quiv = ax.quiver(x1pm[gpi], y1pm[gpi], u[gpi], v[gpi], rpm[gpi])#, scale=2)
    fig.colorbar(quiv, shrink=0.5)
    #plt.quiverkey(quiv, 110000, -770000, 0.05, '0.05 m/s', coordinates='data')
    ax.set_title('Ice drift speed [m/s]')
    ax.add_feature(LAND_50M, edgecolor='black')

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    ax.gridlines(linewidth=2, color='m', 
        draw_labels=True, alpha=0.5, linestyle=':')
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LONGITUDE_FORMATTER)
    t = Template('pm_quiver/pm_quiver_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1ft, n2ft)

    # deformation
    e1, e2, e3, tri_a, tri_p, triangles = get_deformation_nodes(x1pm[gpi], y1pm[gpi], u[gpi], v[gpi])
    # convert deformations to %/day
    e1 *= 8640000
    e2 *= 8640000
    e3 *= 8640000
    for e,t,ttl in zip(
            [e1,e2,e3],
            [Template('shear/shear_${dto1}-${dto2}_${index}.png'),
                Template('divergence/divergence_${dto1}-${dto2}_${index}.png'),
                Template('total-defor/total-defor_${dto1}-${dto2}_${index}.png'),
                ],
            ['Shear [%/day]', 'Divergence [%/day]', 'Total deformation [%/day]']):
        vmin, vmax = np.percentile(e, [20,80])
        fig = plt.figure(figsize=(20,20))
        ax = plt.axes(projection=NS_CRS)
        im = ax.tripcolor(x1pm[gpi], y1pm[gpi],
                          triangles=triangles, facecolors=e,
                          vmin=vmin, vmax=vmax,
                          edgecolors='k', linewidth=1)
        fig.colorbar(im, shrink=0.5)
        ax.set_title(ttl)
        ax.add_feature(LAND_50M, edgecolor='black')

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

        ax.gridlines(linewidth=2, color='m',
            draw_labels=True, alpha=0.5, linestyle=':')
        ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LONGITUDE_FORMATTER)
        save_fig(args, fig, t, index, n1ft, n2ft)

def plot_pm_clean(args, index, n1ft, n2ft, lon1pm, lat1pm, upm, vpm, apm, rpm, hpm):
    # compute ice drift speed [m/s]
    delta_t = (n2ft.time_coverage_start - n1ft.time_coverage_start).total_seconds()
    u = upm / delta_t
    v = vpm / delta_t
    u2 = clean_velo_field(args, u, rpm, hpm)
    v2 = clean_velo_field(args, v, rpm, hpm)

    # start points in stereoprojection
    x1pm, y1pm = NS_PROJ(lon1pm, lat1pm)

    fig = plt.figure()
    ax = fig.add_subplot(121)
    im = ax.imshow(u)
    ax.set_title('x velocity, m/s')
    fig.colorbar(im, shrink=.4)

    ax = fig.add_subplot(122)
    im = ax.imshow(u2)
    ax.set_title('Cleaned x velocity, m/s')
    fig.colorbar(im, shrink=.4)

    t = Template('clean-u/u_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1ft, n2ft)

    # plot final drift on image 1
    res = 2000 #m
    xav = x1pm.mean()
    yav = y1pm.mean()
    xmin, xmax, ymin, ymax = x1pm.min(), x1pm.max(), y1pm.min(), y1pm.max()
    dx = .2*(xmax-xav)
    dy = .2*(ymax-yav)
    xmin -= dx
    xmax += dx
    ymin -= dy
    ymax += dy
    d = Domain(NS_SRS, '-te %f %f %f %f -tr %f %f' % (
        xmin, ymin, xmax, ymax, res, res))
    d_extent = [xmin, xmax, ymin, ymax]
    hh = n1ft[1]
    (vmin,), (vmax,) = Figure(hh).clim_from_histogram(ratio=.9)
    # n1pro = get_projected_hh(n1ft, d)# nearest neighbour
    hh = get_projected_hh(n1ft, d, resample_alg=1)#1:bilinear - better than NN, tps=True
    hh[hh==0] = np.nan

    # plot valid vectors in stereographic projection
    gpi = np.isfinite(u2*v2)
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=NS_CRS)
    im = ax.imshow(hh, vmin=vmin, vmax=vmax, cmap='gray',
            origin='upper', zorder=0, extent=d_extent)
    # plt.colorbar(im, shrink=0.5)
    quiv = ax.quiver(x1pm[gpi], y1pm[gpi], u2[gpi], v2[gpi], rpm[gpi])#, scale=2)
    fig.colorbar(quiv, shrink=0.5)
    #plt.quiverkey(quiv, 110000, -770000, 0.05, '0.05 m/s', coordinates='data')
    ax.set_title('Ice drift speed [m/s]')
    ax.add_feature(LAND_50M, edgecolor='black')

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    ax.gridlines(linewidth=2, color='m', 
        draw_labels=True, alpha=0.5, linestyle=':')
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LONGITUDE_FORMATTER)
    t = Template('clean-quiver/quiver_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1ft, n2ft)

    # deformation
    e1, e2, e3, tri_a, tri_p, triangles = get_deformation_nodes(x1pm[gpi], y1pm[gpi], u2[gpi], v2[gpi])
    # convert deformations to %/day
    e1 *= 8640000
    e2 *= 8640000
    e3 *= 8640000
    for e,t,ttl in zip(
            [e1,e2,e3],
            [
                Template('clean-shear/shear_${dto1}-${dto2}_${index}.png'),
                Template('clean-divergence/divergence_${dto1}-${dto2}_${index}.png'),
                Template('clean-total-defor/total-defor_${dto1}-${dto2}_${index}.png'),
                ],
            ['Shear [%/day]', 'Divergence [%/day]', 'Total deformation [%/day]']):
        vmin, vmax = np.percentile(e, [20,80])
        fig = plt.figure(figsize=(20,20))
        ax = plt.axes(projection=NS_CRS)
        im = ax.tripcolor(x1pm[gpi], y1pm[gpi],
                          triangles=triangles, facecolors=e,
                          vmin=vmin, vmax=vmax,
                          edgecolors='k', linewidth=1)
        fig.colorbar(im, shrink=0.5)
        ax.set_title(ttl)
        ax.add_feature(LAND_50M, edgecolor='black')

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])

        ax.gridlines(linewidth=2, color='m',
            draw_labels=True, alpha=0.5, linestyle=':')
        ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
        ax.yaxis.set_major_formatter(LONGITUDE_FORMATTER)
        save_fig(args, fig, t, index, n1ft, n2ft)

def process_1pair(args, f1, f2, index):
    # create Nansat objects with one band only. 
    n1ft, n1pm = get_nansat(f1)
    n2ft, n2pm = get_nansat(f2)
    print(f'1st image start: {n1ft.time_coverage_start}')
    print(f'2nd image start: {n2ft.time_coverage_start}')

    if args.use_real_features:
        # Run Feature Tracking
        # - get start/end coordinates in the image coordinate system (colums/rows)
        # - works best with spatial mean removed (use `n1ft`, `n2ft`)
        if 1:
            #v1: feature tracking takes HH with spatial mean removed
            c1, r1, c2, r2 = feature_tracking(n1ft, n2ft, nFeatures=100000,
                    ratio_test=0.6, domainMargin=0)
        else:
            #v1: feature tracking takes HH without spatial mean removed
            c1, r1, c2, r2 = feature_tracking(n1pm, n2pm, nFeatures=100000,
                    ratio_test=0.6, domainMargin=0)

        t = Template('npz_files/ft_${dto1}-${dto2}_${index}.npz')
        save_npz(args, t, index, n1ft, n2ft, c1=c1, r1=r1, c2=c2, r2=r2)
        plot_ft(args, index, n1ft, n2ft, c1, r1, c2, r2)
    else:
        # Fake feature tracking to get more starting points for pattern matching
        c1, r1, c2, r2 = fake_feature_tracking(n1pm, n2pm)
    lon1pm, lat1pm = n1ft.get_geolocation_grids(200)

    # Pattern matching
    # lon/lat grids for image_1
    t = Template('npz_files/pm_${dto1}-${dto2}_${index}.npz')
    pm_results = load_npz(args, t, index, n1ft, n2ft)
    if pm_results is None:
        print('Running pattern matching...')
        pm_results = dict()
        # Run Pattern Matching for each element in lon1pm/lat1pm matrix
        # ice displacement upm and vpm are returned in meters in Stereographic projection
        (
                pm_results['upm'], pm_results['vpm'],
                pm_results['apm'], pm_results['rpm'], pm_results['hpm'],
                lon2pm, lat2pm,
                ) = pattern_matching(
                        lon1pm, lat1pm, n1pm, c1, r1, n2pm, c2, r2,
                        min_border=300, max_border=300,
                        img_size=50, srs=NS_SRS)
        pm_results['upm_clean'] = clean_velo_field(args,
                pm_results['upm'], pm_results['rpm'], pm_results['hpm'])
        pm_results['vpm_clean'] = clean_velo_field(args,
                pm_results['vpm'], pm_results['rpm'], pm_results['hpm'])
        save_npz(args, t, index, n1ft, n2ft, **pm_results)

    plot_pm(args, index, n1ft, n2ft,
                lon1pm, lat1pm,
                pm_results['upm'], pm_results['vpm'],
                pm_results['apm'], pm_results['rpm'], pm_results['hpm'],
                )
    plot_pm_clean(args, index, n1ft, n2ft,
                lon1pm, lat1pm,
                pm_results['upm'], pm_results['vpm'],
                pm_results['apm'], pm_results['rpm'], pm_results['hpm'],
                )

def run():
    args = parse_args()
    rs2dir = os.getenv('RS2_dir')

    if args.test:
        # if testing just run 1 example
        f1 = 'RS2_OK37499_PK364900_DK322205_SCWA_20130224_023727_HH_HV_SGF'
        f2 = 'RS2_OK37500_PK364951_DK322250_SCWA_20130225_020811_HH_HV_SGF'
        process_1pair(args,
                os.path.join(rs2dir, f1),
                os.path.join(rs2dir, f2),
                0)
        return

    df = pd.read_csv('out/RS2_pairs.csv',
            names=["File1","File2","Interval","Overlap"], skiprows=1,
            dtype={"File1":str, "File2":str, "Interval":float, "Overlap":float})
    for index, (f1,  f2, interval, overlap) in df.iterrows():
        print(f'Processing pair {index} out of {len(df)}')
        print(f1, f2, sep='\n')
        print(f'Overlap = {overlap}')
        print(f'Time interval = {interval}')
        try:
            process_1pair(args,
                    os.path.join(rs2dir, f1),
                    os.path.join(rs2dir, f2),
                    index)
        except:
            continue

if __name__ == "__main__":
    run()
