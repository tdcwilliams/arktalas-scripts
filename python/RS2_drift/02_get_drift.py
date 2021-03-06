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
import json

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
    parser.add_argument('input_pairs_file', type=str,
            help='csv file with list of image pairs')
    parser.add_argument('outdir', type=str,
            help='Where to save the figures')
    parser.add_argument('-t', '--test', action='store_true',
            help='Just run for one day for testing')
    parser.add_argument('-qt', '--quality-threshold', type=float,
            default=2, help='Threshold on quality, require rpm*hpm > this value')
    parser.add_argument('-f', '--force', action='store_true',
            help='Redo pattern matching even if results are already saved')
    return parser.parse_args()

def get_fixed_n(n, img, time_coverage_start):
    # create new Nansat with one band only
    n2 = Nansat.from_domain(n, img)
    n2.set_metadata('time_coverage_start', time_coverage_start)

    # improve geolocation
    n2.reproject_gcps()
    n2.vrt.tps = True
    return n2

def get_nansat(filename, show=False):
    """
    Returns:
    --------
    nft : nansat.Nansat
        sigma0_HH with spatial mean removed (helps feature tracking and plotting)
    npm: nansat.Nansat
        sigma0_HH without spatial mean removed (so pattern matching doesn't crash)
    """
    n01 = Nansat(filename)
    time_coverage_start = n01.time_coverage_start

    # open with Nansat to find HH bandname
    filename_gdal = f'RADARSAT_2_CALIB:UNCALIB:{filename}/product.xml'
    n02 = Nansat(filename_gdal)

    try:
        band_id = n02.get_band_number({'POLARIMETRIC_INTERP': 'HH'})
    except:
        band_id = 1

    # replace zeros in low values (not near border)
    img = n02[band_id]
    img[100:-100, 100:-100] = np.clip(img[100:-100, 100:-100], 1, 255)
    if show:
        (vmin,), (vmax,) = Figure(img).clim_from_histogram(ratio=.9)
        fig = plt.figure()
        ax = fig.add_subplot(121)
        ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_title('Original')

    # create new Nansat with one band only
    return get_fixed_n(n02, img, time_coverage_start)

def get_img_bbox(geo):
    lon, lat = np.array(json.loads(geo.ExportToJson())['coordinates'][0]).T
    dom = Domain(NS_SRS, '-te 0 0 1000 1000 -tr 1 -1')
    x, y = dom.transform_points(lon, lat, 1)
    return x.min(), y.min(), x.max(), y.max()

def get_domain_extent(n1, n2, resolution=5e3, buffer=50e3):
    geo1 = n1.get_border_geometry()
    geo2 = n2.get_border_geometry()
    inter = geo1.Intersection(geo2)
    interlon, interlat = np.array(json.loads(inter.ExportToJson())['coordinates'][0]).T
    dom = Domain(NS_SRS, '-te 0 0 1000 1000 -tr 1 -1')

    xbrd, ybrd = dom.transform_points(interlon, interlat, 1)
    xmin = np.floor(xbrd.min()/resolution)*resolution - buffer
    xmax = np.ceil(xbrd.max()/resolution)*resolution + buffer
    ymin = np.floor(ybrd.min()/resolution)*resolution - buffer
    ymax = np.ceil(ybrd.max()/resolution)*resolution + buffer
    bbox = xmin, xmax, ymin, ymax

    # plot extent should contain both images
    xmin1, ymin1, xmax1, ymax1 = get_img_bbox(geo1)
    xmin2, ymin2, xmax2, ymax2 = get_img_bbox(geo2)
    bbox_plot = [
        np.min([xmin1, xmin2]) - buffer,
        np.min([ymin1, ymin2]) - buffer,
        np.max([xmax1, xmax2]) + buffer,
        np.max([ymax1, ymax2]) + buffer,
    ]
    return bbox, bbox_plot

def fake_feature_tracking(n1, n2):
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
    print(f'Created {np.sum(gpi)} fake features')
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
    if os.path.exists(fname) and not args.force:
        print(f'Loading {fname}')
        return dict(np.load(fname))

# image in stereographic projection
def get_projected_hh(n, dom):
    n.reproject(dom, addmask=True, resample_alg=0) # 0: NN; 1: Bilinear
    hh = np.ma.array(n[1], dtype=float, mask=1-n[2])
    n.undo()
    return hh

# watermask in specific projection
def get_projected_watermask(n, dom, **kwargs):
    n.reproject(dom, **kwargs)
    wmask = n.watermask()[1]
    n.undo()
    return wmask

def plot_ft(args, index, n1, n2, c1, r1, c2, r2):
    lon1b, lat1b = n1.get_border()
    lon2b, lat2b = n2.get_border()
    x1b, y1b = NS_PROJ(lon1b, lat1b)
    x2b, y2b = NS_PROJ(lon2b, lat2b)
    # convert row/column coordinates of matched features to lon/lat
    lon1ft, lat1ft = n1.transform_points(c1, r1)
    lon2ft, lat2ft = n2.transform_points(c2, r2)
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
    save_fig(args, fig, t, index, n1, n2)

def get_quiver_plot(n, d, d_extent, xq, yq, uq, vq, rpmq):
    xmin, xmax, ymin, ymax = d_extent
    (vmin,), (vmax,) = Figure(n[1]).clim_from_histogram(ratio=.9)
    hh = get_projected_hh(n, d)
    hh = np.ma.array(
            median_filter(hh.data, 3), mask=hh.mask)

    # plot valid vectors in stereographic projection
    fig = plt.figure(figsize=(20,20))
    ax = plt.axes(projection=NS_CRS)
    im = ax.imshow(hh, vmin=vmin, vmax=vmax, cmap='gray',
            origin='upper', zorder=0, extent=d_extent)
    # plt.colorbar(im, shrink=0.5)
    s = slice(None, None, 2) #show every 2nd vector
    quiv = ax.quiver(xq[s], yq[s], uq[s], vq[s], rpmq[s])#, scale=2)
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
    return fig

def plot_pm(args, index, n1, n2, bbox_plot,
        lon1pm, lat1pm, upm, vpm, apm, rpm, hpm, **kwargs):
    # compute ice drift speed [m/s]
    delta_t = (n2.time_coverage_start - n1.time_coverage_start).total_seconds()
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
    im = ax.imshow(u, vmin=-.25, vmax=.25)
    ax.contour(quality, levels=[args.quality_threshold])
    ax.set_title('x velocity, m/s')
    fig.colorbar(im, shrink=.4)
    t = Template('pm_quality/pm_quality_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1, n2)

    # plot final drift on image 1
    res = 1000 #m
    xmin, ymin, xmax, ymax = bbox_plot
    d = Domain(NS_SRS, '-te %f %f %f %f -tr %f %f' % (
        xmin, ymin, xmax, ymax, res, res))
    d_extent = [xmin, xmax, ymin, ymax]

    gpi = (quality > args.quality_threshold)
    xq, yq, uq, vq, rpmq = x1pm[gpi], y1pm[gpi], u[gpi], v[gpi], rpm[gpi]
    fig = get_quiver_plot(n1, d, d_extent, xq, yq, uq, vq, rpmq)
    t = Template('pm_quiver/pm_quiver_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1, n2)

    # plot final drift on image 2
    fig = get_quiver_plot(n2, d, d_extent, xq, yq, uq, vq, rpmq)
    t = Template('pm_quiver2/pm_quiver2_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1, n2)

    # deformation
    e1, e2, e3, tri_a, tri_p, triangles = get_deformation_nodes(xq, yq, uq, vq)
    # convert deformations to %/day
    e1 *= 8640000
    e2 *= 8640000
    e3 *= 8640000
    emax = 10# %/day
    for e,t,ttl,clim in zip(
            [e1,e2,e3],
            [Template('shear/shear_${dto1}-${dto2}_${index}.png'),
                Template('divergence/divergence_${dto1}-${dto2}_${index}.png'),
                Template('total-defor/total-defor_${dto1}-${dto2}_${index}.png'),
                ],
            ['Shear [%/day]', 'Divergence [%/day]', 'Total deformation [%/day]'],
            [(-emax,emax), (-emax,emax), (0,emax)],
            ):
        vmin, vmax = clim
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
        save_fig(args, fig, t, index, n1, n2)

def plot_pm_clean(args, index, n1, n2, bbox_plot,
        lon1pm, lat1pm, upm, vpm, apm, rpm, hpm,
        upm_clean, vpm_clean, **kwargs):
    # compute ice drift speed [m/s]
    delta_t = (n2.time_coverage_start - n1.time_coverage_start).total_seconds()
    u = upm / delta_t
    v = vpm / delta_t
    u2 = upm_clean / delta_t
    v2 = vpm_clean / delta_t

    # start points in stereoprojection
    x1pm, y1pm = NS_PROJ(lon1pm, lat1pm)

    fig = plt.figure()
    ax = fig.add_subplot(121)
    im = ax.imshow(u, vmin=-.25, vmax=.25)
    ax.set_title('x velocity, m/s')
    fig.colorbar(im, shrink=.4)

    ax = fig.add_subplot(122)
    im = ax.imshow(u2, vmin=-.25, vmax=.25)
    ax.set_title('Cleaned x velocity, m/s')
    fig.colorbar(im, shrink=.4)

    t = Template('clean-u/u_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1, n2)

    # plot final drift on image 1
    res = 1000 #m
    xmin, ymin, xmax, ymax = bbox_plot
    d = Domain(NS_SRS, '-te %f %f %f %f -tr %f %f' % (
        xmin, ymin, xmax, ymax, res, res))
    d_extent = [xmin, xmax, ymin, ymax]

    gpi = np.isfinite(u2*v2)
    xq, yq, uq, vq, rpmq = x1pm[gpi], y1pm[gpi], u[gpi], v[gpi], rpm[gpi]
    fig = get_quiver_plot(n1, d, d_extent, xq, yq, uq, vq, rpmq)
    t = Template('clean-quiver/clean-quiver_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1, n2)

    # plot final drift on image 2
    fig = get_quiver_plot(n2, d, d_extent, xq, yq, uq, vq, rpmq)
    t = Template('clean-quiver2/clean-quiver2_${dto1}-${dto2}_${index}.png')
    save_fig(args, fig, t, index, n1, n2)

    # deformation
    e1, e2, e3, tri_a, tri_p, triangles = get_deformation_nodes(x1pm[gpi], y1pm[gpi], u2[gpi], v2[gpi])
    # convert deformations to %/day
    e1 *= 8640000
    e2 *= 8640000
    e3 *= 8640000
    emax = 15 # %/day
    for e,t,ttl,clim in zip(
            [e1,e2,e3],
            [
                Template('clean-shear/shear_${dto1}-${dto2}_${index}.png'),
                Template('clean-divergence/divergence_${dto1}-${dto2}_${index}.png'),
                Template('clean-total-defor/total-defor_${dto1}-${dto2}_${index}.png'),
                ],
            ['Shear [%/day]', 'Divergence [%/day]', 'Total deformation [%/day]'],
            [(-emax,emax), (-emax,emax), (0,emax)],
            ):
        vmin, vmax = clim
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
        save_fig(args, fig, t, index, n1, n2)

def process_1pair(args, f1, f2, index):
    # create Nansat objects with one band only. 
    n1 = get_nansat(f1)
    n2 = get_nansat(f2)
    print(f'1st image start: {n1.time_coverage_start}')
    print(f'2nd image start: {n2.time_coverage_start}')

    # Fake feature tracking to get more starting points for pattern matching
    c1, r1, c2, r2 = fake_feature_tracking(n1, n2)
    plot_ft(args, index, n1, n2, c1, r1, c2, r2)

    # Domain for pattern matching
    res = 5e3
    (xmin, xmax, ymin, ymax), bbox_plot = get_domain_extent(
            n1, n2, resolution=res, buffer=25e3)
    dst_dom = Domain(NS_SRS,
        '-te %d %d %d %d -tr %d %d' %
        (xmin, ymin, xmax, ymax, res, res))
    lon1pm, lat1pm = dst_dom.get_geolocation_grids()
    print(lon1pm.shape)

    # Pattern matching
    # lon/lat grids for image_1
    t = Template('npz_files/pm_${dto1}-${dto2}_${index}.npz')
    pm_results = load_npz(args, t, index, n1, n2)
    if pm_results is None:
        print('Running pattern matching...')
        pm_results = dict(lon1pm=lon1pm, lat1pm=lat1pm,
                water_mask=get_projected_watermask(n1, dst_dom))
        # Run Pattern Matching for each element in lon1pm/lat1pm matrix
        # ice displacement upm and vpm are returned in meters in Stereographic projection
        (
                pm_results['upm'], pm_results['vpm'],
                pm_results['apm'], pm_results['rpm'], pm_results['hpm'],
                pm_results['lon2pm'], pm_results['lat2pm'],
                ) = pattern_matching(
                        lon1pm, lat1pm, n1, c1, r1, n2, c2, r2,
                        min_border=500, max_border=500,
                        srs=NS_SRS, img_size=71,
                        angles=[-10, -5, 0, 5, 10],
                        threads=10,
                        )
        pm_results['upm_clean'] = clean_velo_field(args,
                pm_results['upm'], pm_results['rpm'], pm_results['hpm'])
        pm_results['vpm_clean'] = clean_velo_field(args,
                pm_results['vpm'], pm_results['rpm'], pm_results['hpm'])
        save_npz(args, t, index, n1, n2, **pm_results)

    # count the number of drift vectors
    gpi = np.isfinite(pm_results['upm_clean'] * pm_results['vpm_clean'])
    print(f'Detected {np.sum(gpi)} good drift vectors')

    # diagnostic plots
    plot_pm      (args, index, n1, n2, bbox_plot, **pm_results)
    plot_pm_clean(args, index, n1, n2, bbox_plot, **pm_results)

    # Save images at higher resolution for later plotting
    t = Template('npz_files/hh_${dto1}-${dto2}_${index}.npz')
    fname = get_filename(args, t, index, n1, n2)
    if not os.path.exists(fname):
        print('Generating plot data...')
        res = 1e3
        ext = '-te %f %f %f %f -tr %f %f' %(xmin, ymin, xmax, ymax, res, res)
        dom = Domain(NS_SRS, ext)
        plot_data = dict()
        plot_data['lon_plot'], plot_data['lat_plot'] = dom.get_geolocation_grids()
        plot_data['hh1'] = get_projected_hh(n1, dom)
        plot_data['hh2'] = get_projected_hh(n2, dom)
        print(f'Saving {fname}')
        np.savez(fname, **plot_data)


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
                21)
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
