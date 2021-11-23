#!/usr/bin/env python
import os
import argparse
import numpy as np
from scipy.ndimage import convolve
from netCDF4 import Dataset


class SmoothWRF:
    """ smooth WRF netCDF file """

    def __init__(self, input_file, output_file, smoothing_factor=8):
        """
        Parameters:
        -----------
        input_file : str
            File to be smoothed
        output_file : str
            Smoothed file
        smoothing_factor : float
            Size of averaging kernel
        """
        self.input_file = input_file
        self.output_file = output_file
        self.smoothing_factor = smoothing_factor

    
    @classmethod
    def init_from_cli(cls):
        """
        Init new object from command line inputs

        Returns:
        --------
        obj : cls
            new instance of current class
        """
        ap = cls.get_parser()
        self = cls.__new__(cls)
        self.__init__(**vars(ap.parse_args()))
        return self


    @classmethod
    def get_parser(cls):
        """
        Create parser for command line inputs

        Returns
        -------
        parser : ArgumentParser
            parser with options: output_dir, proc_cfg, nextsim_cfg, date, duration
        """
        parser = argparse.ArgumentParser(description=cls.__doc__)
        parser.add_argument('input_file', type=str, help='File to be smoothed')
        parser.add_argument('output_file', type=str, help='Compressed file')
        parser.add_argument('-sf', '--smoothing-factor', default=8, type=int,
                help='size of averaging kernel')
        return parser


    def get_smoothed_field(self, fld):
        """
        smooth array with coarse graining (scipy.ndimage.convolve)

        Parameters:
        -----------
        fld : numpy.ndarray or numpy.ma.array
            shape (nt, ny, nx)

        Returns:
        --------
        fld_smooth : numpy.ndarray or numpy.ma.array
            - shape (nt, ny, nx)
            - same type as fld
        """
        n = self.smoothing_factor
        # kernel to average over a square of size n in the spatial dimensions
        # (does nothing in the temporal dimension)
        kernel = np.ones((1,n,n))/(n ** 2)
        ma = hasattr(fld, 'mask')
        if ma:
            fld_smooth = fld.filled(0.)
        else:
            fld_smooth = np.array(fld) #copy
        fld_smooth = convolve(fld_smooth, kernel, mode='nearest')
        if ma:
            fld_smooth = np.ma.array(fld_smooth, mask = fld.mask)
        return fld_smooth


    def make_new_netcdf_file(self, src_ds, dst_ds):
        '''
        Copy dataset but compress variables

        Parameters:
        -----------
        src_ds : pynextsimf.netcdf_io.NetcdfArcMFC
            Dataset for moorings file
        dst_ds : pynextsimf.netcdf_io.NetcdfArcMFC
            target dataset
        '''
        # copy global attributes
        dst_ds.setncatts(vars(src_ds))
        # copy dimensions
        for d in src_ds.dimensions.values():
            sz = d.size
            if d.isunlimited():
                sz = None
                time_name = d.name
            dst_ds.createDimension(d.name, sz)
        # copy/smooth variables
        for src_var in src_ds.variables.values():
            print(src_var.name)
            dst_var = dst_ds.createVariable(src_var.name, src_var.dtype,
                    src_var.dimensions, zlib=True) #zlib=True produces the compression
            dst_var.setncatts(vars(src_var))
            fld = src_var[:]
            if src_var.ndim == 3:
                # smooth time-dependent vars spatially
                fld  = self.get_smoothed_field(fld)
            dst_var[:] = fld


    def run(self):
        """ Compress input file """
        nc1 = self.input_file
        nc2 = self.output_file
        print('Size of %s: %i MB' %(nc1, os.stat(nc1).st_size/1e6))
        with Dataset(nc1, 'r') as src_ds, Dataset(nc2, 'w') as dst_ds:
            self.make_new_netcdf_file(src_ds, dst_ds)
        print('Size of %s: %i MB' %(nc2, os.stat(nc2).st_size/1e6))


if __name__ == '__main__':
    obj = SmoothWRF.init_from_cli()
    obj.run()
