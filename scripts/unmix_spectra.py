#####
# unmix_spectra.py performs spectral ummixing on an image file
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import gdal as gdal
import numpy as np
import pysptools.abundance_maps as abm

# create a class to parse out the arguments passed to the main function
class parse_args(arglist):
    
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.infile = ''
        self.outfile = ''
        self.spectral_libs = []
        self.n = 20
        self.bands = []
        self.of = 'GTiff'
        self.normalize = False
        
        # exit if no arguments passed
        if len(arglist) == 1:
            usage(exit=True)
    
        # read arguments from command line
        i = 1
        while i < len(arglist):
            arg = arglist[i]
        
            # check input flag    
            if arg.lower() == '-i':
                i += 1
                arg = arglist[i]
                
                if type(arg) is str:
                    self.infile = arg
                    if not aei.checkFile(infile, quiet = True):
                        usage()
                        aei.checkFile(infile)
                        sys.exit(1)
            
            # check output flag
            if arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                if type(arg) is str:
                    self.outfile = arg
                    outpath = os.path.dirname(outfile)
                    if outpath == '':
                        outpath = '.'
                    if not os.access(outpath, os.W_OK):
                        usage()
                        print("[ ERROR ]: unable to write to output path: %s" % outpath)
                        sys.exit(1)
            
            # check spectral library paths            
            if arg.lower() == "-lib":
                i += 1
                arg = arglist[i]
                libs = arg.split(", ")
                
                # throw an error if only one lib specified
                if len(libs) == 1:
                    usage()
                    print("[ ERROR ]: unable to unmix with one spectral library: %s" % libs[0])
                    sys.exit(1)
                
                # loop through each lib and update spec_lib list
                for j in range(len(libs)):
                    if not aei.checkFile(libs[j], quiet=True):
                        usage()
                        aei.checkFile(libs[j])
                        sys.exit()
                        
                    self.spectral_libs.append(libs[j])
                    
            # check number of iterations
            if arg.lower() == "-n":
                i += 1
                arg = arglist[i]
                
                try:
                    self.n = int(arg)
                except ValueError:
                    usage()
                    print("[ ERROR ]: -n argument is not an integer: %s" % arg)
                    sys.exit(1)
                    
            # check indices
            if arg.lower() == "-bands":
                i += 1
                arg = arglist[i]
                band = arg.split(", ")
                
                # loop through and make sure each is a number
                for j in range(len(band)):
                    
                    try:
                        int(band[j])
                    except ValueError:
                        usage()
                        print("[ ERROR ]: invalid index set: %s" % ind[j])
                        sys.exit(1)
                    
                    self.bands.append(band[j])
                    
            # check normalize flag
            if arg.lower() == "-normalize":
                self.normalize = True
                
            # check output format
            if arg.lower() == "-of":
                i += 1
                arg = arglist[i]
                
                import gdal_types as gdt
                if not arg in gdt.list():
                    usage()
                    print("[ ERROR }: invalid output format: %s" % arg)
                    sys.exit(1)
                    
                self.of = arg

def usage(exit=False):
    """
    describes the aeilas.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ unmix_spectra.py -lib "lib1 lib2 ... libx" [-n n_random_selections]
      [-bands "band1 band2 ... bandx"] [-normalize] [-of output_format]
      [-i] input_file [-o] output_file
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for unmix_spectra.py
    
    syntax: main()
    """
    
    # parse the argument list
    args = parse_args(sys.argv)
    
    # load the spectral libraries
    lib = {}
    lnb = np.zeros(len(args.spectral_libs))
    
    for i in range(len(args.spectral_libs)):
        lib['%s' % i] = aei.readSpeclib(args.spectral_libs)
        lnb[i] = (lib['%s' % i].spectra.shape[-1])
    
    # load the input image and get parameters
    inf = gdal.Open(args.infile)
    ns = inf.RasterXSize
    nl = inf.RasterYSize
    nb = inf.RasterCount
    geo = inf.GetGeoTransform()
    prj = inf.GetProjection()
    n_bundles = lnb.shape
    b1 = inf.GetRasterBand(1)
    
    # if bands were not set in command line, use all bands
    if not args.bands:
        args.bands = range(nb)
    
    # check no data param
    if hasattr(b1, 'GetNoDataValue'):
        nd = b1.GetNoDataValue()
        if nd is None:
            nd = -9999
    
    else:
        nd = -9999
    
    # check that the spectral libraries match the bands
    if not 0 in np.where((lnb-nb) != 0)[0].shape:
        print("[ ERROR ]: number of image bands does not match number of spectral library bands")
        inf = None
        sys.exit(1)
        
    # create an output file
    ouf = gdal.GetDriverByName(args.of).Create(
              args.outfile, ns, nl, n_bundles, gdal.GDT_Float32)
    ouf.SetGeoTransform(geo)
    ouf.SetProjection(prj)
    
    # create an output array 
    arr = np.empty((nl, ns, n_bundles))
    
    # read the image file and flip dimensions
    img = inf.ReadAsArray()
    img = img.transpose([1,2,0])
    
    # find no data vals
    gd = np.where(img[:,:,0] != nd)[0]
    
    # set up unmixing class
    unmixer = abm.FCLS()
    
    # get random indices for each bundle
    
    # loop through each random index and unmix
    for i in range(args.n):
        
                
if __name__ == "__main__":
    main()