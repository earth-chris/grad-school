#/usr/bin/python!
#####
# unmix_spectra.py performs fully-constrained least-squares 
#  spectral ummixing on an image file
#
# unmixing implementation described here:
#  http://pysptools.sourceforge.net/abundance_maps.html
#
# based on the algorithm desribed here:
#  Daniel Heinz, Chein-I Chang, and Mark L.G. Fully Constrained 
#  Least-Squares Based Linear Unmixing. Althouse. IEEE. 1999.
#
# c. 2016 Christopher Anderson
#####

import os
import sys
import aei
import random
import gdal as gdal
import numpy as np
import pysptools.abundance_maps as abm

# create a class to parse out the arguments passed to the main function
class parse_args:
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
                    if not aei.fn.checkFile(self.infile, quiet = True):
                        usage()
                        aei.fn.checkFile(self.infile)
                        sys.exit(1)
            
            # check output flag
            if arg.lower() == '-o':
                i += 1
                arg = arglist[i]
                
                if type(arg) is str:
                    self.outfile = arg
                    outpath = os.path.dirname(self.outfile)
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
                libs = arg.split(" ")
                
                # throw an error if only one lib specified
                if len(libs) == 1:
                    usage()
                    print("[ ERROR ]: unable to unmix with one spectral library: %s" % libs[0])
                    sys.exit(1)
                
                # loop through each lib and update spec_lib list
                for j in range(len(libs)):
                    if not aei.fn.checkFile(libs[j], quiet=True):
                        usage()
                        aei.fn.checkFile(libs[j])
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
                band = arg.split(" ")
                
                # loop through and make sure each is a number
                for j in range(len(band)):
                    
                    try:
                        int(band[j])
                    except ValueError:
                        usage()
                        print("[ ERROR ]: invalid index set: %s" % ind[j])
                        sys.exit(1)
                    
                    self.bands.append(int(band[j]))
                    
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
            
            i += 1

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
    n_libs = len(args.spectral_libs)
    lnb = np.zeros(n_libs)
    
    # get info from each library
    for i in range(n_libs):
        
        # assign libraries to dictionary
        lib['lib_%s' % i] = aei.read.spectralLib(args.spectral_libs[i])
        lnb[i] = (lib['lib_%s' % i].spectra.shape[-1])
        
        # get random indices for each library
        lib['rnd_%s' % i] = random.sample(range(
            lib['lib_%s' % i].spectra.shape[0]), args.n)
            
        # normalize if set
        if args.normalize:
            lib['lib_%s' % i].bn(inds=args.bands)
    
    # load the input image and get parameters
    inf = gdal.Open(args.infile)
    ns = inf.RasterXSize
    nl = inf.RasterYSize
    nb = inf.RasterCount
    geo = inf.GetGeoTransform()
    prj = inf.GetProjection()
    n_bundles = lnb.shape[0]
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
    
    # report a little
    print("[ STATUS ]: Beginning fully-constrained least-squares unmixing")
    print("[ STATUS ]: Input file : %s" % args.infile)
    print("[ STATUS ]: Output file: %s" % args.outfile)
        
    # create an output file
    ouf = gdal.GetDriverByName(args.of).Create(
              args.outfile, ns, nl, n_bundles, gdal.GDT_Float32)
    ouf.SetGeoTransform(geo)
    ouf.SetProjection(prj)
    
    # create an output array 
    arr = np.zeros((nl, ns, n_bundles))
    
    # read the image file and flip dimensions
    img = inf.ReadAsArray()
    img = img.transpose([1,2,0])
    
    # find no data vals
    gd = np.where(img[:,:,0] != nd)
    
    # check that the file is not all no-data
    if gd[0].shape[0] == 0:
        print("[ ERROR ]: No good-data found in input file")
        print("[ ERROR ]: Exiting...")
        sys.exit(1)
    
    # subset the image array to good data only
    img = img[gd[0], gd[1], :]
    
    # normalize the data if set
    if args.normalize:
        img = aei.fn.bn(img, inds = args.bands)
        args.bands = range(len(args.bands))
        
    # add a shallow dimension for unmixing algorithm
    img = np.expand_dims(img[:,args.bands], 0)
    
    # set up unmixing class
    unmixer = abm.FCLS()
    
    # loop through each random index and unmix
    for i in range(args.n):
        
        print("[ STATUS ]: Iteration [%s] of [%s]" % ((i + 1), args.n))
        
        # set up the bundles 
        bundles = np.zeros((n_libs,len(args.bands)))
        for j in range(n_libs):
            bundles[j,:] = lib['lib_%s' % j].spectra[lib['rnd_%s' % j][i],args.bands]
            
        # perform the unmixing
        arr[gd[0],gd[1],:] += unmixer.map(img, bundles).squeeze()
        
    # divide by n iterations to get average response
    arr /= args.n
    
    # report completion
    print("[ STATUS ]: Completed unmixing iterations")
    print("[ STATUS ]: Writing data to file")
    
    # and write the results to the output file
    for i in range(n_libs):
        
        band = ouf.GetRasterBand(i+1)
        band.WriteArray(arr[:,:,i])
        band.SetNoDataValue(0.0)
        band.FlushCache()
                
if __name__ == "__main__":
    main()