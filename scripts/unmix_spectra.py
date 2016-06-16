#####
# unmix_spectra.py performs spectral ummixing on an image file
#####
import os
import sys
import aei
import gdal as gdal
import numpy as np
import pysptools.abundance_maps as abm

def usage(exit=False):
    """
    describes the aeilas.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ unmix_spectra.py -lib "lib1 lib2 ... libx" [-n n_random_selections]
      [-inds "ind1 ind2 ... indx"] [-normalize] [-of output_format]
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
    
    # set up main variables to parse
    infile = ''
    outfile = ''
    spectral_libs = []
    n = 20
    inds = []
    of = 'GTiff'
    normalize = False
    
    if len(sys.argv) == 1:
        usage(exit=True)

    # read arguments from command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
    
        # check input flag    
        if arg.lower() == '-i':
            i += 1
            arg = sys.argv[i]
            
            if type(arg) is str:
                infile = arg
                if not aei.checkFile(infile, quiet = True):
                    usage()
                    aei.checkFile(infile)
                    sys.exit(1)
        
        # check output flag
        if arg.lower() == '-o':
            i += 1
            arg = sys.argv[i]
            
            if type(arg) is str:
                outfile = arg
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
            arg = sys.argv[i]
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
                    
                spectral_libs.append(libs[j])
                
        # check number of iterations
        if arg.lower() == "-n":
            i += 1
            arg = sys.argv[i]
            
            try:
                n = int(arg)
            except ValueError:
                usage()
                print("[ ERROR ]: -n argument is not an integer: %s" % arg)
                sys.exit(1)
                
        # check indices
        if arg.lower() == "-inds":
            i += 1
            arg = sys.argv[i]
            ind = arg.split(", ")
            
            # loop through and make sure each is a number
            for j in range(len(ind)):
                
                try:
                    int(ind[j])
                except ValueError:
                    usage()
                    print("[ ERROR ]: invalid index set: %s" % ind[j])
                    sys.exit(1)
                
                inds.append(ind[j])
                
        # check normalize flag
        if arg.lower() == "-normalize":
            normalize = True
            
        # check output format
        if arg.lower() == "-of":
            i += 1
            arg = sys.argv[i]
            
            import gdal_types as gdt
            if not arg in gdt.list():
                usage()
                print("[ ERROR }: invalid output format: %s" % arg)
                sys.exit(1)
                
            of = arg
                
if __name__ == "__main__":
    main()