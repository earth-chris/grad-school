#####
# creates raster files with unique values for each unique attribute
#  in an input shape file
#
# c. 2016 Christopher Anderson
#####

import aei
import gdal as gdal
import ogr as ogr
import osr as osr

# set up errr handling
ogr.UseExceptions()
gdal.UseExceptions()

class parse_args:
    def __init__(self, arglist):
    
        # set up main variables and defaults to parse
        self.infile = ''
        self.outfile = ''
        self.te = []
        self.tr = []
        self.prj = ''
        
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
                    if not aei.checkFile(self.infile, quiet = True):
                        usage()
                        aei.checkFile(self.infile)
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
        
def usage(exit=False):
    """
    describes the rasterize_by_attribute.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ rasterize_by_attribute.py [-tr xres yres] [-te xmin ymin xmax ymax] 
      [-a attribute] [-use_values] [-tf target_file]
      [-i] input_file [-o] output_file
        """
        )
    if exit:
        sys.exit(1)

def main():
    """
    the main program for rasterize_by_attribute.py
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    
    # open the input file read-only and get info
    shape = ogr.Open(args.infile, 0)
    layer = shape.GetLayer(0)
    nFeatures = layer.GetFeatureCount()
    fields = []
    
    # test that attribute set actually works. if not, exit
    feature = layer.GetFeature(0)
    
    try:
        dummy = feature.GetField(args.a)
    except Exception, e:
        usage()
        print("[ ERROR ]: unable to parse -a attribute: %s" % args.a)
        print("[ ERROR ]: %s" % e)
        sys.exit(1)
    
    # loop through the vector and get the attribute info for each feature
    for i in range(nFeatures):
        feature = layer.GetFeature(i)
        fields.append(feature.GetField(args.a))
        
    # re-cast as a set to isolate only unique values, then re-cast as list
    fields = list(set(fields)).sort()
    
    # get count of unique items / number of loops to iterate
    nFields = len(fields)
    
    # create array of unique values to burn into raster. default
    #  is to just count from 1 to nFields, but can assign the attribute
    #  value if it contains a series of numbers
    if args.use_values:
        
        burnVals = np.asarray(fields)
        
        if not np.isreal(burnVals):
            usage()
            print("[ ERROR ]: unable to use -use_values")
            print("[ ERROR ]: the attribute set (%s) does not contain numerical values" % args.a)
            sys.exit(1)
            
    else:
        burnVals = np.arange(nFields)+1
        
    # loop through each field and burn it in to the output raster
    for i in range(nFields):
        
        # stuff happens here
        dummy = True
    
if __name__ == "__main__":
    main()