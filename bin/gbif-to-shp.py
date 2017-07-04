import csv
import os as os
import ogr as ogr
import osr as osr
import aei as aei
import sys as sys

# create a class to parse out the arguments passed to the main function
class parse_args:
    def __init__(self, arglist):
        
        # set up main variables and defaults to parse
        self.infile = ''
        self.outfile = ''
        self.filterFlag = False
        self.filterType = ''
        self.filterName = ''
        
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
                      
            # check if filtering by kingdom            
            if arg.lower() == '-kingdom':
                i += 1
                arg = arglist[i]
                
                if not arg.lower() in ['animalia','plantae','fungi','protista','monera']:
                    usage()
                    print("[ ERROR ]: Unrecognized kingdom specified: %s" % arg)
                    sys.exit(1)
                    
                self.filterFlag = True
                self.filterType = 'kingdom'
                self.filterName = arg.lower().capitalize()
                
            # check if filtered by phyllum
            if arg.lower() == '-phylum':
                i += 1
                arg = arglist[i]
                
                # too many options to filter, just set variables
                self.filterFlag = True
                self.filterType = 'phylum'
                self.filterName = arg.lower().capitalize()
                
            # check if filtered by class
            if arg.lower() == '-class':
                i += 1
                arg = arglist[i]
                
                # too many options to filter, just set variables
                self.filterFlag = True
                self.filterType = 'class'
                self.filterName = arg.lower().capitalize()
                
            # check if filtered by order
            if arg.lower() == '-order':
                i += 1
                arg = arglist[i]
                
                # too many options to filter, just set variables
                self.filterFlag = True
                self.filterType = 'order'
                self.filterName = arg.lower().capitalize()
                
            # check if filtered by family
            if arg.lower() == '-family':
                i += 1
                arg = arglist[i]
                
                # too many options to filter, just set variables
                self.filterFlag = True
                self.filterType = 'family'
                self.filterName = arg.lower().capitalize()
                
            # check if filtered by genus
            if arg.lower() == '-genus':
                i += 1
                arg = arglist[i]
                
                # too many options to filter, just set variables
                self.filterFlag = True
                self.filterType = 'genus'
                self.filterName = arg.lower().capitalize()
                
            # check if filtered by country
            if arg.lower() == '-country':
                i += 1
                arg = arglist[i]
                
                # make sure country code is legit
                if not arg.upper() in get_country_codes():
                    usage()
                    print("[ ERROR ]: Unable to parse country code: %s" % arg)
                    print("[ ERROR ]: Must use 2-letter ISO country code")
                    sys.exit(1)
                    
                self.filterFlag = True
                self.filterType = 'countrycode'
                self.filterName = arg.upper()

            i += 1

# return list of the parameters set by gbif
def get_params():
    """
    returns the parameter keys for GBIF files
    usage: var = get_params()
    """
    # uncomment the fields you want to retain
    params = [
             'gbifid',
             'kingdom',
             'phylum',
             'class',
             'order',
             'family',
             'genus',
             'species',
             'infraspecificepithet',
             #'taxonrank',
             #'scientificname',
             'countrycode',
             'locality',
             #'publishingorgkey',
             'decimallatitude',
             'decimallongitude',
             #'coordinateuncertaintyinmeters',
             'day',
             'month',
             'year',
             'basisofrecord',
             #'institutioncode',
             #'collectioncode',
             #'catalognumber',
             #'recordnumber',
             'identifiedby',
             #'eventdate',
             #'rights',
             #'rightsholder',
             'recordedby',
             'typestatus',
             'establishmentmeans',
             #'mediatype',
             #'lastinterpreted',
             'issue'
             ]
        
    return params
    
def get_country_codes():
    """
    returns the country codes used in the GBIF database
    source: http://www.iso.org/iso/country_codes/iso_3166_code_lists/country_names_and_code_elements.htm
    """
    
    codes = [
        'AF',        'AX',        'AL',        'DZ',        'AS',        'AD',        'AO',        'AI',
        'AQ',        'AG',        'AR',        'AM',        'AW',        'AU',        'AT',        'AZ',
        'BS',        'BH',        'BD',        'BB',        'BY',        'BE',        'BZ',        'BJ',
        'BM',        'BT',        'BO',        'BQ',        'BA',        'BW',        'BV',        'BR',
        'IO',        'BN',        'BG',        'BF',        'BI',        'CV',        'KH',        'CM',
        'CA',        'KY',        'CF',        'TD',        'CL',        'CN',        'CX',        'CC',
        'CO',        'KM',        'CD',        'CG',        'CK',        'CR',        'CI',        'HR',
        'CU',        'CW',        'CY',        'CZ',        'DK',        'DJ',        'DM',        'DO',
        'EC',        'EG',        'SV',        'GQ',        'ER',        'EE',        'ET',        'FK',
        'FO',        'FJ',        'FI',        'FR',        'GF',        'PF',        'TF',        'GA',
        'GM',        'GE',        'DE',        'GH',        'GI',        'GR',        'GL',        'GD',
        'GP',        'GU',        'GT',        'GG',        'GN',        'GW',        'GY',        'HT',
        'HM',        'VA',        'HN',        'HK',        'HU',        'IS',        'IN',        'ID',
        'IR',        'IQ',        'IE',        'IM',        'IL',        'IT',        'JM',        'JP',
        'JE',        'JO',        'KZ',        'KE',        'KI',        'KP',        'KR',        'KW',
        'KG',        'LA',        'LV',        'LB',        'LS',        'LR',        'LY',        'LI',
        'LT',        'LU',        'MO',        'MK',        'MG',        'MW',        'MY',        'MV',
        'ML',        'MT',        'MH',        'MQ',        'MR',        'MU',        'YT',        'MX',
        'FM',        'MD',        'MC',        'MN',        'ME',        'MS',        'MA',        'MZ',
        'MM',        'NA',        'NR',        'NP',        'NL',        'NC',        'NZ',        'NI',
        'NE',        'NG',        'NU',        'NF',        'MP',        'NO',        'OM',        'PK',
        'PW',        'PS',        'PA',        'PG',        'PY',        'PE',        'PH',        'PN',
        'PL',        'PT',        'PR',        'QA',        'RE',        'RO',        'RU',        'RW',
        'BL',        'SH',        'KN',        'LC',        'MF',        'PM',        'VC',        'WS',
        'SM',        'ST',        'SA',        'SN',        'RS',        'SC',        'SL',        'SG',
        'SX',        'SK',        'SI',        'SB',        'SO',        'ZA',        'GS',        'SS',
        'ES',        'LK',        'SD',        'SR',        'SJ',        'SZ',        'SE',        'CH',
        'SY',        'TW',        'TJ',        'TZ',        'TH',        'TL',        'TG',        'TK',
        'TO',        'TT',        'TN',        'TR',        'TM',        'TC',        'TV',        'UG',
        'UA',        'AE',        'GB',        'UM',        'US',        'UY',        'UZ',        'VU',
        'VE',        'VN',        'VG',        'VI',        'WF',        'EH',        'YE',        'ZM',
        'ZW']
        
    return codes
        
# set a function to return how it works
def usage(exit=False):
    """
    describes the gbif2shp.py procedure in case of incorrect parameter calls
    
    syntax: usage()
    """
    print(
        """
$ gbif2shp.py [-i] input_file [-o] output_file [-country xx] [-kingdom xx]
    [-phylum xx] [-class xx] [-order xx] [-family xx] [-genus xx]
  where:
      input_file   = the input GBIF-format csv file
      output_file  = the output *.shp file
      country flag = the 2-letter ISO country code to filter by
      taxon flags  = filters the results to output only a specific taxon
      
      note: only one taxon flag will be read and filtered (i.e. you cannot 
            set -kingdom and -phylum simultaneously)
        """
        )
    if exit:
        sys.exit(1)

# main program
def main():
    """
    the main program for gbif2shp.py
    
    syntax: main()
    """
    # parse the argument list
    args = parse_args(sys.argv)
    #args = parse_args(argl)
    
    # check that an input file has been set
    if not aei.fn.checkFile(args.infile, quiet = True):
        usage()
        aei.fn.checkFile(args.infile)
        sys.exit(1)
        
    # get the GBIF header parameters
    params = get_params()
    
    # set up the output driver
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.CreateDataSource(args.outfile)
    
    # set the projection info
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    
    # create the output layer
    layer = data_source.CreateLayer("GBIF species occurrence", srs, ogr.wkbPoint)
    
    # add the necessary fields to the DBF
    for field in params:
        
        # set up integer, float or string by param
        if field in ['gbfid','day','month','year']:
            field_id = ogr.FieldDefn(field, ogr.OFTInteger)
            
        elif field in ['decimallatitude','decimallongitude']:
            field_id = ogr.FieldDefn(field, ogr.OFTReal)
            
        else:
            field_id = ogr.FieldDefn(field, ogr.OFTString)
            
        layer.CreateField(field_id)
    
    # open the file for reading
    reader = csv.DictReader(open(args.infile, "rb"),
        delimiter = '\t', quoting = csv.QUOTE_NONE)
    
    # debug
    print(args.filterFlag)
    print(args.filterType)
    print(args.filterName)
        
    # process the csv file and update fields per necessary
    j = long(0)
    for row in reader:
        
        # report on the status every 100,000 lines
        if j in range(0,long(1e8),100000l): 
            print("[ STATUS ]: Processing line %s") % (j+1)
        
        # iterate counter
        j += 1l
        
        # move on if there is no lat/lon info
        if row['decimallatitude'] == '':
            continue
        if row['decimallongitude'] == '':
            continue
        
        # move on if there are any filtering operations
        if args.filterFlag:
            if row[args.filterType] != args.filterName:
                continue
        
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        
        # set attributes from csv
        for field in params:
            feature.SetField(field, row[field])
        
        # set point geometry from WKT
        wkt = "POINT(%f %f)" %  (float(row['decimallongitude']) , float(row['decimallatitude']))
         
        # Create the point from the well known text
        point = ogr.CreateGeometryFromWkt(wkt)
        feature.SetGeometry(point)
        
        # Create the feature in the layer (shapefile)
        layer.CreateFeature(feature)
        
        # Destroy the feature to free resources
        feature.Destroy()
        
    print(j)
    # Destroy the data source
    data_source.Destroy()
    
if __name__ == "__main__":
    main()    