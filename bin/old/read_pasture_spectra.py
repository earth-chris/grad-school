import os
import csv
import pickle
import numpy as np

# set the paths for the input/output files
spectral_file = "/home/cba/cba/global/spectra/ND01_PASTURE_SPECTRA_1154/data/AllSpectra_Values.csv"
biophysical_file = "/home/cba/cba/global/spectra/ND01_PASTURE_SPECTRA_1154/data/Pasture_spectra_biophysical_measurements.csv"
output_file = os.path.join(os.environ['AEI_GS'],'data','ND01_pasture_spectra.pickle')

# get number of lines in each file
slines = sum(1 for line in open(spectral_file))
blines = sum(1 for line in open(biophysical_file))

# state number of lines to ignore
signore = 12
bignore = 15

# read in data as csv reader objects
sfile = open(spectral_file, "r")
bfile = open(biophysical_file, "r")
spectral_csv = csv.reader(sfile)
biophysical_csv = csv.reader(bfile)

# get through the metadata for each file
#  spectral file has 11 lines of metadata to exclude
for i in range(0,signore-1):
    sdummy = spectral_csv.next()

#  biophysical file has 15 lines of metadata to exclude
for i in range(0,bignore):
    bdummy = biophysical_csv.next()
    
# get wavelengths from next line of spectral cv file
wavelengths = spectral_csv.next()[2:]
                                     
# set number of records to read
#n_recs = 68

# set up variables to read in as arrays
reflectance = {
    'point_id' : [],
    'biophysical' : [],
    'spectra' : np.zeros((slines-signore, len(wavelengths)))
    }
    
metadata = {
    'point_id' : [],
    'species' : [],
    'biomass_live' : [],
    'biomass_senesced' : [],
    'biomass_total' : [],
    'water_content' : []
    }
    
# read each line of the csv file and assign variables
#for i in range(0,n_recs-1):
#    sline = spectral_csv.next()
#    bline = biophysical_csv.next()
#    reflectance[i] = sline[2:]
#    metadata['point_id'].append(bline[0])
#    metadata['species'].append(bline[3])
#    metadata['biomass_live'].append(bline[4])
#    metadata['biomass_senesced'].append(bline[5])
#    metadata['biomass_total'].append(bline[6])
#    metadata['water_content'].append(bline[7])

i = 0    
for row in spectral_csv:
    reflectance['point_id'].append(row[0])
    reflectance['biophysical'].append(row[1])
    reflectance['spectra'][i] = row[2:]
    i += 1
    
for row in biophysical_csv:
    metadata['point_id'].append(row[0])
    metadata['species'].append(row[3])
    metadata['biomass_live'].append(row[4])
    metadata['biomass_senesced'].append(row[5])
    metadata['biomass_total'].append(row[6])
    metadata['water_content'].append(row[7])
    
# close files after reading
sfile.close()
bfile.close()

# write output file
with open(output_file, 'w') as f:
    pickle.dump([reflectance,metadata], f)