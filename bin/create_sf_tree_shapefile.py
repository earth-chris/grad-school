# takes the csv file with the lat/lon coords for each tree in sf and
# creates a point-format shape file
#####

import csv
import ogr as ogr
import osr as osr
import aei as aei

# set up input file
infile = '/home/cba/cba/san_francisco/Street_Tree_Map.csv'

# set up output file
outfile = 'sf_tree_locations.shp'
driver = ogr.GetDriverByName("ESRI Shapefile")
data_source = driver.CreateDataSource(outfile)

# set the projection info
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)

# create the output layer
layer = data_source.CreateLayer("sf_tree_locations", srs, ogr.wkbPoint)

# add the necessary fields
treeID = ogr.FieldDefn("Tree ID", ogr.OFTInteger)
treeID.SetWidth(6)
species = ogr.FieldDefn("Sp. Name", ogr.OFTString)
species.SetWidth(30)
commonName = ogr.FieldDefn("CommonName", ogr.OFTString)
commonName.SetWidth(30)
plantType = ogr.FieldDefn("Plant Type", ogr.OFTString)
plantType.SetWidth(10)
careTaker = ogr.FieldDefn("Caretaker", ogr.OFTString)
careTaker.SetWidth(15)
DBH = ogr.FieldDefn("DBH", ogr.OFTInteger)
DBH.SetWidth(3)
lat = ogr.FieldDefn("Latitude", ogr.OFTReal)
lon = ogr.FieldDefn("Longitude", ogr.OFTReal)

# add to layer
layer.CreateField(treeID)
layer.CreateField(species)
layer.CreateField(commonName)
layer.CreateField(plantType)
layer.CreateField(careTaker)
layer.CreateField(DBH)
layer.CreateField(lat)
layer.CreateField(lon)

# open the file for reading
reader = csv.DictReader(open(infile, "rb"),
    delimiter = ',', quoting = csv.QUOTE_NONE)
    
# process the csv file and update fields per necessary
for row in reader:
    
    # move on if there is no lat/lon info
    if row['Latitude'] == '':
        continue
    if row['Longitude'] == '':
        continue
    
    # create the feature
    feature = ogr.Feature(layer.GetLayerDefn())
    
    # set attributes from csv
    feature.SetField("Tree ID", row["TreeID"])
    feature.SetField("Sp. Name", row["qSpecies"].split("::")[0].strip())
    feature.SetField("CommonName", row["qSpecies"].split("::")[1].strip())
    feature.SetField("Plant Type", row["PlantType"])
    feature.SetField("Caretaker", row["qCaretaker"])
    feature.SetField("DBH", row["DBH"])
    feature.SetField("Latitude", row["Latitude"])
    feature.SetField("Longitude", row["Longitude"])
    
    # set point geometry from WKT
    wkt = "POINT(%f %f)" %  (float(row['Longitude']) , float(row['Latitude']))
  
    # Create the point from the Well Known Txt
    point = ogr.CreateGeometryFromWkt(wkt)
    feature.SetGeometry(point)
    
    # Create the feature in the layer (shapefile)
    layer.CreateFeature(feature)
    
    # Destroy the feature to free resources
    feature.Destroy()
    
# Destroy the data source
data_source.Destroy()