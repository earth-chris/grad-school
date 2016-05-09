#####
# this script uses PyPROSAIL and Py6S to model canopy reflectance
# for a suite of vegetation types for use in unmixing models and in
# canopy reflectance inversions
#
# c. 2016 Christopher Anderson
#####

import numpy as np
import pyprosail
from Py6S import *

# load sixs
s = SixS()

#####
# set up output files and processing parameters
#####

# set output file name
output_file = 'tropical_atmosphere_landsat8_veg_bundles.csv'

# set number of random iterations to produce
n_iterations = 10

#####
# set up the spectral output parameters
#####

# specify the output sensor configuration for modeled spectra
#  options include ali, aster, er2_mas, gli, landsat_etm, landsat_mss, landsat_oli,
#  landsat_tm, meris, modis, polder, spot_hrv, spot_vgt, vnir, whole_range, and custom
target_sensor = 'landsat_oli'

# set up output wavelengths (in um) if using custom wl range (ignored for pre-defined sensors)
wl_start = 0.4
wl_end = 2.5
wl_interval = 0.01
wl = np.arange(wl_start, wl_end, wl_interval)
wl_sixs = Wavelength(wl_start, wl_end)

# set up output type (an option from s.outputs)
#  examples include pixel_radiance, pixel_reflectance, apparent_radiance, apparent_reflectance, etc.
output_type = 'apparent_reflectance'

#####
# set up atmospheric modeling parameters
#####

# select the atmospheric profile (an option from AtmosProfile)
#  examples include Tropical, MidlatitudeSummer, MidlatitudeWinter, etc.
atmos_profile = AtmosProfile.Tropical

# select the aerosol profile to use (an option from AeroProfile)
#  examples include BiomassBurning, Continental, Desert, Maritime, NoAerosols, Urban, etc.
aero_profile = AeroProfile.Continental

# select ground reflectance method (an option from GroundReflectance)
#  examples include GreenVegetation, Hetero(/Homo)geneousLambertian, HotSpot, LeafDistPlanophile, etc
ground_reflectance = GroundReflectance.HomogeneousLambertian

# set altitudes for sensor (ground) and target (target altitude for modeled reflectance)
#  Landsat 4, 5, 7, 8 altitude - 705 km 
altitudes = Altitudes()
#altitudes.set_sensor_custom_altitude(705)
altitudes.set_sensor_satellite_level()
altitudes.set_target_sea_level()

# set aerosol optical thickness (550 nm)
aot = 0.5

# select viewing geometry parameters
geo = Geometry.User()
geo.solar_a = 0
geo.solar_z = 30
geo.view_a = 0
geo.view_z = 10

#####
# set up the leaf and canopy modeling parameters
#####

# structural coefficient (arbitrary units)
#  range 1.3 - 2.5 from Rivera et al. 2013 http://dx.doi.org/10.3390/rs5073280
N = 1.5

# total chlorophyll content (ug/cm^2)
#  range ~ 5 - 75 from Rivera et al. 2013
chloro = 40

# total carotenoid content (ug/cm^2)
caroten = 8

# brown pigment content (arbitrary units)
brown = 0

# equivalent water thickness (cm)
#  range 0.002 - 0.05 from Rivera et al. 2013
EWT = 0.01

# leaf mass per area (g/cm^2)
#  global range 0.0022 - 0.0365 (median 0.01)
#  from Asner et al. 2011 http://dx.doi.org/10.1016/j.rse.2011.08.020
LMA = 0.01

# soil reflectance metric (wet soil = 0, dry soil = 1)
soil_reflectance = 1

# leaf area index (unitless, cm^2 leaf area/cm^2 ground area)
#  range 0.01 - 18.0 (5.5 mean) globally
#  range 0.2 - 8.7 (3.6 mean) for crops
#  range 0.6 - 2.8 (1.3 mean) for desert plants
#  range 0.5 - 6.2 (2.6 mean) for boreal broadleaf forest
#  range 0.5 - 8.5 (4.6 mean) for boreal needle forest
#  range 0.8 - 11.6 (5.1 mean) for temperate broadleaf forest
#  range 0.01 - 15.0 (5.5 mean) for temperate needle forest
#  range 0.6 - 8.0 (4.8 mean) for tropical broadleaf forest
#  range 0.3 - 5.0 (1.7 mean) for grasslands
#  range 1.6 - 18.0 (8.7 mean) for plantations
#  range 0.4 - 4.5 (2.1 mean) for shrublands
#  range 0.2 - 5.3 (1.9 mean) for tundra
#  range 2.5 - 8.4 (6.3 mean) for wetlands
#  from Asner, Scurlock and Hicke 2003 http://dx.doi.org/10.1046/j.1466-822X.2003.00026.x
LAI = 3

# hot spot parameter (derived from brdf model)
#  range 0.05-0.5 from Rivera et al. 2013
hot_spot = 0.01

# solar zenith angle (degrees)
solar_zenith = geo.solar_z

# solar_azmiuth angle (degrees)
solar_azimuth = geo.solar_a

# view zenith angle (degrees)
view_zenith = geo.view_z

# view azimuth angle (degrees)
view_azimuth = geo.view_a

# leaf distribution function parameter. options include:
# Planophile, Erectophile, Plagiophile, Extremophile, Spherical, Uniform
LIDF = pyprosail.Planophile

#####
# set up if statements to set parameters based on sesnsor type
#####
if target_sensor == 'custom':
    run_sixs_params = SixSHelpers.Wavelengths.run_wavelengths
    num_bands = len(wl)
    good_bands = np.arange(num_bands)
elif target_sensor == 'ali':
    run_sixs_params = SixSHelpers.Wavelengths.run_ali
    num_bands = 7
    good_bands = np.arange(num_bands)
elif target_sensor == 'aster':
    run_sixs_params = SixSHelpers.Wavelengths.run_aster
    num_bands = 9
    good_bands = np.arange(num_bands)
elif target_sensor == 'er2_mas':
    run_sixs_params = SixSHelpers.Wavelengths.run_er2_mas
    num_bands = 7
    good_bands = np.arange(num_bands)
elif target_sensor == 'gli':
    run_sixs_params = SixSHelpers.Wavelengths.run_gli
    num_bands = 30
    good_bands = np.arange(num_bands)
elif target_sensor == 'landsat_etm':
    run_sixs_params = SixSHelpers.Wavelengths.run_landsat_etm
    num_bands = 6
    good_bands = np.arange(num_bands)
elif target_sensor == 'landsat_mss':
    run_sixs_params = SixSHelpers.Wavelengths.run_landsat_mss
    num_bands = 4
    good_bands = np.arange(num_bands)
elif target_sensor == 'landsat_oli':
    run_sixs_params = SixSHelpers.Wavelengths.run_landsat_oli
    num_bands = 10
    # use only the 6 traditional optical bands
    good_bands = np.arange(6)+1
elif target_sensor == 'landsat_tm':
    run_sixs_params = SixSHelpers.Wavelengths.run_landsat_tm
    num_bands = 6
    good_bands = np.arange(num_bands)
elif target_sensor == 'meris':
    run_sixs_params = SixSHelpers.Wavelengths.run_meris
    num_bands = 16
    good_bands = np.arange(num_bands)
elif target_sensor == 'modis':
    run_sixs_params = SixSHelpers.Wavelengths.run_modis
    num_bands = 8
    good_bands = np.arange(num_bands)
elif target_sensor == 'polder':
    run_sixs_params = SixSHelpers.Wavelengths.run_polder
    num_bands = 8
    good_bands = np.arange(num_bands)
elif target_sensor == 'spot_hrv':
    run_sixs_params = SixSHelpers.Wavelengths.run_spot_hrv
    num_bands = 4
    good_bands = np.arange(num_bands)
elif target_sensor == 'spot_vgt':
    run_sixs_params = SixSHelpers.Wavelengths.run_spot_vgt
    num_bands = 4
    good_bands = np.arange(num_bands)
elif target_sensor == 'vnir':
    run_sixs_params = SixSHelpers.Wavelengths.run_vnir
    num_bands = 200
    good_bands = np.arange(num_bands)
elif target_sensor == 'whole_range':
    run_sixs_params = SixSHelpers.Wavelengths.run_whole_range
    num_bands = 380
    # use only 400-2500 nm range
    good_bands = np.arange(210)+20
else:
    raise OSError('Unsupported sensor configuration')

#####
# load prosail and run the canopy model
#####

spectrum = pyprosail.run(N, chloro, caroten, brown, EWT, LMA, soil_reflectance,
			            LAI, hot_spot, solar_zenith, solar_azimuth,
			            view_zenith, view_azimuth, LIDF)

# set the ground reflectance component to the canopy model
s.ground_reflectance = ground_reflectance(spectrum)

#####
# update the sixs parameters for running
#####

# update aerosol, altitude, aot, atmosphere, geometry and wavelength info
s.aero_profile = AeroProfile.PredefinedType(aero_profile)
s.altitudes = altitudes
s.aot550 = aot
s.atmos_profile = AtmosProfile.PredefinedType(atmos_profile)
s.geometry = geo

#####
# now that our parameters are set, run the model for each iteration
#####

# generate the output spectrum using specified wavelengths
if target_sensor == 'custom':
    wavelengths, results = run_sixs_params(s, wl, output_name = output_type)
else:
    wavelengths, results = run_sixs_params(s, output_name = output_type)

#convert output to array for ease of output
wavelengths = np.asarray(wavelengths)
results = np.asarray(results)

# plot the result
SixSHelpers.Wavelengths.plot_wavelengths(wavelengths[good_bands], results[good_bands], output_type)
