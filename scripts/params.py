def get_sensor_bands(target_sensor):
    """
    get_sensor_bands returns a list of the band definitions from
    the PredfinedWavelengths function for a list of different sensors
    
    availalbe definitions for target_sensor are:
        custom, ali, aster, er2_mas, gli, landsat_etm, *_mss, *_oli, *_tm,
        meris, modis, polder, spot_hrv, spot_vgt
    """
    if target_sensor == 'custom':
        run_sixs_params = SixSHelpers.Wavelengths.run_wavelengths
    elif target_sensor == 'ali':
        run_sixs_params = SixSHelpers.Wavelengths.run_ali
        bands = [PredefinedWavelengths.ALI_B1,
                 PredefinedWavelengths.ALI_B2,
                 PredefinedWavelengths.ALI_B3,
                 PredefinedWavelengths.ALI_B4,
                 PredefinedWavelengths.ALI_B5,
                 PredefinedWavelengths.ALI_B7]
    elif target_sensor == 'aster':
        run_sixs_params = SixSHelpers.Wavelengths.run_aster
        bands = [PredefinedWavelengths.ASTER_B1,
                 PredefinedWavelengths.ASTER_B2,
                 PredefinedWavelengths.ASTER_B3B,
                 PredefinedWavelengths.ASTER_B4,
                 PredefinedWavelengths.ASTER_B5,
                 PredefinedWavelengths.ASTER_B6,
                 PredefinedWavelengths.ASTER_B7,
                 PredefinedWavelengths.ASTER_B8,
                 PredefinedWavelengths.ASTER_B9]
    elif target_sensor == 'er2_mas':
        run_sixs_params = SixSHelpers.Wavelengths.run_er2_mas
        bands = [PredefinedWavelengths.ER2_MAS_B1,
                 PredefinedWavelengths.ER2_MAS_B2,
                 PredefinedWavelengths.ER2_MAS_B3,
                 PredefinedWavelengths.ER2_MAS_B4,
                 PredefinedWavelengths.ER2_MAS_B5,
                 PredefinedWavelengths.ER2_MAS_B6,
                 PredefinedWavelengths.ER2_MAS_B7]
    elif target_sensor == 'gli':
        run_sixs_params = SixSHelpers.Wavelengths.run_gli
        bands = [PredefinedWavelengths.GLI_B1,
                 PredefinedWavelengths.GLI_B2,
                 PredefinedWavelengths.GLI_B3,
                 PredefinedWavelengths.GLI_B4,
                 PredefinedWavelengths.GLI_B5,
                 PredefinedWavelengths.GLI_B6,
                 PredefinedWavelengths.GLI_B7,
                 PredefinedWavelengths.GLI_B8,
                 PredefinedWavelengths.GLI_B9,
                 PredefinedWavelengths.GLI_B10,
                 PredefinedWavelengths.GLI_B11,
                 PredefinedWavelengths.GLI_B12,
                 PredefinedWavelengths.GLI_B13,
                 PredefinedWavelengths.GLI_B14,
                 PredefinedWavelengths.GLI_B15,
                 PredefinedWavelengths.GLI_B16,
                 PredefinedWavelengths.GLI_B17,
                 PredefinedWavelengths.GLI_B18,
                 PredefinedWavelengths.GLI_B19,
                 PredefinedWavelengths.GLI_B20,
                 PredefinedWavelengths.GLI_B21,
                 PredefinedWavelengths.GLI_B22,
                 PredefinedWavelengths.GLI_B23,
                 PredefinedWavelengths.GLI_B24,
                 PredefinedWavelengths.GLI_B25,
                 PredefinedWavelengths.GLI_B26,
                 PredefinedWavelengths.GLI_B27,
                 PredefinedWavelengths.GLI_B28,
                 PredefinedWavelengths.GLI_B29,
                 PredefinedWavelengths.GLI_B30]
    elif target_sensor == 'landsat_etm':
        run_sixs_params = SixSHelpers.Wavelengths.run_landsat_etm
        bands = [PredefinedWavelengths.LANDSAT_ETM_B1,
                 PredefinedWavelengths.LANDSAT_ETM_B2,
                 PredefinedWavelengths.LANDSAT_ETM_B3,
                 PredefinedWavelengths.LANDSAT_ETM_B4,
                 PredefinedWavelengths.LANDSAT_ETM_B5,
                 PredefinedWavelengths.LANDSAT_ETM_B7]
    elif target_sensor == 'landsat_mss':
        run_sixs_params = SixSHelpers.Wavelengths.run_landsat_mss
        bands = [PredefinedWavelengths.LANDSAT_MSS_B1,
                 PredefinedWavelengths.LANDSAT_MSS_B2,
                 PredefinedWavelengths.LANDSAT_MSS_B3,
                 PredefinedWavelengths.LANDSAT_MSS_B4]
    elif target_sensor == 'landsat_oli':
        run_sixs_params = SixSHelpers.Wavelengths.run_landsat_oli
        bands = [PredefinedWavelengths.LANDSAT_OLI_B1,
                 PredefinedWavelengths.LANDSAT_OLI_B2,
                 PredefinedWavelengths.LANDSAT_OLI_B3,
                 PredefinedWavelengths.LANDSAT_OLI_B4,
                 PredefinedWavelengths.LANDSAT_OLI_B5,
                 PredefinedWavelengths.LANDSAT_OLI_B6,
                 PredefinedWavelengths.LANDSAT_OLI_B7,
                 PredefinedWavelengths.LANDSAT_OLI_B8,
                 PredefinedWavelengths.LANDSAT_OLI_B9]
    elif target_sensor == 'landsat_tm':
        run_sixs_params = SixSHelpers.Wavelengths.run_landsat_tm
        bands = [PredefinedWavelengths.LANDSAT_TM_B1,
                 PredefinedWavelengths.LANDSAT_TM_B2,
                 PredefinedWavelengths.LANDSAT_TM_B3,
                 PredefinedWavelengths.LANDSAT_TM_B4,
                 PredefinedWavelengths.LANDSAT_TM_B5,
                 PredefinedWavelengths.LANDSAT_TM_B7]
    elif target_sensor == 'meris':
        run_sixs_params = SixSHelpers.Wavelengths.run_meris
        bands = [PredefinedWavelengths.MERIS_B1,
                 PredefinedWavelengths.MERIS_B2,
                 PredefinedWavelengths.MERIS_B3,
                 PredefinedWavelengths.MERIS_B4,
                 PredefinedWavelengths.MERIS_B5,
                 PredefinedWavelengths.MERIS_B6,
                 PredefinedWavelengths.MERIS_B7,
                 PredefinedWavelengths.MERIS_B8,
                 PredefinedWavelengths.MERIS_B9,
                 PredefinedWavelengths.MERIS_B10,
                 PredefinedWavelengths.MERIS_B11,
                 PredefinedWavelengths.MERIS_B12,
                 PredefinedWavelengths.MERIS_B13,
                 PredefinedWavelengths.MERIS_B14,
                 PredefinedWavelengths.MERIS_B15]
    elif target_sensor == 'modis':
        run_sixs_params = SixSHelpers.Wavelengths.run_modis
        bands =[PredefinedWavelengths.MODIS_B1,
                PredefinedWavelengths.MODIS_B2,
                PredefinedWavelengths.MODIS_B3,
                PredefinedWavelengths.MODIS_B4,
                PredefinedWavelengths.MODIS_B5,
                PredefinedWavelengths.MODIS_B6,
                PredefinedWavelengths.MODIS_B7,
                PredefinedWavelengths.MODIS_B8]
    elif target_sensor == 'polder':
        run_sixs_params = SixSHelpers.Wavelengths.run_polder
        bands = [PredefinedWavelengths.POLDER_B1,
                 PredefinedWavelengths.POLDER_B2,
                 PredefinedWavelengths.POLDER_B3,
                 PredefinedWavelengths.POLDER_B4,
                 PredefinedWavelengths.POLDER_B5,
                 PredefinedWavelengths.POLDER_B6,
                 PredefinedWavelengths.POLDER_B7,
                 PredefinedWavelengths.POLDER_B8]
    elif target_sensor == 'spot_hrv':
        run_sixs_params = SixSHelpers.Wavelengths.run_spot_hrv
        bands = [PredefinedWavelengths.SPOT_HRV2_B1,
                 PredefinedWavelengths.SPOT_HRV2_B2,
                 PredefinedWavelengths.SPOT_HRV2_B3]
    elif target_sensor == 'spot_vgt':
        run_sixs_params = SixSHelpers.Wavelengths.run_spot_vgt
        bands = [PredefinedWavelengths.SPOT_VGT_B1,
                 PredefinedWavelengths.SPOT_VGT_B2,
                 PredefinedWavelengths.SPOT_VGT_B3,
                 PredefinedWavelengths.SPOT_VGT_B4]
    elif target_sensor == 'vnir':
        run_sixs_params = SixSHelpers.Wavelengths.run_vnir
    elif target_sensor == 'whole_range':
        run_sixs_params = SixSHelpers.Wavelengths.run_whole_range
    else:
        raise OSError('Unsupported sensor configuration')
        
        