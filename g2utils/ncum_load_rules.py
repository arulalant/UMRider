'''
This module supports for 'unknown' variables while loading ncum output 
(fieldsfiles).
Authour : Arulalan.T <arulalan@ncmrwf.gov.in>
'''
from cf_units import Unit
from iris.fileformats.grib import GribWrapper
from iris.coords import DimCoord, AuxCoord
import iris.analysis.cartography as icart
import cartopy.crs as ccrs
from numpy import meshgrid, array

regular2DLatitude, regular2DLongitude = None, None
warning = False
windCount = 0

ncumSTASH_Vs_cf = {
'm01s01i202': ('surface_net_downward_shortwave_flux', None, 'W m-2', None),
'm01s01i216': ('surface_diffuse_downwelling_shortwave_flux_in_air', None, 'W m-2', None),
'm01s03i229': ('water_evaporation_flux_from_soil', None, 'kg m-2', None),
'm01s30i403': ('atmosphere_mass_content_of_dust_dry_aerosol_particles', None, 'kg m-2', None),
# though m01s05i233 already standard_name exist, but wrongly assinged as 
# mass_fraction_of_convective_cloud_liquid_water_in_air in um_cf_map.py. So 
# in the below we are correcting the cf_standard_name of m01s05i233.
'm01s05i233': ('atmosphere_convective_available_potential_energy_wrt_surface', None, 'J kg-1', None),
'm01s05i234': ('atmosphere_convective_inhibition_wrt_surface', None, 'J kg-1', None),
# 50meter B-Grid U component wind 
'm01s15i212': ('x_wind', None, 'm s-1', 50),
# 50meter B-Grid V component wind   
'm01s15i213': ('y_wind', None, 'm s-1', 50),
}

duplicateSTASH_vs_cf = {
# duplicate STASH but cell_methods are different, variables comes here.
# STASH : {cell_method1 : (standard_name, long_name, unit, height),
#          cell_method2 : (standard_name, long_name, unit, height)}
'm01s03i236': {'maximum': (None, 'air_temperature_maximum', 'K', None),
               'minimum': (None, 'air_temperature_minimum', 'K', None)}
}

G2Param_vs_cf = {
# whose  surface level type as top of atmosphere (8)
# grib version, discipline, parameter category, parameter no, typeOfFirstFixedSurface
(2, 0, 5, 4, 8): ('toa_outgoing_longwave_flux', None, 'W m-2'),
(2, 0, 4, 8, 8): ('toa_outgoing_shortwave_flux', None, 'W m-2'), 
(2, 0, 4, 11, 8): ('toa_outgoing_shortwave_flux_assuming_clear_sky', None, 'W m-2'),
(2, 0, 4, 7, 8): ('toa_incoming_shortwave_flux', None, 'W m-2'),   
# whose  surface level type as tropopause (7)
(2, 0, 3, 0, 7): ('tropopause_air_pressure', None, 'Pa'),
(2, 0, 0, 0, 7): ('tropopause_air_temperature', None, 'K'),
(2, 0, 3, 6, 7): ('tropopause_altitude', None, 'm'),
# Need to load from NCMRWF Local Table entries Begin
(2, 0, 4, 198, 8): ('toa_outgoing_shortwave_flux_assuming_clear_sky', None, 'W m-2'),
(2, 0, 5, 195, 8): ('toa_outgoing_longwave_flux_assuming_clear_sky', None, 'W m-2'), 
(2, 0, 1, 192, 1): ('fog_area_fraction', None, '%'),
# grib version, discipline, parameter category, parameter no, typeOfFirstFixedSurface, scaleFactorOfFirstFixedSurface
(2, 3, 1, 192, 1, 7): (None, 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol_at_0.38um', '1'),
(2, 3, 1, 193, 1, 7): (None, 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol_at_0.44um', '1'), 
(2, 3, 1, 194, 1, 7): (None, 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol_at_0.55um', '1'), 
(2, 3, 1, 195, 1, 7): (None, 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol_at_0.67um', '1'), 
(2, 3, 1, 196, 1, 7): (None, 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol_at_0.87um', '1'), 
(2, 3, 1, 197, 1, 7): (None, 'atmosphere_optical_thickness_due_to_dust_ambient_aerosol_at_1.02um', '1'),
# Need to load from NCMRWF Local Table entries end 
}
    
def update_cf_standard_name(cube, field, filename):
    """
    This call back function used to update the cubes standard_name if not 
    exists (i.e unknown) or to correct the cf_standard_name of particular
    cube (if wrongly assigned in um_cf_map.py) while loading cubes via iris.
    
    >>> import iris
    >>> cubes = iris.load('umglaa_pf000', callback=update_cf_standard_name)
    >>> cubes = iris.load('um_ana.grib2', callback=update_cf_standard_name)
    
    Arulalan.T
    12-Feb-2016
    """ 
    
    if cube.attributes:
        # loading from pp/fieldsfiles
        varSTASH = str(cube.attributes['STASH'])
        if varSTASH in ncumSTASH_Vs_cf:
            # get correct standard_name and units
            sname, lname, unit, height = ncumSTASH_Vs_cf[varSTASH]
            # update cube's standard_name and its units
            cube.standard_name = sname
            cube.long_name = lname
            cube.units = Unit(unit)
            if height:
                heightAx = DimCoord(array([float(height)]), standard_name='height',
                             units=Unit('m'), attributes={'positive': 'up'})
                cube.add_aux_coord(heightAx) 
        elif varSTASH in duplicateSTASH_vs_cf:
            ccm = cube.cell_methods
            if ccm:
                ccmm = ccm[0].method 
                if ccmm in duplicateSTASH_vs_cf[varSTASH]:
                    # get correct standard_name, units and height 
                    sname, lname, unit, height = duplicateSTASH_vs_cf[varSTASH][ccmm]
                    # update cube's standard_name and its units
                    cube.standard_name = sname
                    cube.long_name = lname
                    cube.units = Unit(unit)
                    if height:
                        heightAx = DimCoord(array([float(height)]), standard_name='height',
                                     units=Unit('m'), attributes={'positive': 'up'})
                        cube.add_aux_coord(heightAx)                        
            # end of if ccm:
        # end of if varSTASH in ncumSTASH_Vs_cf:
    elif isinstance(field, GribWrapper):
        # loading from grib file
        if field.editionNumber == 2 and (field.parameterNumber > 191 or 
                             field.typeOfFirstFixedSurface in [7, 8] or 
                             field.typeOfSecondFixedSurface in [7]):        
            # loaded from grib2 file 
            # loaded tropopause or toa or local table entries 
            discipline = field.discipline
            pcategory = field.parameterCategory
            pnumber = field.parameterNumber
            fsfc = field.typeOfFirstFixedSurface           
            ssfc = field.typeOfSecondFixedSurface
            if ssfc in [7]:
                G2ParamKey = (2, discipline, pcategory, pnumber, fsfc, ssfc)
            else:
                G2ParamKey = (2, discipline, pcategory, pnumber, fsfc)
            if G2ParamKey in G2Param_vs_cf:
                # get correct standard_name and units
                sname, lname, unit = G2Param_vs_cf[G2ParamKey]
                # update cube's standard_name and its units
                cube.standard_name = sname
                cube.long_name = lname
                cube.units = Unit(unit)
            # end of if G2ParamKey in G2Param_vs_cf:
# end of def update_cf_standard_name_if_not_exists(cube, field, filename):

def setCubeRegularLatLon(cube, field, filename):
        
    # https://github.com/SciTools/iris/issues/448
    glon = cube.coord('grid_longitude')
    rotated_lons = glon.points
    glat = cube.coord('grid_latitude')
    rotated_lats = glat.points
    pole_lat = glat.coord_system.grid_north_pole_latitude
    pole_lon = glat.coord_system.grid_north_pole_longitude

    rotated_lats, rotated_lons = meshgrid(rotated_lons, rotated_lats)

    lons, lats = icart.unrotate_pole(rotated_lons, rotated_lats, pole_lon, pole_lat)

    ## the below lat is 2D array
    regular2DLatitude = AuxCoord(lons, standard_name='longitude', units='degree_north')
    ## the below lon is 2D array
    regular2DLongitude = AuxCoord(lats, standard_name='latitude', units='degree_east')
    
    cube.remove_coord('grid_latitude')
    cube.remove_coord('grid_longitude')
    
    cube.add_aux_coord(regular2DLatitude, [0, 1])
    cube.add_aux_coord(regular2DLongitude, [0, 1])
    ## https://groups.google.com/forum/#!searchin/scitools-iris/unrotate/scitools-iris/m5B2752zRjQ/Xmcp5yaqAQAJ
    ## Andrew Dawson's reply : If you really want the data on a Plate Carree projection with 1D dimension coordinates (lat/lon) then you will need to regrid/reproject to move the actual data points onto a regular grid in lat/lon space.
    
    ## So after this, if we regrid to our regular lat, lon then regridded cube lat, lon will become two 1D array.
    ## Hopefully !!!!
        
# end of def setCubeRegularLatLon(cube, field, filename):
