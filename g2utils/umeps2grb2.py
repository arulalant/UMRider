#!/usr/bin/env python

__author__ = 'arulalant'
__version__ = 'v1.0.1'
__long_name__ = 'NCUM Ensembles Parallel Rider'

"""
Inputs: NCUM fieldsfile / pp format files

Outputs: WMO-NCEP Grib2 format files

Disclaimers (if any!)
This is just test code as of now and is meant for a specific purpose only!

Copyright: ESSO-NCMRWF, MoES, 2015-2016, 2016-2017.

Author : Arulalan.T
initial code : 15-Mar-2016
latest Update : 30-Mar-2016
"""

# -- Start importing necessary modules
import os, sys, time, subprocess, errno
import numpy 
import iris
import gribapi
from cf_units import Unit
import multiprocessing as mp
import multiprocessing.pool as mppool       
# We must import this multiprocessing.pool explicitly, it is not imported
# by the top-level multiprocessing module.
import datetime
from iris.time import PartialDateTime
from cubeutils import cubeAverager, cubeSubtractor
from ncum_load_rules import update_cf_standard_name
from um2grb2 import (createDirWhileParallelRacing, getCubeData, myLog, 
             __getAnlFcstFileNameIdecies__, __genAnlFcstOutFileName__, 
            getCubeAttr, _NoDaemonProcess, _MyPool)

# End of importing business

# We have to make sure that strict_grib_load as False, since we have to 
# read the cubes from grib2 to re-order the variables. True throws an error
# while reading for tweaked_messages (say pf varibles)
iris.FUTURE.strict_grib_load = False
iris.FUTURE.netcdf_promote = True
iris.FUTURE.netcdf_no_unlimited = True
iris.FUTURE.cell_datetime_objects = True
# -- Start coding
# create global _lock_ object
_lock_ = mp.Lock()

# global path variables
g2ctl = "/gpfs2/home/umtid/Softwares/grib2ctl/g2ctl.pl"
grib2ctl = "/gpfs2/home/umtid/Softwares/grib2ctl/grib2ctl.pl"
gribmap = "/gpfs1/home/Libs/GNU/GRADS/grads-2.0.2.oga.1/Contents/gribmap"
cnvgrib = "/gpfs1/home/Libs/INTEL/cnvgrib-1.4.0/cnvgrib"

# other global variables
__LPRINT__ = False
__utc__ = '00'
__outFileType__ = 'ana'
# start and step hour in short forecast files
__anl_step_hour__ = 6
# start and step hour in long forecast files
__start_step_long_fcst_hour__ = 6
# maximum long forecast hours produced by model
__max_long_fcst_hours__ = 240
# grib1 file suffix
__grib1FilesNameSuffix__ = '.grib1'
# flag for removing grib2 files after grib1 has been converted 
__removeGrib2FilesAfterGrib1FilesCreated__ = False

# Defining default out grib2 file name structure for analysis 
__anlFileNameStructure__ = ('um_ana', '_', '*HHH*', 'hr', '_', 
                            '*YYYYMMDD*', '_', '*ZZ*', 'Z', '.grib2')

# Defining default out grib2 file name structure for forecast                             
__fcstFileNameStructure__ = ('um_prg', '_', '*HHH*', 'hr', '_', 
                            '*YYYYMMDD*', '_', '*ZZ*', 'Z', '.grib2')
                            
# the _convertVars_ is global list which should has final variables list of 
# tuples (varName, varSTASH) will be converted, otherwise default variables 
# of this module will be converted!
_convertVars_ = []
_removeVars_ = []   # used to store temporary vars 
_current_date_ = None
_startT_ = None
_tmpDir_ = None
_inDataPath_ = None
_opPath_ = None
_doRegrid_ = False
_targetGrid_ = None
_targetGridRes_ = None
_requiredLat_ = None
_requiredLon_ = None
_requiredPressureLevels_ = None
_preExtension_ = '_unOrdered'
_createGrib2CtlIdxFiles_ = True
_createGrib1CtlIdxFiles_ = False
_convertGrib2FilestoGrib1Files_ = False
__setGrib2TableParameters__ = None
_ensemble_count_ = 44
# store out grib2 files name for the purpose of creat ctl files in parallel
__outg2files__ = []
# global ordered variables (the order we want to write into grib2)
_orderedVars_ = {'PressureLevel': [
## Pressure Level Variable names & STASH codes
('geopotential_height', 'm01s16i202'),        
('x_wind', 'm01s15i243'), 
('y_wind', 'm01s15i244'),   
('upward_air_velocity', 'm01s15i242'),
('air_temperature', 'm01s16i203'),
('relative_humidity', 'm01s16i256'),
('specific_humidity', 'm01s30i205')],

## Non Pressure Level Variable names & STASH codes
'nonPressureLevel': [
('tropopause_altitude', 'm01s30i453'),
('tropopause_air_temperature', 'm01s30i452'),
('tropopause_air_pressure', 'm01s30i451'),
('surface_air_pressure', 'm01s00i409'),
('air_pressure_at_sea_level', 'm01s16i222'),
('surface_temperature', 'm01s00i024'),
('relative_humidity', 'm01s03i245'), 
('specific_humidity', 'm01s03i237'),
('air_temperature', 'm01s03i236'),
('dew_point_temperature', 'm01s03i250'),
('atmosphere_convective_available_potential_energy_wrt_surface', 'm01s05i233'), # CAPE
('atmosphere_convective_inhibition_wrt_surface', 'm01s05i234'), #CIN
('high_type_cloud_area_fraction', 'm01s09i205'),
('medium_type_cloud_area_fraction', 'm01s09i204'),
('low_type_cloud_area_fraction', 'm01s09i203'), 
('x_wind', 'm01s03i209'), 
('y_wind', 'm01s03i210'),    
('visibility_in_air', 'm01s03i247'),
('precipitation_amount', 'm01s05i226'),
('stratiform_snowfall_amount', 'm01s04i202'),
('convective_snowfall_amount', 'm01s05i202'),
('stratiform_rainfall_amount', 'm01s04i201'),
('convective_rainfall_amount', 'm01s05i201'),
('rainfall_flux', 'm01s05i214'),
('snowfall_flux', 'm01s05i215'),
('precipitation_flux', 'm01s05i216'),
('fog_area_fraction', 'm01s03i248'),
('toa_incoming_shortwave_flux', 'm01s01i207'), 
('toa_outgoing_shortwave_flux', 'm01s01i205'),
('toa_outgoing_shortwave_flux_assuming_clear_sky', 'm01s01i209'),  
('toa_outgoing_longwave_flux', 'm01s02i205'),
('toa_outgoing_longwave_flux_assuming_clear_sky', 'm01s02i206'),   
('surface_upward_latent_heat_flux', 'm01s03i234'),
('surface_upward_sensible_heat_flux', 'm01s03i217'),
('surface_downwelling_shortwave_flux_in_air', 'm01s01i235'),
('surface_downwelling_longwave_flux', 'm01s02i207'),
('surface_net_downward_longwave_flux', 'm01s02i201'), 
('surface_net_downward_shortwave_flux', 'm01s01i202'),
('atmosphere_boundary_layer_thickness', 'm01s00i025'),
('atmosphere_optical_thickness_due_to_dust_ambient_aerosol', 'm01s02i422'),
('moisture_content_of_soil_layer', 'm01s08i223'),  # 4 layers 
# single layer, this must be after 4 layers as in order 
('soil_moisture_content', 'm01s08i208'),       # single layer  
## though moisture_content_of_soil_layer and volumetric_moisture_of_soil_layer
## has same STASH code, but we must include seperate entry here.
('volumetric_moisture_of_soil_layer', 'm01s08i223'), # 4 layers
# single layer, this must be after 4 layers as in order 
('volumetric_moisture_of_soil_layer', 'm01s08i208'), # single layer
('soil_temperature', 'm01s03i238'),  
('land_binary_mask', 'm01s00i030'),
('sea_ice_area_fraction', 'm01s00i031'),
('sea_ice_thickness', 'm01s00i032'),
# the snowfall_amount might be changed as 
# liquid_water_content_of_surface_snow by convert it into
# water equivalent of snow amount, before re-ordering itself.
('liquid_water_content_of_surface_snow', 'm01s00i023'),
# the below one is for orography which presents only in analysis 00 file.
# so we must keep this as the last one in the ordered variables!
('surface_altitude', 'm01s00i033')],
}

# Define _accumulationVars_
# The following variables should be 6-hourly accumulated, but model
# produced as 1-hourly accumulation. So we need to sum of 6-hours data to 
# get 6-hourly accumulation.
# rainfall_flux, snowfall_flux, precipitation_flux are not accumulated 
# vars, since those are averaged rain rate (kg m-2 s-1). 
# But the following vars unit is (kg m-2), accumulated vars.  
_accumulationVars_ = []                      

## Define _ncfilesVars_
## the following variables need to be written into nc file, initially for 
## some reason (like either it may contain soil_model_level_number or 
## duplicate grib param where typeOfFirstFixedSurface [i.e toa, tropopause]
## is not implemented yet), but at the end of the program (after re-ordering)
## these intermediate nc files will be deleted automatically! 
## while storing into nc file, there wont be much problem and reading from
## nc file also wont be problem to load cf_standard_name into cube and 
## followed by storing into grib2 file. All because we need to write variables
## in the way we need (i.e. ordered always)!!
_ncfilesVars_ = []
                 
## Define _ncmrGrib2LocalTableVars_
## the following variables need to be set localTableVersion no as 1 and
## master table version no as 255 (undefined), since WRF grib2 table doesnt
## support for the following variables. So we created our own local table.
_ncmrGrib2LocalTableVars_ = ['fog_area_fraction',
                            'toa_outgoing_longwave_flux_assuming_clear_sky',   
                            'toa_outgoing_shortwave_flux_assuming_clear_sky',
                  'atmosphere_optical_thickness_due_to_dust_ambient_aerosol']

## Define _maskOverOceanVars_
## the following variables need to be set mask over ocean because the original
## model itself producing mask over ocean. but when we are doing regrid it 
## couldnt retain the mask ! dont know why ! So using land_binary_mask 
## variable, we are resetting mask over ocean for the following vars.
_maskOverOceanVars_ = []

## Define dust aerosol optical thickness of model pseudo level with its 
## corresponding micron / micro wavelength. We need to tweak with following 
## information before writing into final grib2 file.
_aod_pseudo_level_var_ = {}

## Define _depedendantVars_ where A is key, B is value. A is depedendant on B,
## B is not. B not necessarily to be written in out file. User may just specify
## only A in var.cfg configure file.
_depedendantVars_ = {}
                                            

# start definition #2
def getVarInOutFilesDetails(inDataPath, fname, hr):
    """
    This definition module gets the required variables from the passed
    cube as per the WRF-Variables.txt file.

    :param inDataPath: data path which contains data and hour.
    :param fname: filename of the fieldsfile that has been passed as a string.

    :return: varNamesSTASH: a list of tuples (Variable name and its STASH code) 
    :return: fcstHours: Time slices of the cube as an array/scalar - integer (number)
    :return: doMultiHourlyMean: Logical expression as either True or False, indicating
                            whether the field is instantaneous or accumulated
    :return: infile: It returns absolute path of infile by inDataPath and fname.
                     Also it updates inDataPath yesterday, hour for analysis pf files
    :return: outfile: It returns outfile absolute path with ana or fcst type 
                      along with date and hour.

    """
    global __anl_step_hour__, __start_step_long_fcst_hour__
    
    hr = int(hr)
    
    infile = os.path.join(inDataPath, fname)    
            
    ##### ENSEMBLES FILE BEGIN     
    if 'pb' in fname:                   
        varNamesSTASH = [('air_pressure_at_sea_level', 'm01s16i222'),
            ('air_temperature', 'm01s03i236'),
            ('relative_humidity', 'm01s03i245'),
            ('specific_humidity', 'm01s03i237'),
            ('precipitation_amount', 'm01s05i226'),
            ('x_wind', 'm01s03i209'), 
            ('y_wind', 'm01s03i210'),]
        # the cube contains Instantaneous data at every 6-hours.        
                
        if __start_step_long_fcst_hour__ == 6:
            # applicable only for 6 hour instantaneous/intervals
            fcstHours = numpy.array([6, 12, 18, 24]) + hr        
        # we are extracting at particular instantaneous value, so no need to 
        # do hourly mean.
        doMultiHourlyMean = False
        
    elif 'xbfti' in fname and fname.endswith('.pp0'):            
        # consider variable
        varNamesSTASH = [('air_pressure_at_sea_level', 'm01s16i222'),
                    ('air_temperature', 'm01s16i203'),  
                    ('geopotential_height', 'm01s16i202'),                    
                    ('relative_humidity', 'm01s16i256'), 
                    ('surface_air_pressure', 'm01s00i409'),                   
                    ('x_wind', 'm01s15i243'),
                    ('y_wind', 'm01s15i244'),
                    ('upward_air_velocity', 'm01s15i242')]
        # the cube contains Instantaneous data at every 24-hours.        
        if __start_step_long_fcst_hour__ == 24:
            # applicable only for 24 hour instantaneous/intervals
            fcstHours = numpy.arange(0, 241, 24)       
        # we are extracting at particular instantaneous value, so no need to 
        # do hourly mean.      
        doMultiHourlyMean = False    
        
    elif 'xbfti' in fname and fname.endswith('.pp2'):            
        # consider variable
        varNamesSTASH = [('air_temperature_maximum', 'm01s03i236'),  
                    ('air_temperature_minimum', 'm01s03i236'),  
                    ('precipitation_amount', 'm01s05i226')]              
        # the cube contains accumulated/min/max data at every 24-hours.        
        if __start_step_long_fcst_hour__ == 24:
            # applicable only for 24 hour intervals
            fcstHours = numpy.arange(12, 240, 24)       
        # we are extracting at particular instantaneous value, so no need to 
        # do hourly mean.      
        doMultiHourlyMean = False    
        
    ##### FORECAST FILE END
    else:
        raise ValueError("Filename '%s' not implemented yet!" % fname)
    # end if-loop

    return varNamesSTASH, fcstHours, doMultiHourlyMean, infile
# end of def getVarInOutFilesDetails(inDataPath, fname, hr):

def packEnsembles(arg):
    
    global _targetGrid_, _targetGridRes_,  _startT_, _inDataPath_, _opPath_, \
            _preExtension_, _ncfilesVars_, _requiredLat_, _requiredLon_, \
            _doRegrid_, __utc__, _requiredPressureLevels_, __LPRINT__, \
            __outg2files__, _lock_
           
    infiles, varNamesSTASHFcstHour = arg
    varName, varSTASH, fhr = varNamesSTASHFcstHour
    
    if (varName, varSTASH) in [('precipitation_amount', 'm01s05i226')]:
        # update the forecast hour, since precipitation_amount is accumulated
        # var, not instantaneous one.
        fhr -= 3
    # end of if (varName, varSTASH) in [('precipitation_amount', 'm01s05i226')]:
    
    simulated_hr = __utc__
    # define variable name constraint
    varConstraint = iris.Constraint(name=varName)
    # define varibale stash code constraint
    STASHConstraint = iris.AttributeConstraint(STASH=varSTASH)    
    # 
    forecast_period_constraint = iris.Constraint(forecast_period=fhr)            
     # define (simulated_hr) forecast_reference_time constraint
    fcstRefTimeConstraint = iris.Constraint(forecast_reference_time=PartialDateTime(hour=simulated_hr))
    # Define default lat, lon, pressure contraint (None just bring model global data)
    latConstraint, lonConstraint, pressureConstraint = None, None, None
    if _requiredLat_: 
        # make constraint of required latitude
        latConstraint = iris.Constraint(latitude=lambda cell: 
                                _requiredLat_[0] <= cell <= _requiredLat_[-1])
    if _requiredLon_:
        # make constraint of required longitude
        lonConstraint = iris.Constraint(longitude=lambda cell: 
                                _requiredLon_[0] <= cell <= _requiredLon_[-1])
    
    if _requiredPressureLevels_:
        # make constraint of required pressure 
        # To slice out particular pressure levels (like 850, 200, 1000 alone)
        # then the following way is essential.       
        pressureConstraint = iris.Constraint(pressure=lambda cell: 
                                int(cell.point) in _requiredPressureLevels_)
    
    loadConstraints = varConstraint & STASHConstraint & forecast_period_constraint & latConstraint & lonConstraint   
    ensembleData, controlVar, dshape = None, None, None
    print "packEnsembles Started using", infiles
    for idx, infile in enumerate(infiles):
        print "extracting ensemble data", infile
        # load ensemble cube with all constraints
        ensCube = getCubeData(infile, constraints=loadConstraints)[0]        
        # extract pressure levels
        if pressureConstraint and ensCube.coords('pressure'): 
            ensCube = ensCube.extract(pressureConstraint)
        # ene of if pressureConstraint and tmpCube.coords('pressure'): 

        if ensCube.has_lazy_data():
            print "Loaded", ensCube.standard_name, "into memory",
            ## By accessing tmpCube.data (even for printing), the full 
            ## data has been loaded into memory instead of being lazy 
            ## data. Especially for dust aod, we must make it as fully 
            ## loaded otherwise full data will be treated as zeros only 
            ## instead of 6 pseudo_level data.
            print "- min", ensCube.data.min(), "max", ensCube.data.max(),
            print "has_lazy_data =", ensCube.has_lazy_data()
        # end of if ensCube.has_lazy_data():

        # interpolate it as per targetGridResolution deg resolution by 
        # setting up sample points based on coord        
        if _doRegrid_:
            print "\n Regridding data to %sx%s degree spatial resolution \n" % (_targetGridRes_, _targetGridRes_)
            if __LPRINT__: print "From shape", ensCube.shape 
            try:
                # This lienar interpolate will do extra polate over ocean even 
                # though original data doesnt have values over ocean and wise versa.
                # So lets be aware of this.
                # DO NOT APPLY iris.analysis.Linear(extrapolation_mode='mask'), 
                # which writes nan every where for the snowfall_flux,  
                # rainfall_flux, precipitation_flux. So donot apply that.               
                regdCube = ensCube.interpolate(_targetGrid_, iris.analysis.Linear())
            except Exception as e:
                print "ALERT !!! Error while regridding!! %s" % str(e)
                print " So skipping this without saving data"
                continue
            # end of try:      
        else:
            # do not apply regrid. this is temporary fix. 
            regdCube = ensCube
        # end of if _doRegrid_:
        
        if not idx:
            # assumming firs file points to control run 
            controlVar = regdCube
            # introduce ensemble dimension at first axis 
            dshape = list(controlVar.data.shape)
            dshape.insert(0, 1)    
            ensembleData = controlVar.data.reshape(dshape)
            # make memory free 
            del regdCube
        else:
            # introduce ensemble dimension at first axis 
            edata = regdCube.data.reshape(dshape)           
            # pack all ensembles next to each other in first dimension
            ensembleData = numpy.concatenate((ensembleData, edata))
            # make memory free 
            del edata, regdCube, ensCube
        # end of if not idx:
    # end of for idx, infile in enumerate(infiles):
    print "taken into memory of all ensembles", ensembleData.shape 
    # convert data into masked array
    ensembleData = numpy.ma.masked_array(ensembleData, dtype=numpy.float64)
    # http://www.cpc.ncep.noaa.gov/products/wesley/g2grb.html
    # Says that 9.999e+20 value indicates as missingValue in grib2
    # by default g2ctl.pl generate "undefr 9.999e+20", so we must 
    # keep the fill_value / missingValue as 9.999e+20 only.
    numpy.ma.set_fill_value(ensembleData, 9.999e+20)
            
    totEns = len(ensembleData)
    # create ensemble coordinate
    enscoord = iris.coords.DimCoord(numpy.arange(totEns, dtype=numpy.int32), 
                         standard_name='realization', units=Unit('no_unit'), 
                                                long_name='ensemble_member')
    # get list of dimension coordinates
    dim_coords = list(controlVar.dim_coords)
    # insert ensemble dimension at first axis 
    dim_coords.insert(0, enscoord)
    # generate list of tuples contain index and coordinate
    dim_coords = [(coord, i) for i,coord in enumerate(dim_coords)]
    # get all other dimensions
    aux_coords = list(controlVar.aux_coords)
    aux_factories = controlVar.aux_factories
    t = controlVar.coords('time')[0]
    fp = controlVar.coords('forecast_period')[0]
    ft = controlVar.coords('forecast_reference_time')[0]
    # create ensemble packed cubes 
    ensembleData = iris.cube.Cube(ensembleData, controlVar.standard_name, 
                             controlVar.long_name, controlVar.var_name,
                               controlVar.units, controlVar.attributes, 
                                   controlVar.cell_methods, dim_coords)
    # add all time coordinates
    print "setting aux_coords to", ensembleData.shape, varName, fhr 
    ensembleData.add_aux_coord(fp)
    ensembleData.add_aux_coord(ft)
    ensembleData.add_aux_coord(t)
    # create cell method for ensembles
    cm = iris.coords.CellMethod('realization', ('realization',), 
                           intervals=('1',), comments=(' ENS',))
    # add cell_methods to the ensembleData                        
    if controlVar.cell_methods:             
        ensembleData.cell_methods = (cm, controlVar.cell_methods[0])
    else:
        ensembleData.cell_methods = (cm,)
    print ensembleData
    # make memory free 
    del controlVar  
    
    print "To ensembleData shape", ensembleData.shape  

    # get the regridded ensembles meta data 
    varName, varSTASH, fcstTm, refTm, lat1, lon1 = getCubeAttr(ensembleData)

    if fcstTm.bounds is not None:                
        # this is needed for forecast 00th simulated_hr
        # get the last hour from bounds
        hr = str(int(fcstTm.bounds[-1][-1]))
        if __LPRINT__: print "Bounds comes in ", hr, fcstTm.bounds, fileName                        
    else:
        # get the fcst time point 
        # this is needed for analysis/forecast 00th simulated_hr
        hr = str(int(fcstTm.points))
        if __LPRINT__: print "points comes in ", hr, fileName 
    # end of if fcstTm.bounds:
    
    outFileNameStructure = __fcstFileNameStructure__
    # get the out fileName Structure based on pre / user defined indecies                       
    outFnIndecies = __getAnlFcstFileNameIdecies__(outFileNameStructure)
    # get the file name extension
    fileExtension = outFileNameStructure[-1]  
                              
    # generate the out file name based on actual informations                                 
    outFn = __genAnlFcstOutFileName__(outFileNameStructure, 
                         outFnIndecies, _current_date_, hr, 
                                   __utc__, _preExtension_) 
    # get the file full name except last extension, for the purpose
    # of writing intermediate nc files
    ofname = outFn.split(fileExtension)[0]                    
    ncfile = False
    print "outFn = ", outFn
    if (varName, varSTASH) in _ncfilesVars_:
        # other than soil_model_level_number, few variables may be 
        # need to write into nc file and then convert to grib2. why 
        # because of duplicate grib param id (but actually not, if 
        # we implement typeOfFirstFixedSurface). so we are stoing into 
        # nc file, then load into memory (cf_standard_name) while 
        # re-ordering, followed by save into grib2 file. cf -> grib2 
        # dictionary may not throw error, due to different key cfname.
        ncfile = True
        outFn = varSTASH + '_'+ ofname + '.nc'
    # end of if (varName, varSTASH) in _ncfilesVars_:
    
    # append out grib2 files for the purpose of creating ctl files.
    if not outFn in __outg2files__: __outg2files__.append(outFn)
    print "__outg2files__ = ", __outg2files__
    outFn = os.path.join(_opPath_, outFn)
    print "Going to be save into ", outFn
    print ensembleData
                
    try:                
        # _lock_ other threads / processors from being access same file 
        # to write other variables. 
        _lock_.acquire()
        if ncfile:
            iris.fileformats.netcdf.save(ensembleData, outFn)  # save nc file 
        else:
            # before save it, tweak the cubes by setting centre no and 
            # address other temporary issues before saving into grib2.
            iris.fileformats.grib.save_messages(tweaked_messages([ensembleData,]), 
                                                           outFn, append=True) # save grib2 file 
        # release the _lock_, let other threads/processors access this file.
        _lock_.release()
    except Exception as e:
        print "ALERT !!! Error while saving!! %s" % str(e)
        print " So skipping this without saving data"        
    # end of try:
    print "saved"
    # make memory free 
    del ensembleData
# end of def packEnsembles(arg):                     

def packEnsemblesInParallel(arg):

    global  _startT_, _inDataPath_, __start_step_long_fcst_hour__, __LPRINT__, \
            _opPath_, _ensemble_count_
   
    fpname, fpext, hr = arg 
    
    start_step_fcst_hour = __start_step_long_fcst_hour__
            
    if start_step_fcst_hour == 6:
        # 044_pb120
        ensembleFiles = [os.path.join(_inDataPath_, str(ens).zfill(3)+'_'+fpname+hr.zfill(3)) 
                                                for ens in range(0, _ensemble_count_+1, 1)]
        fileName = '000_' + fpname + '000'
    elif start_step_fcst_hour == 24:
        # xbfti_044.pp2
        ensembleFiles = [os.path.join(_inDataPath_, fpname+'_'+str(ens).zfill(3)+fpext) 
                                                for ens in range(0, _ensemble_count_+1, 1)]
        fileName = fpname + '_000' + fpext
    
    fname = os.path.join(_inDataPath_, fileName)
    
    # get variable indices
    varNamesSTASH, fcstHours, doMultiHourlyMean, infile = getVarInOutFilesDetails(_inDataPath_,
                                                                                 fileName, hr)
    
    for fname in ensembleFiles:
        if not os.path.isfile(fname): 
            print "The file doesn't exists: %s.. \n" %fname
            return  
        # end of if not os.path.isfile(fname): 
    # end of for fname in ensembleFiles:
    
    if _convertVars_:
        # load only needed variables from this file as well sort as per user specified ordered vars!
        varNamesSTASH = [vns for vns in _convertVars_ if vns in varNamesSTASH]
    
    if not varNamesSTASH:
        print "No varibale selected to load from the file '%s' " % fname
        if __LPRINT__: 
            print "Because global variable _convertVars_ doesn't contain any one of the following variables"
            print "\n".join([str(i+1)+' : ' + str(tu) for i, tu in enumerate(varNamesSTASH)])
        return None
    else:
        print "The following variables are going to be converted from file ", fname
        print "\n".join([str(i+1)+' : ' + str(tu) for i, tu in enumerate(varNamesSTASH)])
    # end of if not varNamesSTASH:
    
    print "Started Processing the file:  \n"     
    ensembleFiles_allConstraints_list = []
    for varName, varSTASH in varNamesSTASH:        
        for fhr in fcstHours:            
            allConstraints = [varName, varSTASH, fhr]     
            ensembleFiles_allConstraints_list.append((ensembleFiles, allConstraints))
    # end of for varName, varSTASH in varNamesSTASH:      
            
    ## get the no of childs process to create fcst ensemble files  
    nchild = len(ensembleFiles_allConstraints_list)     
    # create the no of child parallel processes
    inner_pool = mp.Pool(processes=nchild)
    print "Creating %i (daemon) workers and jobs in child." % nchild

    print "parallel ensemble begins for", varName, varSTASH
    # pass the (ensemblefileslist, allConstraints, pressureConstraint) as 
    # argument to take one fcst ensemble file per process / core to regrid it.
    results = inner_pool.map(packEnsembles, ensembleFiles_allConstraints_list)
    # closing and joining child pools      
    inner_pool.close() 
    inner_pool.join()
    # parallel end
    # end of for varName, varSTASH in varNamesSTASH:        
                                   
    print "Time taken to convert the file: %8.5f seconds \n" %(time.time()-_startT_)
    print "Finished converting file: %s into grib2 format for fcst file: %s \n" %(fpname, hr)
# end of def packEnsemblesInParallel(arg):

def createGrib2CtlIdxFiles(g2filepath, ftype='fcst'):
    ## g2ctl.pl usage option refer the below link 
    ## https://tuxcoder.wordpress.com/2011/08/31/how-to-install-g2ctl-pl-and-wgrib2-in-linux/
    ## though options says -verf for forecast end time, -0 for analysis time 
    ## -b for forecast start time, its all about setting reference in ctl file.
    ## Nothing more than that. We already set correct reference and forecast time bounds 
    ## in analysis files (whichever variables are actually taken from previous short forecast 0-6 hours).
    ## so here we no need to pass any options like -0 or -b. 
    ## By default g2ctl takes -verf option, same option we are passing 
    ## here to make sure that in future it will not affect.
    ctlfile = open(g2filepath+'.ctl', 'w')
    if ftype == 'fcst':
        # create ctl & idx files for forecast file
        # by default -verf as passed which takes end time of fcst bounds to set as base time.
        subprocess.call([g2ctl, '-ts6hr', '-verf', g2filepath], stdout=ctlfile)
        subprocess.call([gribmap, '-i', g2filepath+'.ctl'])  
        print "Successfully created control and index file using g2ctl !", g2filepath+'.ctl'       
# end of def createCtlIdxFiles(g2filepath):

def tweaked_messages(cubeList):
    global _ncmrGrib2LocalTableVars_, __setGrib2TableParameters__
    
    for cube in cubeList:
        for cube, grib_message in iris.fileformats.grib.as_pairs(cube):
            print "Tweaking begin ", cube.standard_name
            print cube 
            # post process the GRIB2 message, prior to saving
            gribapi.grib_set_long(grib_message, "centre", 29) # RMC of India
            gribapi.grib_set_long(grib_message, "subCentre", 0) # No subcentre
            print "reset the centre as 29"
            
            # ensembles tweak begin
            # http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-3.shtml 
            # 4 points ensemble forecast
            gribapi.grib_set_long(grib_message, "typeOfGeneratingProcess", 4)
            
            if cube.standard_name == 'precipitation_amount':
                print "precipitation_amount"
                print cube.coord("forecast_period"), cube.coord("forecast_period").bounds
            if cube.coord("forecast_period").bounds is None:       
                #http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-0.shtml 
                # template 01 would be better
                gribapi.grib_set(grib_message, "productDefinitionTemplateNumber", 1)
            else:
                # template 11 would be better
                gribapi.grib_set(grib_message, "productDefinitionTemplateNumber", 11)                   
                # if we set bounds[0][0] = 0, wgrib2 gives error for 0 fcst time.
                # so we need to set proper time intervals 
                # (typeOfTimeIncrement) as 2 as per below table.
                # http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-11.shtml
                # fileformats/grib/_save_rules.py-> set_forecast_time() ->
                # _non_missing_forecast_period() returns 'fp' as bounds[0][0]. 
                # but mean while lets fix by setting typeOfTimeIncrement=2.
                # http://www.cosmo-model.org/content/model/documentation/grib/pdtemplate_4.11.htm 
                gribapi.grib_set(grib_message, "typeOfTimeIncrement", 2)           
                print 'reset typeOfTimeIncrement as 2 for', cube.standard_name
            # end of if cube.coord("forecast_period").bounds is None:       
            
            # setting ensemble no   
            gribapi.grib_set(grib_message, "perturbationNumber",
                         int(cube.coord('realization').points[0]))
            # no encoding at present in Iris, set to missing
            gribapi.grib_set(grib_message, "numberOfForecastsInEnsemble", 255)
            #http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-6.shtml 
            # 3 would be better, since we keep on increasing ensmble points 
            # from 0 to 44
            gribapi.grib_set(grib_message, "typeOfEnsembleForecast", 3)
            # ensembles tweak end
            
            if cube.standard_name or cube.long_name:
                if cube.standard_name:
                    loc_longname = None
                    if cube.standard_name.startswith('air_pressure_at_sea_level'):
                        # we have to explicitly re-set the type of first fixed
                        # surfcae as Mean sea level (101)
                        gribapi.grib_set(grib_message, "typeOfFirstFixedSurface", 101)                     
                    # end of if cube.standard_name.startswith('tropopause'): 
                # end of if cube.standard_name:

                if cube.long_name:                     
                    loc_longname = [1 for lname in _ncmrGrib2LocalTableVars_ if cube.long_name.startswith(lname)]
                # end of if cube.long_name: 
                
                if cube.standard_name in _ncmrGrib2LocalTableVars_ or loc_longname:
                    # We have to enable local table version and disable the 
                    # master table only the special variables.
                    # http://www.cosmo-model.org/content/model/documentation/grib/grib2keys_1.htm 
                    # Above link says that tablesVersion must be set to 255, 
                    # then only local table will be enabled.
                    gribapi.grib_set_long(grib_message, "tablesVersion", 255)
                    # http://apt-browse.org/browse/debian/wheezy/main/i386/libgrib-api-1.9.16/1.9.16-2%2Bb1/file/usr/share/grib_api/definitions/grib2/section.1.def (line no 42)
                    # Above link says versionNumberOfGribLocalTables is alias 
                    # of LocalTablesVersion.        
                    # Set local table version number as 1 as per 
                    # ncmr_grib2_local_table standard.
                    gribapi.grib_set_long(grib_message, "versionNumberOfGribLocalTables", 1)
                # end of if cube.standard_name in _ncmrGrib2LocalTableVars_:
            # end of if cube.standard_name or ...:
            if __setGrib2TableParameters__:
                # This user defined parameters must be at last of this function!
                for key, val in __setGrib2TableParameters__:
                    gribapi.grib_set_long(grib_message, key, val)
                    print "set user defined grib2table parameter ('%s', %s)" % (key, val)
            # end of if __setGrib2TableParameters__:
            print "Tweaking end ", cube.standard_name
            
            yield grib_message
        # end of for cube, grib_message in iris.fileformats.grib.as_pairs(cube):
    # end of for cube in cubeList:
# end of def tweaked_messages(cube):

# Start the convertFilesInParallel function
def convertFilesInParallel(fname, fext, ftype):
    """
    convertFilesInParallel function calling all the sub-functions
    :param fnames: a simple filename as argument in a string format
    """
    
    global _startT_, _tmpDir_, _opPath_, __max_long_fcst_hours__,\
           __start_step_long_fcst_hour__, _createGrib2CtlIdxFiles_

    if __start_step_long_fcst_hour__ == 6:
        # here max fcst hours goes upto 240 only, not 241. why ??
        # because 216 long fcst hours contains upto 240th hour fcst.
        # and 240th long fcst contains upto 264th hour fcst.
        # so here no need to add +1 to __max_long_fcst_hours__.
        fcst_times = [str(hr).zfill(3) for hr in range(0, __max_long_fcst_hours__, 24)]
    elif __start_step_long_fcst_hour__ == 24:
        fcst_times = ['000', ]
    # end of if __start_step_long_fcst_hour__ == 6:
    
    fcst_filenames = [(fname, fext, hr) for hr in fcst_times]
    print fcst_filenames
    ## get the no of files and 
    nprocesses = len(fcst_filenames)
    # lets create no of parallel process w.r.t no of files.
    
    # parallel begin - 1 
    pool = _MyPool(nprocesses)
    print "Creating %d (non-daemon) workers and jobs in convertFilesInParallel process." % nprocesses        
    if ftype in ['fcst', 'forecast']:        
        results = pool.map(packEnsemblesInParallel, fcst_filenames)
    else:
        raise ValueError("Unknown file type !")
    # end of if ftype in ['anl', 'analysis']:    

    # closing and joining master pools
    pool.close()     
    pool.join()
    # parallel end - 1     
    print "Total time taken to convert %d files was: %8.5f seconds \n" %(len(fcst_filenames),(time.time()-_startT_))
    
    return
# end of def convertFilesInParallel(fnames):

def doShuffleVarsInOrder(fpath):
    """
    order the variables and create new grib2 files;
    delete the older shuffled variables grib2 files.
    create ctl, idx files using g2ctl.pl, gribmap scripts for the ordered grib2 files.
    
    Arulalan/T
    11-12-2015
    """
    global  _orderedVars_, _preExtension_, _ncfilesVars_, _inDataPath_, \
           _maskOverOceanVars_, _aod_pseudo_level_var_, _createGrib2CtlIdxFiles_, \
           _createGrib1CtlIdxFiles_, _convertGrib2FilestoGrib1Files_, \
           _convertVars_, __outFileType__, __grib1FilesNameSuffix__, \
           __removeGrib2FilesAfterGrib1FilesCreated__, _removeVars_, \
           __start_step_long_fcst_hour__, g2ctl, grib2ctl, gribmap, cnvgrib
    
    print "re-order take place: going to access ", fpath
    try:        
        f = iris.load(fpath)
    except gribapi.GribInternalError as e:
        if str(e) == "Wrong message length":
            print "ALERT!!!! ERROR!!! Couldn't read grib2 file to re-order", e
        else:
            print "ALERT!!! ERROR!!! couldn't read grib2 file to re-order", e
        return 
    except Exception as e:
        print "ALERT!!! ERROR!!! couldn't read grib2 file to re-order", e
        return 
    # end of try:
    # get only the pressure coordinate variables
    unOrderedPressureLevelVarsList = [i for i in f if len(i.coords('pressure')) == 1]
    # get only the non pressure coordinate variables
    unOrderedNonPressureLevelVarsList = list(set(f) - set(unOrderedPressureLevelVarsList))
    
    # generate dictionary (standard_name, STASH) as key and cube variable as value
    unOrderedPressureLevelVars = {}
    for i in unOrderedPressureLevelVarsList:
        name = i.standard_name if i.standard_name else i.long_name
        unOrderedPressureLevelVars[name] = i
        
    unOrderedNonPressureLevelVars = {}
    for i in unOrderedNonPressureLevelVarsList:
        name = i.standard_name if i.standard_name else i.long_name
        unOrderedNonPressureLevelVars[name] = i   
    
    # need to store the ordered variables in this empty list
    orderedVars = []    
    if _convertVars_:
        # user has passed their own ordered and limited vars 
        orderedVarsList = _convertVars_
    else:
        # use inbuilt ordered list from this module itself
        orderedVarsList = _orderedVars_['PressureLevel'] + _orderedVars_['nonPressureLevel']
        
    for (varName, varSTASH) in orderedVarsList:
        # skip if user specified var not in pressure level vars list 
        if not (varName, varSTASH) in _orderedVars_['PressureLevel']: continue
        # got pressure vars, add to ordered final vars list  
        if varName in unOrderedPressureLevelVars: orderedVars.append(unOrderedPressureLevelVars[varName])
    # end of for name, STASH in _orderedVars_['PressureLevel']:
            
    ncloaddic = {}
    ncloadedfiles = []
    for (varName, varSTASH) in orderedVarsList:
        # skip if user specified var not in non-pressure level vars list 
        if not (varName, varSTASH) in _orderedVars_['nonPressureLevel']: continue
        # got non-pressure vars, add to ordered final vars list  
        if varName in unOrderedNonPressureLevelVars: 
            orderedVars.append(unOrderedNonPressureLevelVars[varName])
        elif (varName, varSTASH) in _ncfilesVars_:            
            ## generate nc file name 
            ncfpath = varSTASH + '_' + '.'.join(fpath.split('.')[:-1]) + '.nc'
            if not os.path.isfile(ncfpath): continue
            if varSTASH not in ncloaddic:
                try:
                    ncloaddic['varSTASH'] = iris.load(ncfpath)
                    ncloadedfiles.append(ncfpath)
                except Exception as e:
                    print "ALERT!!! ERROR!!! couldn't read nc file to re-order", e
                    return 
                # end of try:             
            # end of if varSTASH not in ncloaddic:
            # define variable name constraint
            varConstraint = iris.Constraint(name=varName)
            # define varibale stash code constraint
            STASHConstraint = iris.AttributeConstraint(STASH=varSTASH)
            var = ncloaddic['varSTASH'].extract(varConstraint & STASHConstraint)
            
            if var: 
                # apped the ordered / corrected vars into the list, which will  
                # be going to saved into grib2 files by tweaking it further!                
                if varName in _aod_pseudo_level_var_:
                    for plev, pval in _aod_pseudo_level_var_[varName]:
                        pvar = var[0].extract(iris.Constraint(pseudo_level=plev))
                        # set standard_name as None
                        pvar.standard_name = None
                        # set long_name as standard_name + '_at_micronwavelength'
                        # so that _grib_cf_map will be able to identify local 
                        # table grib2 param code. 
                        pvar.long_name = varName + '_at_%sum' % pval
                        orderedVars.append(pvar)
                else:
                    orderedVars.append(var[0])
            # end of if var:
    # end of for (varName, STASH) in orderedVarsList:
    
    # store the land_binary_mask data into temporary variable
    land_binary_mask_var = [var for var in orderedVars 
                            if var.standard_name == 'land_binary_mask']
                                
    if _maskOverOceanVars_ and land_binary_mask_var:
        
        land_binary_mask = land_binary_mask_var[0].data < 1
        # here we are masking less than 1. we can do just simply == 0 also, 
        # but somehow it retains fraction values between 0 to 1. To get 
        # ride out of this fraction values, just mask out < 1.
        # get the shapes
        lsh = land_binary_mask.shape
        # Define constraint to extract latitude from 60S to 60N
        lat_60S_60N = iris.Constraint(latitude=lambda cell: -60 < cell < 60)
        lat_30S_30N = iris.Constraint(latitude=lambda cell: -30 < cell < 30)
        for vidx, var in enumerate(orderedVars):
            vname = var.standard_name if var.standard_name else var.long_name
            if vname in _maskOverOceanVars_:    
                # Lets reset zero values lies within 60S to 60N band
                # with 0.0051, before ocean region has been masked.
                # Now extract data only lies between 60S to 60N
                var_60S_60N = var.extract(lat_60S_60N)
                if vname == 'volumetric_moisture_of_soil_layer':
                    # reset the minimum values as 0.0051
                    var_60S_60N.data[var_60S_60N.data < 0.005] = 0.0051
                elif vname == 'soil_temperature':
                    # We should assign min of extra tropical band !
                    # because polar minimum might have been assigned 
                    # while extracting data (function regridAnlFcstFiles).
                    # So lets re-set here!
                    # zero will not make sense when temperature unit is Kelvin
                    # get tropical data 
                    var_30S_30N = var.extract(lat_30S_30N) 
                    # find next mean value of tropical data 
                    nmean = numpy.ma.masked_less_equal(var_30S_30N.data, var_30S_30N.data.min()).mean()
                    # make memory free 
                    del var_30S_30N
                    # set mean value (of tropical data 30S to 30N) as 
                    # min value full extra tropical data (60S to 60N).
                    # This will solve the abnormal temperature values over 
                    # small islands in tropical ocean regions.
                    var_60S_60N.data[var_60S_60N.data <= var_60S_60N.data.min()] = nmean
                # end of if vname == 'volumetric_moisture_of_soil_layer':
                
                # extract latitude coords of 60S to 60N data 
                lat_60S_60N_points = var_60S_60N.coords('latitude')[0].points
                # get its start and end lat values of 60S and 60N
                lat_60S_start, lat_60N_end = lat_60S_60N_points[0], lat_60S_60N_points[-1]
                # get the original global data latitude points 
                originalLat = var.coords('latitude')[0].points.tolist()
                # find the 60S index in original global latitude
                lat_60S_index = originalLat.index(lat_60S_start)
                # find the 60N index in original global latitude
                lat_60N_index = originalLat.index(lat_60N_end) + 1
                # Lets insert the updated data (0.0051) within the 
                # original global data itself.
                # Now Lets do masking over Ocean regions!
                vsh = var.shape
                if lsh != vsh:
                    # first dimension points 4 layer depth_below_land_surface
                    # so second dimension points latitude.
                    var.data[:, lat_60S_index: lat_60N_index, :] = var_60S_60N.data                        
                    # get the ocean mask by masking 0s of land_binary_mask 
                    # (0-sea, 1-land) and set it to the required variables. 
                    land_binary_mask_grown = land_binary_mask.reshape(1, lsh[0], lsh[-1])
                    land_binary_mask_grown = land_binary_mask_grown.repeat(vsh[0], axis=0)                    
                    var.data = numpy.ma.masked_where(land_binary_mask_grown, var.data)
                else:    
                    # single layer only. so first dimension points latitude
                    var.data[lat_60S_index: lat_60N_index, :] = var_60S_60N.data
                    # get the ocean mask by masking 0s of land_binary_mask 
                    # (0-sea, 1-land) and set it to the required variables. 
                    var.data = numpy.ma.masked_where(land_binary_mask, var.data)
                # end of if lsh != vsh:
                print "updated ocean masked vars"
        # end of for vidx, var in enumerate(orderedVars):
    # end of if _maskOverOceanVars_:
    
    # removing land_binary_mask_var from out files if it is forecast grib2 file
    # why do we need to repeat the same static variables in all the 
    # forecast files... So removing it, but keeps in analysis file.
    if __outFileType__ in ['prg', 'fcst'] and land_binary_mask_var and \
                                  __start_step_long_fcst_hour__ in [6]: 
        # remove only for 6 hourly ncum post prodction. Not for others!
        # say for 3 hourly hycom model input landsea binary mask needed in all
        # forecast files.
        orderedVars.remove(land_binary_mask_var[0])
    # But still we have to use land_binary_mask variable to set 
    # ocean mask for the soil variables. Thats why we included it in vars list.
    
    oidx = None
    if ('surface_upwelling_shortwave_flux_in_air', 'None') in _convertVars_:
        # find the index in _convertVars_
        idx = _convertVars_.index(('surface_upwelling_shortwave_flux_in_air', 'None'))
        # adjust the current index by subtract 1, because in previous insertion 
        # causes order index increased by 1.
        idx = idx-1 if (idx and oidx is None) else idx
        oidx = idx
        # store the surface_net_downward_shortwave_flux data into temporary variable
        surface_net_downward_shortwave_flux = [var for var in orderedVars 
               if var.standard_name == 'surface_net_downward_shortwave_flux']    
        if not surface_net_downward_shortwave_flux:
            raise ValueError("Can not calculate surface_upwelling_shortwave_flux, because unable to load surface_net_downward_shortwave")    
        # store the surface_downwelling_shortwave_flux_in_air data into temporary variable
        surface_downwelling_shortwave_flux = [var for var in orderedVars 
               if var.standard_name == 'surface_downwelling_shortwave_flux_in_air']
        if not surface_downwelling_shortwave_flux:
            raise ValueError("Can not calculate surface_upwelling_shortwave_flux, because unable to load surface_downwelling_shortwave_flux")
        # calculate 'surface_upwelling_shortwave_flux' by subtract 'surface_net_downward_shortwave_flux'
        # from 'surface_downwelling_shortwave_flux'       
        surface_upwelling_shortwave_flux = cubeSubtractor(surface_downwelling_shortwave_flux[0], 
                                           surface_net_downward_shortwave_flux[0], 
                       standard_name='surface_upwelling_shortwave_flux_in_air',
                                                              removeSTASH=True)
        # store the 'surface_upwelling_shortwave_flux' into orderedVars
        orderedVars.insert(idx, surface_upwelling_shortwave_flux)
    # end of if ('surface_upwelling_shortwave_flux_in_air', 'None') in _convertVars_:
    
    if ('surface_upwelling_longwave_flux_in_air', 'None') in _convertVars_:
        # find the index in _convertVars_
        idx = _convertVars_.index(('surface_upwelling_longwave_flux_in_air', 'None'))
        # adjust the current index by subtract 1, because in previous insertion 
        # causes order index increased by 1.
        idx = idx-1 if (idx and oidx is None) else idx
        oidx = idx
        # store the surface_net_downward_longwave_flux data into temporary variable
        surface_net_downward_longwave_flux = [var for var in orderedVars 
               if var.standard_name == 'surface_net_downward_longwave_flux'] 
        if not surface_net_downward_longwave_flux:
            raise ValueError("Can not calculate surface_upwelling_longwave_flux, because unable to load surface_downwelling_longwave_flux")           
        # store the surface_downwelling_longwave_flux data into temporary variable
        surface_downwelling_longwave_flux = [var for var in orderedVars 
               if var.standard_name == 'surface_downwelling_longwave_flux']
        if not surface_downwelling_longwave_flux:
            raise ValueError("Can not calculate surface_upwelling_longwave_flux, because unable to load surface_downwelling_shortwave_flux")
        # calculate 'surface_upwelling_longwave_flux' by subtract 'surface_net_downward_longwave_flux'
        # from 'surface_downwelling_longwave_flux'       
        surface_upwelling_longwave_flux = cubeSubtractor(surface_downwelling_longwave_flux[0], 
                                           surface_net_downward_longwave_flux[0], 
                       standard_name='surface_upwelling_longwave_flux_in_air',
                                                              removeSTASH=True)
        # store the 'surface_upwelling_longwave_flux' into orderedVars
        orderedVars.insert(idx, surface_upwelling_longwave_flux)
    # end of if ('surface_upwelling_longwave_flux_in_air', 'None') in _convertVars_:
    
    # remove temporary variables from ordered vars list 
    for dvar, dSTASH in _removeVars_:
        for ovar in orderedVars:
            # removed temporary vars from ordered vars list        
            if ovar.standard_name == dvar: orderedVars.remove(ovar)
    # end of for dvar in _removeVars_:

    # generate correct file name by removing _preExtension_
    g2filepath = fpath.split(_preExtension_)
    g2filepath = g2filepath[0] + g2filepath[-1]
    
    # now lets save the ordered variables into same file
    try:   
        # before save it, tweak the cubes by setting centre no and 
        # address other temporary issues before saving into grib2.
        iris.fileformats.grib.save_messages(tweaked_messages(orderedVars), 
                                                g2filepath, append=True)
    except Exception as e:
        print "ALERT !!! Error while saving orderd variables into grib2!! %s" % str(e)
        print " So skipping this without saving data"
        return 
    # end of try:
    
    # make memory free 
    del orderedVars
    
    # remove the older file 
    os.remove(fpath)
    for ncf in ncloadedfiles: os.remove(ncf)
        
    print "Created the variables in ordered fassion and saved into", g2filepath
                
    if _convertGrib2FilestoGrib1Files_:
        g1filepath = '.'.join(g2filepath.split('.')[:-1])
        g1filepath = g1filepath if g1filepath else g2filepath[:-1]
        if __grib1FilesNameSuffix__: g1filepath += str(__grib1FilesNameSuffix__)
        
        if os.path.isfile(g1filepath): os.remove(g1filepath)
        
        cmd = [cnvgrib, '-g21', g2filepath, g1filepath]
        subprocess.call(cmd, shell=False)
        cmd = ['chmod', '644', g1filepath]
        subprocess.call(cmd, shell=False)
        print "Converted grib2 to grib1 file : -", g1filepath
        
        if  _createGrib1CtlIdxFiles_:
            ## grib2ctl.pl usage option refer the below link 
            ## https://tuxcoder.wordpress.com/2011/04/11/how-to-install-grib2ctl-pl-and-wgrib-in-linux/    
            ctlfile = open(g1filepath+'.ctl', 'w')
            if __outFileType__ in ['ana', 'anl']:
                # create ctl & idx files for analysis file 
                subprocess.call([grib2ctl, '-ts6hr', g1filepath], stdout=ctlfile)
                subprocess.call([gribmap, '-ts6hr', '-0', '-i', g1filepath+'.ctl'])
            elif __outFileType__ in ['prg', 'fcst']:
                # create ctl & idx files for forecast file
                subprocess.call([grib2ctl, '-ts6hr', '-verf', g1filepath], stdout=ctlfile)
                subprocess.call([gribmap, '-i', g1filepath+'.ctl'])
            else:
                raise ValueError("unknown file type while executing grib2ctl.pl!!")
            
            print "Successfully created control and index file using grib2ctl !", g1filepath+'.ctl'
    # end of if _createGrib1CtlIdxFiles_:
    # end of if _convertGrib2FilestoGrib1Files_:
    
    if __removeGrib2FilesAfterGrib1FilesCreated__ and _convertGrib2FilestoGrib1Files_:
        # grib1 files are converted. so we can remove grib2 files.
        os.remove(g2filepath)
        print "deleted grib2 file", g2filepath        
    elif _createGrib2CtlIdxFiles_:
        ## g2ctl.pl usage option refer the below link 
        ## https://tuxcoder.wordpress.com/2011/08/31/how-to-install-g2ctl-pl-and-wgrib2-in-linux/
        ## though options says -verf for forecast end time, -0 for analysis time 
        ## -b for forecast start time, its all about setting reference in ctl file.
        ## Nothing more than that. We already set correct reference and forecast time bounds 
        ## in analysis files (whichever variables are actually taken from previous short forecast 0-6 hours).
        ## so here we no need to pass any options like -0 or -b. 
        ## By default g2ctl takes -verf option, same option we are passing 
        ## here to make sure that in future it will not affect.
        ctlfile = open(g2filepath+'.ctl', 'w')
        # create ctl & idx files for forecast file
        # by default -verf as passed which takes end time of fcst bounds to set as base time.
        subprocess.call([g2ctl, '-ts6hr', '-verf', g2filepath], stdout=ctlfile)
        subprocess.call([gribmap, '-i', g2filepath+'.ctl'])                
        print "Successfully created control and index file using g2ctl !", g2filepath+'.ctl'
    # end of if __removeGrib2FilesAfterGrib1FilesCreated__:    
# end of def doShuffleVarsInOrder(fpath):

def doShuffleVarsInOrderInParallel(ftype, simulated_hr):
            
    global _current_date_, _opPath_, _preExtension_, __max_long_fcst_hours__, \
           __anlFileNameStructure__, __fcstFileNameStructure__, \
           __max_long_fcst_hours__, __start_step_long_fcst_hour__, __utc__
            
    print "Lets re-order variables for all the files!!!"
    #####
    ## 6-hourly Files have been created with extension.
    ## Now lets do re-order variables within those individual files, in parallel mode.  
            
    # get current working directory
    current_dir = os.getcwd()
    # lets change current working directory as out path
    os.chdir(_opPath_)
    if ftype in ['fcst', 'forecast']:
        ## generate all the forecast filenames w.r.t forecast hours 
        # get the out fileName Structure based on pre / user defined indecies                       
        outFnIndecies = __getAnlFcstFileNameIdecies__(__fcstFileNameStructure__)
        fcstFiles = []
        for fcsthr in range(__start_step_long_fcst_hour__, 
                   __max_long_fcst_hours__+1, __start_step_long_fcst_hour__):            
            # generate the out file name based on actual informations                                 
            outFn = __genAnlFcstOutFileName__(__fcstFileNameStructure__, 
                                  outFnIndecies, _current_date_, fcsthr, 
                                           simulated_hr, _preExtension_)  
            fcstFiles.append(outFn)
        # end of for fcsthr in range(...):

        ## get the no of created fcst files  
        nprocesses = len(fcstFiles)        
        # parallel begin - 3
        pool = _MyPool(nprocesses)
        print "Creating %d (non-daemon) workers and jobs in doShuffleVarsInOrder process." % nprocesses
        results = pool.map(doShuffleVarsInOrder, fcstFiles)   
        
        # closing and joining master pools
        pool.close()     
        pool.join()
        # parallel end - 3    
    elif ftype in ['anl', 'analysis']:
        ## generate the analysis filename w.r.t simulated_hr
        # get the out fileName Structure based on pre / user defined indecies                       
        outFnIndecies = __getAnlFcstFileNameIdecies__(__anlFileNameStructure__) 
        anlFiles = []
        simulated_hr = int(__utc__)
        # since ncum producing analysis files 00, 06, 12, 18 utc cycles and 
        # its forecast time starting from 0 and reference time based on utc.
        # so we should calculate correct hour as below.
        for fcsthr in range(0+simulated_hr, 6+simulated_hr, __anl_step_hour__):            
            # generate the out file name based on actual informations                                 
            outFn = __genAnlFcstOutFileName__(__anlFileNameStructure__, 
                                  outFnIndecies, _current_date_, fcsthr, 
                                           simulated_hr, _preExtension_)  
            anlFiles.append(outFn)
        # end of for fcsthr in range(...):
        ## get the no of created anl files  
        nprocesses = len(anlFiles)        
        # parallel begin - 3
        pool = _MyPool(nprocesses)
        print "Creating %d (non-daemon) workers and jobs in doShuffleVarsInOrder process." % nprocesses
        results = pool.map(doShuffleVarsInOrder, anlFiles)   
        
        # closing and joining master pools
        pool.close()     
        pool.join()
        # parallel end - 3        
    # end of if ftype in ['fcst', 'forecast']: 
    print "Total time taken to convert and re-order all files was: %8.5f seconds \n" % (time.time()-_startT_)
    # reset current working directory
    os.chdir(current_dir)
    return 
# end of def doShuffleVarsInOrderInParallel(arg):
    
def _checkInFilesStatus(path, ftype, pfname, fext):
    
    global __max_long_fcst_hours__, _ensemble_count_
    
    if ftype in ['ana', 'anl']:
        fhrs = ['000'] 
    elif ftype in ['fcst', 'prg']:
        # here max fcst hours goes upto 240 only, not 241. why ??
        # because 216 long fcst hours contains upto 240th hour fcst.
        # and 240th long fcst contains upto 264th hour fcst.
        # so here no need to add +1 to __max_long_fcst_hours__.
        fhrs = [str(hr).zfill(3) for hr in range(0, __max_long_fcst_hours__, 24)]
    
    fileNotExistList = []

    for ehr in range(0, _ensemble_count_+1, 1):
        if fext and pfname == 'xbfti':
            # eg : xbfti_044.pp0
            fname = pfname + '_' + str(ehr).zfill(3) + fext
            fpath = os.path.join(path, fname)
            if not os.path.isfile(fpath): fileNotExistList.append(fpath)
        elif pfname == 'pb':
            for fhr in fhrs:
                # construct the correct fileName from partial fileName and hours
                # add hour only if doenst have any extension on partial filename.
                fname = str(ehr).zfill(3) + '_' + pfname + fhr
                fpath = os.path.join(path, fname)
                if not os.path.isfile(fpath): fileNotExistList.append(fpath)
    # end of for ehr in range(0, _ensemble_count_+1, 1):
    status = False if fileNotExistList else True
    if status is False:    
        print "The following infiles are not exists!\n"
        print "*" * 80
        print "\n".join(fileNotExistList)
        print "*" * 80
        
    return status
# end of def _checkInFilesStatus(path, ftype, pfnames):

def _checkOutFilesStatus(path, ftype, date, utc, overwrite):
    
    global _preExtension_, __max_long_fcst_hours__, __anlFileNameStructure__,\
           __fcstFileNameStructure__, __start_step_long_fcst_hour__, \
           __anl_step_hour__, __utc__
           
    if ftype in ['ana', 'anl']:
        outFileNameStructure = __anlFileNameStructure__
        fhrs = [utc] # ana_hour (short forecast hour) is same as simulated_hr (i.e. utc)
        simulated_hr = int(__utc__)
        # since ncum producing analysis files 00, 06, 12, 18 utc cycles and 
        # its forecast time starting from 0 and reference time based on utc.
        # so we should calculate correct hour as below.
        fhrs = range(0+simulated_hr, 6+simulated_hr, __anl_step_hour__)     
    elif ftype in ['fcst', 'prg']:
        outFileNameStructure = __fcstFileNameStructure__
        fhrs = range(__start_step_long_fcst_hour__, __max_long_fcst_hours__+1, 
                                                __start_step_long_fcst_hour__)
    
    # get the out fileName Structure based on pre / user defined indecies
    outFnIndecies = __getAnlFcstFileNameIdecies__(outFileNameStructure)
    status = None
    for fhr in fhrs:
        # generate the out file name based on actual informations.
        # here preExtension is empty string to create final needed out file name                        
        fname = __genAnlFcstOutFileName__(outFileNameStructure, outFnIndecies,  
                                                               date, fhr, utc)
        fpath = os.path.join(path, fname)        
        for ext in ['', '.ctl', '.idx']:
            fpath = os.path.join(path, fname+ext)
            if os.path.isfile(fpath):
                print "Out File already exists", fpath,
                if overwrite: 
                    os.remove(fpath)
                    status = 'FilesRemoved'
                    print ", but overwrite option is True. So removed it!"
                else:
                    status = 'FilesExist' 
            else:
                print "\nOut File does not exists", fpath 
                if status in [None, 'FilesRemoved']:
                    status = 'FilesDoNotExist' 
                    continue
                elif status is 'FilesExist':
                    status = 'PartialFilesExist'
                    break
        # end of for ext in ['', '.ctl', '.idx']:
    # end of for fhr in fhrs:
    
    ifiles = [fname for fname in os.listdir(path) if _preExtension_ in fname]
    if ifiles:        
        print "Intermediate files are exists in the outdirectory.", path
        for ifile in ifiles:        
            if outFileNameStructure[0] in ifile and utc in ifile and _preExtension_ in ifile:
                os.remove(os.path.join(path, ifile))
                status = 'IntermediateFilesExist'
                print "removed intermediate nc file"                
        # end of for ncfile in ncfiles:        
    # end of if ncfiles:
    if status in ['PartialFilesExist', 'IntermediateFilesExist']:
        # partial files exist, so make overwrite option as True and do 
        # recursive call one time to remove all output files.
        print "Partial/Intermediate out files exist, so going to overwrite all files"
        return _checkOutFilesStatus(path, ftype, date, utc, overwrite=True)
    else:
        return status
# end of def _checkOutFilesStatus(path, ftype, date, hr, overwrite):

def convertFcstFiles(inPath, outPath, tmpPath, **kwarg):
           
    global _targetGrid_, _targetGridRes_, _current_date_, _startT_, _tmpDir_, \
       _inDataPath_, _opPath_, _doRegrid_, _convertVars_, _requiredLat_, \
       _requiredLon_, _createGrib2CtlIdxFiles_, _createGrib1CtlIdxFiles_, \
       _convertGrib2FilestoGrib1Files_, __fcstFileNameStructure__, \
       __LPRINT__, __utc__, __start_step_long_fcst_hour__, \
       __max_long_fcst_hours__, __outFileType__, __grib1FilesNameSuffix__, \
       __removeGrib2FilesAfterGrib1FilesCreated__, _depedendantVars_, \
       _removeVars_, _requiredPressureLevels_, __setGrib2TableParameters__, \
        __outg2files__
     
    # load key word arguments
    targetGridResolution = kwarg.get('targetGridResolution', 0.25)
    date = kwarg.get('date', time.strftime('%Y%m%d'))
    utc = kwarg.get('utc', '00')
    overwrite = kwarg.get('overwrite', False)
    lprint = kwarg.get('lprint', False)
    convertVars = kwarg.get('convertVars', None)
    latitude = kwarg.get('latitude', None)
    longitude = kwarg.get('longitude', None)
    pressureLevels = kwarg.get('pressureLevels', None)
    start_step_long_fcst_hour = kwarg.get('start_step_long_fcst_hour', 6)
    max_long_fcst_hours = kwarg.get('max_long_fcst_hours', 240)
    fcstFileNameStructure = kwarg.get('fcstFileNameStructure', None)
    createGrib2CtlIdxFiles = kwarg.get('createGrib2CtlIdxFiles', True)
    createGrib1CtlIdxFiles = kwarg.get('createGrib1CtlIdxFiles', False)
    convertGrib2FilestoGrib1Files = kwarg.get('convertGrib2FilestoGrib1Files', False)
    grib1FilesNameSuffix = kwarg.get('grib1FilesNameSuffix', '1')
    removeGrib2FilesAfterGrib1FilesCreated = kwarg.get('removeGrib2FilesAfterGrib1FilesCreated', False)
    callBackScript = kwarg.get('callBackScript', None)
    setGrib2TableParameters = kwarg.get('setGrib2TableParameters', None)
    
    # assign out file type in global variable
    __outFileType__ = 'fcst'
    # assign the convert vars list of tuples to global variable
    if convertVars: _convertVars_ = convertVars
    # assign the analysis file name structure
    if fcstFileNameStructure: __fcstFileNameStructure__ = fcstFileNameStructure
    # set print variables details options
    __LPRINT__ = lprint    
    # update global variables
    __utc__ = utc
    __start_step_long_fcst_hour__ = start_step_long_fcst_hour
    __max_long_fcst_hours__ = max_long_fcst_hours
    __removeGrib2FilesAfterGrib1FilesCreated__ = removeGrib2FilesAfterGrib1FilesCreated
    __grib1FilesNameSuffix__ = grib1FilesNameSuffix
    _targetGridRes_ = str(targetGridResolution)
    _requiredLat_ = latitude
    _requiredLon_ = longitude
    _requiredPressureLevels_ = pressureLevels    
    _createGrib2CtlIdxFiles_ = createGrib2CtlIdxFiles
    _createGrib1CtlIdxFiles_ = createGrib1CtlIdxFiles
    _convertGrib2FilestoGrib1Files_ = convertGrib2FilestoGrib1Files
    __setGrib2TableParameters__ = setGrib2TableParameters
    # forecast filenames partial name
    if __start_step_long_fcst_hour__ == 6:
        fcst_fname = 'pb'
        fext = ['',]
    elif __start_step_long_fcst_hour__ == 24:
        fcst_fname = 'xbfti'
        fext = ['.pp0', '.pp2']        
        
    # get the current date in YYYYMMDD format
    _tmpDir_ = tmpPath
    _current_date_ = date
    print "\n _current_date_ is %s" % _current_date_
    logpath = os.path.join(_tmpDir_, _current_date_)
    createDirWhileParallelRacing(logpath)
    logfile = 'um2grb2_fcst_stdout_'+ _current_date_ +'_' + utc +'Z.log'
    sys.stdout = myLog(os.path.join(logpath, logfile))
    
    # start the timer now
    _startT_ = time.time()

    # set-up base folders    
    _inDataPath_ = os.path.join(inPath, _current_date_)
    if not os.path.exists(_inDataPath_):
        raise ValueError("In datapath does not exists %s" % _inDataPath_)
    # end of if not os.path.exists(_inDataPath_):
    
    if convertVars:
        # check either depedendant vars are need to be loaded 
        for var, dvars in _depedendantVars_.iteritems():
            if var in convertVars:
                for dvar in dvars:
                    if dvar not in convertVars:               
                        _convertVars_.append(dvar)  # include depedendant var
                        _removeVars_.append(dvar)   # remove depedendant var at last
        # end of for var, dvar in _depedendantVars_.iteritems():
        if fcst_fname == 'xbfti':            
            # load only required file names to avoid unnneccessary computations
            # by cross checking with user defined variables list.
            for ext in fext[:]:   
                # loop through copy of fext[:], because fext list 
                # will may change within this loop.
                hr = utc.zfill(3)
                ## if fileName has some extension, then do not add hr to it.
                fileName = fcst_fname + '_' + hr + ext
                varNamesSTASH, _, _, _ = getVarInOutFilesDetails(_inDataPath_, fileName, hr)
                # check either user requires this file or not!
                if not set(varNamesSTASH).intersection(convertVars):
                    # remove the ext from fext, because user didn't 
                    # require variabels from this ext file.
                    fext.remove(ext)
                    print "removed %s from list of files" % ext 
            print "Final fext list :", fext
        # end of if fcst_fname == 'xbfti':            
    # end of if convertVars:    
    
    for ext in fext:
        # check either infiles are exist or not!
        status = _checkInFilesStatus(_inDataPath_, 'prg', fcst_fname, ext)
        print "in status+++++++++++++++++++++++++++", status
        if not status:
            raise ValueError("In datapath does not contain the above valid infiles")
        # end of if not instatus:
    # end of for ext in fext:
    
    _opPath_ = os.path.join(outPath, _current_date_)
    createDirWhileParallelRacing(_opPath_) 
        
    if targetGridResolution is None:
        _doRegrid_ = False  
    else:
        if not isinstance(targetGridResolution, (int, float)):
            raise ValueError("targetGridResolution must be either int or float")
        # define default global lat start, lon end points
        slat, elat = (-90., 90.)
        # define default global lon start, lon end points 
        slon, elon = (0., 360.)
        # define user defined custom lat & lon start and end points
        if latitude: (slat, elat) = latitude
        if longitude: (slon, elon) = longitude
        # target grid as 0.25 deg (default) resolution by setting up sample points 
        # based on coord    
        _targetGrid_ = [('latitude', numpy.arange(slat, 
                          elat+targetGridResolution, targetGridResolution)),
                        ('longitude', numpy.arange(slon, 
                          elon+targetGridResolution, targetGridResolution))]
        _doRegrid_ = True  
    # end of if targetGridResolution is None:
    
    # check either files are exists or not. delete the existing files in case
    # of overwrite option is True, else return without re-converting files.
    status = _checkOutFilesStatus(_opPath_, 'prg', _current_date_, utc, overwrite)
    if status is 'FilesExist': 
        print "All files are already exists. So skipping convert Fcst files porcess"
        return # return back without executing conversion process.
    elif status in [None, 'FilesDoNotExist', 'FilesRemoved']:
        print "Going to start convert Fcst files freshly"
    # end of if status is 'FilesExists': 
    
    # do convert for forecast files     
    for ext in fext: convertFilesInParallel(fcst_fname, ext, ftype='fcst')  
    
    # do re-order variables within files in parallel
    doShuffleVarsInOrderInParallel('fcst', utc)
    
    if callBackScript:
        callBackScript = os.path.abspath(callBackScript)
        if not os.path.exists(callBackScript): 
            print "callBackScript '%s' doenst exist" % callBackScript
            return 
        kwargs = ' --date=%s --outpath=%s --oftype=forecast --utc=%s' % (_current_date_, _opPath_, utc)
        scriptExecuteCmd = callBackScript + ' ' + kwargs
        # execute user defined call back script with keyword arguments
        subprocess.call(scriptExecuteCmd, shell=True)
    # end of if callBackScript:
# end of def convertFcstFiles(...):



# -- End code
