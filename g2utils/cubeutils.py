import numpy, iris  

def cubeAverager(tmpCube, action='mean', dt='1 hour', 
                actionIntervals='6 hour', tpoint='cbound', fpoint='cbound'):
    """
    :param tmpCube:     The temporary cube data (in Iris format) with non-singleton time dimension
    :param action:      mean| sum (accumulated fields are summed and instantaneous are averaged).
    :param dt:   A standard string representing forecast step duration/intervals.
    :param actionIntervals: A non standard string to add inside cell_methods comments section.
    :param tpoint: cbound | lbound | rbound time point 
    :param fpoint: cbound | lbound | rbound forecast period
    :return: meanCube:  An Iris formatted cube date containing the resultant data either as
                        averaged or summed.
    Arulalan.T
    16-Nov-2015    
    """

    # extract time points 
    tpoints = tmpCube.coord('forecast_period').points
    # assign first time point data to mean data 
    meanCube = tmpCube.extract(iris.Constraint(forecast_period=tpoints[0]))
    # get the time coord of first time cube and set to mean
    timeAxFirst = meanCube.coords('time')[0]
    # get the fcst time coord of first time cube and set to mean
    fcstAxFirst = meanCube.coords('forecast_period')[0]
    
    # loop through remaining time points 
    for tp in tpoints[1:]:
        # extract remaining time points and add to meanCube
        lastCube = tmpCube.extract(iris.Constraint(forecast_period=tp))
        meanCube = iris.analysis.maths.add(meanCube, lastCube) 
    # end of for tp in tpoints[1:]:
    
    # get the time coord of last time cube and set to mean
    timeAxLast = lastCube.coords('time')[0]
    # get the fcst time coord of last time cube and set to mean
    fcstAxLast = tmpCube[-1].coords('forecast_period')[0]
    
    if action == 'mean':
        # to compute mean value of meanCube divide it by length of time
        # points which we added before 
        meanCube = iris.analysis.maths.divide(meanCube, float(len(tpoints)))
        print "Converted cube to %s mean : %s" % (actionIntervals, tmpCube.standard_name)
    elif action == 'sum':
        # for sum action, we no need to do here anything, because already 
        # we computed accumulation only!
        print "Converted cube to %s accumulation : %s" % (actionIntervals, tmpCube.standard_name)
    else:
        raise ValueError('argument "%s" not support' % action)
    # end of if not isAccumulation:

    # get the bounds and time points from two extremes    
    bounds = [timeAxFirst.bounds[0][0], timeAxLast.bounds[-1][-1]]
    
    if tpoint == 'cbound':
        #### THE CENTRE POINT OF REFERENCE TIME BOUNDS
        timepoint = [bounds[-1] + ((bounds[-1] - bounds[0]) / 2.0)]
    elif tpoint == 'lbound':
       ###### THE START BOUNDS OF REFERENCE TIME POINT
        timepoint = [bounds[0]]        
    elif tpoint == 'rbound':    
        ###### THE END BOUNDS OF REFERENCE TIME POINT
        timepoint = [bounds[-1]]   
    # end of if tpoint == 'cbound':
    
    # update the time coordinate with new time point and time bounds
    timeAxFirst.points = timepoint
    timeAxFirst.bounds = bounds
    # add the updated time coordinate to the meanCube
    meanCube.add_aux_coord(timeAxFirst)
    # get the bounds and time points from two extremes
    bounds = [fcstAxFirst.bounds[0][0], fcstAxLast.bounds[-1][-1]]
    
    if action is 'sum' and bounds[0] != 0:
        # this change is required only for _accumulationVars_ vars, since its 
        # hourly accumulation, which we converting to 6-hourly accumulation.
        # Instead of cross check by _accumulationVars_, here we are checking
        # by action is 'sum', since sum arg passed only to _accumulationVars_.
        bounds = [fcstAxFirst.bounds[0][1], fcstAxLast.bounds[-1][-1]]
    # end of if ...:    
    
    if fpoint == 'cbound':
        #### THE CENTRE POINT OF FORECAST TIME BOUNDS
        fcstpoint = [bounds[0] + ((bounds[-1] - bounds[0]) / 2.0)]
    elif fpoint == 'lbound':           
        ###### THE START BOUNDS OF FORECAST TIME POINT
        fcstpoint = [bounds[0]]         
    elif fpoint == 'rbound':           
        ###### THE END BOUNDS OF FORECAST TIME POINT
        fcstpoint = [bounds[-1]] 
    # end of if fpoint == 'cbound':   
    
    # update the time coordinate with new fcst time point and fcst time bounds
    fcstAxFirst.points = fcstpoint
    fcstAxFirst.bounds = bounds
    # add the updated fcst time coordinate to the meanCube
    meanCube.add_aux_coord(fcstAxFirst)
    # add attributes back to meanCube
    meanCube.attributes = tmpCube.attributes  
    # add standard_name
    meanCube.standard_name = tmpCube.standard_name
    meanCube.long_name = tmpCube.long_name if tmpCube.long_name else tmpCube.standard_name
    
    if action == 'mean':
        cm = iris.coords.CellMethod('mean', ('time',), intervals=(dt,), 
                                     comments=(actionIntervals+' mean',))
    elif action == 'sum':
        cm = iris.coords.CellMethod('sum', ('time',), intervals=(dt,), 
                                     comments=(actionIntervals+' accumulation',))
                                     
    # add cell_methods to the meanCube                                     
    meanCube.cell_methods = (cm,) 

    # make memory free
    del tmpCube
    
    # return mean cube 
    return meanCube
# end of def cubeAverager(...):


def cubeSubtractor(cube, otherCube, standard_name=None, 
                        long_name=None, removeSTASH=True):
    '''
    iris.analysis.maths.subtract doesnt serving out purpose.
    So this cubeSubtractor subtract two cubes and set all necessary meta 
    information.

    Arulalan.T
    12-Fen-2016
    '''
    
    timeDimension = True if len(cube.shape) == 3 else False
    attr = cube.attributes
    if removeSTASH and 'STASH' in attr: attr.pop('STASH')
    lname = long_name if long_name else cube.long_name
    sname = standard_name if standard_name else cube.standard_name
    cm = cube.cell_methods[0]  
    unit = cube.units 
    # do the simple substraction
    subtracted = cube.data - otherCube.data
    # just set proper masked array
    subtracted = numpy.ma.masked_less(subtracted, 1e-15)
    subtracted  = numpy.ma.masked_greater(subtracted, 1e+15)                              
    numpy.ma.set_fill_value(subtracted, 9.999e+20) 
    
    subtracted = iris.cube.Cube(data=subtracted, units=unit, standard_name=sname, 
                       long_name=lname, attributes=attr, cell_methods=(cm,))
    idx = 0
    if timeDimension:
        subtracted.add_dim_coord(cube.coords('time')[0], idx)
        idx += 1        
    subtracted.add_dim_coord(cube.coords('latitude')[0], idx)
    subtracted.add_dim_coord(cube.coords('longitude')[0], idx+1)
    for axc in cube.aux_coords: subtracted.add_aux_coord(axc)    
    # return the updated, subtracted cube
    return subtracted
# end of def cubeSubtractor(...):
