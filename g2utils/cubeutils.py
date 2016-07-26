import numpy, iris  

def cubeAverager(tmpCube, action='mean', dt='1 hour', 
                actionIntervals='6 hour', tpoint='cbound', fpoint='cbound',
                tbounds=True, fbounds=True):
    """
    :param tmpCube:     The temporary cube data (in Iris format) with non-singleton time dimension
    :param action:      mean| sum (accumulated fields are summed and instantaneous are averaged).
    :param dt:   A standard string representing forecast step duration/intervals.
    :param actionIntervals: A non standard string to add inside cell_methods comments section.
    :param tpoint: cbound | lbound | rbound time point 
    :param fpoint: cbound | lbound | rbound forecast period
    :param tbounds: True | False (False removes time bounds)
    :param fbounds: True | False (False removes fcst bounds)
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

    # get the reference time bounds and time points from two extremes    
    rbounds = [timeAxFirst.bounds[0][0], timeAxLast.bounds[-1][-1]]
    
    if tpoint == 'cbound':
        #### THE CENTRE POINT OF REFERENCE TIME BOUNDS
        timepoint = [rbounds[-1] + ((rbounds[-1] - rbounds[0]) / 2.0)]
    elif tpoint == 'lbound':
       ###### THE START BOUNDS OF REFERENCE TIME POINT
        timepoint = [rbounds[0]]
    elif tpoint == 'rbound':    
        ###### THE END BOUNDS OF REFERENCE TIME POINT
        timepoint = [rbounds[-1]]
    # end of if tpoint == 'cbound':
    
    # update the time coordinate with new time point and time bounds
    timeAxFirst.points = timepoint
    # set time bounds if enabled
    if tbounds: 
        timeAxFirst.bounds = rbounds
    else:
        timeAxFirst.bounds = None
    # add the updated time coordinate to the meanCube
    meanCube.add_aux_coord(timeAxFirst)
    
    # get the forecat time bounds and time points from two extremes 
    fbounds = [fcstAxFirst.bounds[0][0], fcstAxLast.bounds[-1][-1]]

    if action is 'sum' and fbounds[0] != 0:
        # this change is required only for _accumulationVars_ vars, since its 
        # hourly accumulation, which we converting to 6-hourly accumulation.
        # Instead of cross check by _accumulationVars_, here we are checking
        # by action is 'sum', since sum arg passed only to _accumulationVars_.
        fbounds = [fcstAxFirst.bounds[0][1], fcstAxLast.bounds[-1][-1]]
    # end of if ...:    
    
    if fpoint == 'cbound':
        #### THE CENTRE POINT OF FORECAST TIME BOUNDS
        fcstpoint = [fbounds[0] + ((fbounds[-1] - fbounds[0]) / 2.0)]
    elif fpoint == 'lbound':           
        ###### THE START BOUNDS OF FORECAST TIME POINT
        fcstpoint = [fbounds[0]]         
    elif fpoint == 'rbound':           
        ###### THE END BOUNDS OF FORECAST TIME POINT
        fcstpoint = [fbounds[-1]]        
    # end of if fpoint == 'cbound':   
    
    # update the time coordinate with new fcst time point and fcst time bounds
    fcstAxFirst.points = fcstpoint
    # set fcst bounds if enabled
    if fbounds: 
        fcstAxFirst.bounds = fbounds
    else:
        fcstAxFirst.bounds = None
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
    if fbounds: meanCube.cell_methods = (cm,)   
    
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
    cm = cube.cell_methods[0] if cube.cell_methods else None 
    unit = cube.units 
    # do the simple substraction
    subtracted = cube.data - otherCube.data
    # set fill_value 
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

def cubeRealizationAverager(tmpCube, action='mean', dr='1 ENS'):
    """
    :param tmpCube:     The temporary cube data (in Iris format) with non-singleton time dimension
    :param action:      mean| sum (accumulated fields are summed and instantaneous are averaged).
    :param dr:   A standard string representing realization step duration/intervals.
    
    :return: meanCube:  An Iris formatted cube date containing the resultant data either as
                        averaged or summed.
    Arulalan.T
    18-July-2016
    """

    # extract realization points 
    rpoints = tmpCube.coord('realization').points
    # get the time coord of first time cube and set to mean
    timeAxFirst = tmpCube.coords('time')[0]
    fcstAxFirst = tmpCube.coords('forecast_period')[0]
    # assign first realization point data to mean data 
    meanCube = tmpCube.extract(iris.Constraint(realization=rpoints[0]))
       
    # loop through remaining realization points 
    for rp in rpoints[1:]:
        # extract remaining realization points and add to meanCube
        lastCube = tmpCube.extract(iris.Constraint(realization=rp))
        meanCube = iris.analysis.maths.add(meanCube, lastCube) 
    # end of for rp in rpoints[1:]:
           
    if action == 'mean':
        # to compute mean value of meanCube divide it by length of time
        # points which we added before 
        meanCube = iris.analysis.maths.divide(meanCube, float(len(rpoints)))
        print "Converting cube to realization mean : %s" % (tmpCube.standard_name)
    elif action == 'sum':
        # for sum action, we no need to do here anything, because already 
        # we computed accumulation only!
        print "Converting cube to realization accumulation : %s" % (tmpCube.standard_name)
    else:
        raise ValueError('argument "%s" not support' % action)
    # end of if not isAccumulation:    
       
    if action == 'mean':
        cm = iris.coords.CellMethod('mean', ('realization',), intervals=(dr,), 
                                     comments=('mean of (1 control run + %d ensemble members)' % (len(rpoints)-1)))
    elif action == 'sum':
        cm = iris.coords.CellMethod('sum', ('realization',), intervals=(dr,), 
                                     comments=('accumulation of (1 control run + %d ensemble members)' % (len(rpoints)-1)))
    
    meanCube = iris.cube.Cube(data=meanCube.data, units=tmpCube.units, 
                       standard_name=tmpCube.standard_name, 
                       long_name=tmpCube.long_name, 
                       attributes=tmpCube.attributes,
                       cell_methods=(cm,))
        
    idx = 0
    for ax in tmpCube.coords(): 
        name = ax.standard_name if ax.standard_name else ax.long_name
        if name == 'realization': continue
        if isinstance(ax, iris.coords.DimCoord):      
            print ax
            if ax.standard_name in ['forecast_period', 'time'] and len(ax.points) == 1:
                meanCube.add_aux_coord(ax)
            else:
                meanCube.add_dim_coord(ax, idx)
                idx += 1
        elif isinstance(ax, iris.coords.AuxCoord):
            meanCube.add_aux_coord(ax)
    # end of for ax in tmpCube.coords(): 
    print meanCube    
    
    # make memory free
    del tmpCube
    
    # return mean cube 
    return meanCube
# end of def cubeRealizationAverager(...):
