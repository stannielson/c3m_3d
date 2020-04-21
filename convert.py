"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          convert
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Series of functions to perform conversion of values for
                provided sequences or values.
-------------------------------------------------------------------------------
"""


import math


def all_float(args):
    
    # Converts list/tuple/set numbers to float and returns a tuple

    return tuple(float(i) if i else i for i in args)


def all_radians(args):

    # Converts list/tuple/set numbers to radians and returns a tuple

    return tuple(math.radians(i) if i else i for i in args)


def all_degrees(args):

    # Converts list/tuple/set numbers to degrees and returns a tuple

    return tuple(math.degrees(i) if i else i for i in args)


def all_float2(kwargs):

    # Converts dictionary value numbers to float and returns a dictionary

    return dict((i, float(j)) if j else (i, j) for i, j in kwargs.items())


def all_radians2(kwargs):

    # Converts dictionary value numbers to radians and returns a dictionary

    return dict((i, math.radians(j)) if j else (i, j) for i, j in
                kwargs.items())


def all_degrees2(kwargs):

    # Converts dictionary value numbers to degrees and returns a dictionary

    return dict((i, math.degrees(j)) if j else (i, j) for i, j in
                kwargs.items())


def m_to_ft(meters):
    # Converts meters to feet
    return meters * 3750.0 / 1143.0


def ft_to_m(feet):

    # Converts feet to meters

    return feet * 1143.0 / 3750.0


def coord_dms(value, ang_meas):

    # Converts coordinate values from decimal degrees to degrees-minutes-
    # seconds for display
    
    conds = [(i, j) for i in ['+', '-'] for j in ['lat', 'lon']]
    dirs = ['N', 'E', 'S', 'W']
    _ = list(zip(conds, dirs))
    crit = {i: j for i, j in _}
    pos = abs(value)
    if value >= 0:
        oper = '+'
    else:
        oper = '-'
    deg = int(pos)
    min_ = int(60 * (pos - int(pos)))
    sec = ((60 * (pos - int(pos))) - int(60 * (pos - int(pos)))) * 60
    dir_ = crit.get((oper, ang_meas))
    
    return u'{}\u00B0{}\'{}"{}'.format(deg, min_, sec, dir_)


def dd_dms(value, ang_meas):

    # Converts decimal degrees to degrees-minutes-seconds for display
    
    pos = abs(value)
    deg = int(value)
    min_ = int(60 * (pos - int(pos)))
    sec = ((60 * (pos - int(pos))) - int(60 * (pos - int(pos)))) * 60
    dir_ = crit.get((oper, ang_meas))
    return u'{}\u00B0{}\'{}"'.format(deg, min_, sec)
