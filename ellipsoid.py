"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          ellipsoid
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Calculates the ellipsoid values for a given datum for use in
                cadastral measurement calculations
-------------------------------------------------------------------------------
"""


import math
from c3m_3d.validate import missing


class ellipsoid():


    def __init__(self, **kwargs):

        """
        Parameters
        -----------------------------------------------------------------------
        a           Semi-major axis of the ellipsoid
        b           Semi-minor axis of the ellipsoid
        f           Flattening of the ellipsoid
        -----------------------------------------------------------------------

        Calculated Values
        -----------------------------------------------------------------------
        a           Semi-major axis of the ellipsoid (if not in parameters)
        b           Semi-minor axis of the ellipsoid (if not in parameters)
        f           Flattening of the ellipsoid (if not in parameters)
        epsilon     Linear eccentiricity of the ellipsoid
        e           First eccenstricity of the ellipsoid
        e_prime     Second eccentricity of the ellipsoid
        -----------------------------------------------------------------------

        Usage
        -----------------------------------------------------------------------
        Instantiation of an ellpsoid object will calculate eccentricity values
        based on provided parameters.

        Use of the curvature.data() method provides parameters and calculated
        values.
        
        Manual calculation is also available through object attribute
        assignment and use of the ellipsoid.build() method.
        -----------------------------------------------------------------------
        """
        
        self.params = kwargs
        self.__build()
        

    def __build(self):

        self.params.update(**self.__ellipsoid(**self.params))
        return


    def __ellipsoid(self, a=None, b=None, f=None, **kwargs):
        
        if a and b and not f:
            f = (a - b) / a
            kwargs['f'] = f
        elif a and f and not b:
            b = a * (1 - f)
            kwargs['b'] = b
        elif b and f and not a:
            a = b / (1 - f)
            kwargs['a'] = a
        if a and b:
            a2, b2 = a**2, b**2
            epsilon = math.sqrt(a2 - b2)
            e, e_prime = epsilon / a, epsilon / b
            kwargs['epsilon'], kwargs['e'], kwargs['e_prime'], = (
                epsilon, e, e_prime)
        return kwargs
        

    def data(self):
        
        return self.params


    def build(self):

        self.__build()
        
        return


    
