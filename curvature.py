"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          curvature
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Calculates the curvature values for origin and destination
                points on a ellipsoid for use in cadastral measurement
                calculations.
-------------------------------------------------------------------------------
"""


import math
from c3m_3d.validate import missing


class curvature():


    def __init__(self, **kwargs):

        """
        Parameters
        -----------------------------------------------------------------------
        phi_0       Latitude of the origin point (in radians)
        phi         Latitude of the destination point (in radians)
        a           Semi-major axis of the ellipsoid
        e           First eccentricity of the ellispoid
        -----------------------------------------------------------------------

        Calculated Values
        -----------------------------------------------------------------------
        M_0         Meridian radius of curvature for the origin point
        M           Meridian radius of curvature for the destination point
        N_0         Radius of curvature for the origin point
        N           Radius of curvature for the destination point
        -----------------------------------------------------------------------

        Usage
        -----------------------------------------------------------------------
        Instantiation of a curvature object will calculate origin and
        destination curvature values based on provided parameters.

        Use of the curvature.data() method provides parameters and calculated
        values.

        Manual calculation is also available through object attribute
        assignment and use of the curvature.build() method.
        -----------------------------------------------------------------------
        """
        
        self.params = kwargs
        self.__build()


    def __build(self):

        self.params.update(**self.__curvature(**self.params))
        
        return


    def __curvature(self, phi_0=None, M_0=None, N_0=None, phi=None,
                    M=None, N=None, a=None, e=None, **kwargs):
        
        # Calculates curvature values for origin and destination coordinate
        # points based on geodetic latitude, semi-major axis, and first
        # eccentricity of the ellipsoid
        
        if not missing.any(a, e, phi_0):
            if not M_0:
                kwargs['M_0'] = self.__M(a, e, phi_0)
            if not N_0:
                kwargs['N_0'] = self.__N(a, e, phi_0)
        if not missing.any(a, e, phi):
            if not M:
                kwargs['M'] = self.__M(a, e, phi)
            if not N:
                kwargs['N'] = self.__N(a, e, phi)
                
        return kwargs


    def __M(self, a, e, phi):
        
        # Calculates the meridian radius of curvature for a point based on
        # geodetic latitude, along with the semi-major axis and first
        # eccentricity of the ellipsoid

        e2, sin_phi = e**2, math.sin(phi)
        sin2_phi = sin_phi**2
        M_num = a * (1 - e2)
        M_denom = (1 - e2 * sin2_phi)**(3.0/2.0)
        
        return M_num / M_denom


    def __N(self, a, e, phi):
        
        # Calculates the radius of curvature for a point based on geodetic
        # latitude, along with the semi-major axis and first eccentricity of the
        # ellipsoid
        
        e2, sin_phi = e**2, math.sin(phi)
        sin2_phi = sin_phi**2
        _1_e2_sin2_phi = 1 - e2 * sin2_phi
        
        return a / math.sqrt(_1_e2_sin2_phi)


    def data(self):

        return self.params


    def build(self):

        self.__build()
        
        return



