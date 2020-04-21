"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          local_polar
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Performs local-polar and polar-local transformation of origin
                and destination points.
-------------------------------------------------------------------------------
"""


import math
from c3m_3d.validate import missing


class local_polar():


    def __init__(self, **kwargs):

        """
        Parameters
        -----------------------------------------------------------------------
        xL_0        Local plane x-value of the origin point
        xL          Local plane x-value of the destination point
        yL_0        Local plane y-value of the origin point
        yL          Local plane y-value of the destination point
        zL_0        Local plane z-value of the origin point
        zL          Local plane z-value of the destination point
        alphaL_0    Local plane forward azimuth from the origin point (in
                    radians)
        zetaL_0     Local plane zenith angle at the origin point (in radians)
        sL          Local plane distance between origin and desintation points
        -----------------------------------------------------------------------

        Calculated Values
        -----------------------------------------------------------------------
        xL          Local plane x-value of the destination point (if not in
                    parameters)
        yL          Local plane y-value of the destination point (if not in
                    parameters)
        zL          Local plane z-value of the destination point (if not in
                    parameters)
        alphaL_0    Local plane forward azimuth from the origin point (in
                    radians; if not in parameters)
        alphaL      Local plane back azimuth to the origin point (in radians;
                    if not in parameters)
        zetaL_0     Local plane zenith angle at the origin point (in radians;
                    if not in parameters)
        sL          Local plane distance between origin and desintation points
                    (if not in parameters)
        -----------------------------------------------------------------------

        Usage
        -----------------------------------------------------------------------
        Instantiation of a local_polar object will perform calculations based
        on provided parameters.

        Use of the local_polar.data() method provides parameters and calulcated
        values.
        
        Manual calculation is also available through object attribute
        assignment and use of the local_polar.build() method.
        -----------------------------------------------------------------------
        """

        self.params = kwargs
        self.__build()


    def __build(self):
        
        self.params.update(**self.__local_to_polar(**self.params))
        self.params.update(**self.__polar_to_local(**self.params))
        return


    def __polar_to_local(self, xL=None, yL=None, zL=None, sL=None,
                         alphaL_0=None, zetaL_0=None, **kwargs):
        # Performs polar-local transform for local destination point
        if missing.all(xL, yL, zL) and not missing.any(sL, alphaL_0, zetaL_0):
            xL, yL, zL = self.__local_destination(sL, alpha_0, zeta_0)
            kwargs['xL'], kwargs['yL'], kwargs['zL'] = xL, yL, zL
        return kwargs


    def __local_to_polar(self, sL=None, alphaL_0=None, alphaL=None,
                         zetaL_0=None, xL_0=None, yL_0=None, zL_0=None,
                         xL=None, yL=None, zL=None, **kwargs):
        # Performs local-polar transform for vector between local origin and
        # destination points
        if missing.all(sL, alphaL_0, zetaL_0) and not missing.any(
            xL_0, yL_0, zL_0, xL, yL, zL):
            alphaL_0, zetaL_0, sL = self.__polar_vector(
                xL_0, yL_0, zL_0, xL, yL, zL)
            kwargs['alphaL_0'], kwargs['zetaL_0'], kwargs['sL'] = (
                alphaL_0, zetaL_0, sL)
        if alphaL_0 and not alphaL:
            kwargs['alphaL'] = self.__back_azimuth(alphaL_0)
        return kwargs


    def __local_destination(self, sL, alphaL_0, zetaL_0):
        # Local destination point in three dimensions based on vector
        sin_alphaL_0 = math.sin(alphaL_0)
        cos_alphaL_0 = math.cos(alphaL_0)
        sin_zetaL_0 = math.sin(zetaL_0)
        cos_zetaL_0 = math.cos(zetaL_0)
        xL = sL * cos_alphaL_0 * sin_zetaL_0
        yL = sL * sin_alphaL_0 * sin_zetaL_0
        zL = sL * cos_zetaL_0
        return xL, yL, zL


    def __polar_vector(self, xL_0, yL_0, zL_0, xL, yL, zL):
        sL = self.__ground_distance(xL_0, yL_0, zL_0, xL, yL, zL)
        alphaL_0 = self.__forward_azimuth(xL_0, yL_0, xL, yL)
        zetaL_0 = self.__zenith_angle(zL_0, zL, sL)
        return alphaL_0, zetaL_0, sL


    def __forward_azimuth(self, xL_0, yL_0, xL, yL):
        # Azimuth from the origin point, clockwise from north
        delta_xL, delta_yL = xL - xL_0, yL - yL_0
        alphaL_0 = math.atan2(delta_xL, delta_yL)
        if alphaL_0 < 0: # Ensures azimuth is positive
            alphaL_0 += math.pi * 2
        return alphaL_0


    def __back_azimuth(self, alphaL_0):
        # Azimuth from the destination point, clockwise from north
        if alphaL_0 <= math.pi:
            alphaL = alphaL_0 + math.pi
        else:
            alphaL = alphaL_0 - math.pi
        return alphaL


    def __ground_distance(self, xL_0, yL_0, zL_0, xL, yL, zL):
        # Three-dimensional distance from origin to destination points
        delta_xL, delta_yL, delta_zL = xL - xL_0, yL - yL_0, zL - zL_0
        delta_xL2, delta_yL2, delta_zL2 = delta_xL**2, delta_yL**2, delta_zL**2
        sum_delta2 = delta_xL2 + delta_yL2 + delta_zL2
        return math.sqrt(sum_delta2)


    def __zenith_angle(self, zL_0, zL, sL):
        # Angle between zenith of origin point and vector to destination point
        delta_zL = zL - zL_0
        delta_zL_sL = delta_zL / sL
        return math.acos(delta_zL_sL)


    def data(self):

        return self.params


    def build(self):

        self.__build()
        return
