"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          cartesian_local
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Performs cartesian-local and local-cartesian transformation of
                origin and destination points.
-------------------------------------------------------------------------------
"""


import math
from c3m_3d.validate import missing


class cartesian_local():


    def __init__(self, **kwargs):

        """
        Parameters
        -----------------------------------------------------------------------
        phi_0       Latitude of the origin point (in radians)
        phi         Latitude of the destination point (in radians)
        lambda__0   Longitude of the origin point (in radians)
        lambda_     Longitude of the origin point (in radians)
        x_0         Cartesian x-value of the origin point
        x           Cartesian x-value of the destination point
        y_0         Cartesian y-value of the origin point
        y           Cartesian y-value of the destination point
        z_0         Cartesian z-value of the origin point
        z           Cartesian z-value of the destination point
        xL_0        Local plane x-value of the origin point
        xL          Local plane x-value of the destination point
        yL_0        Local plane y-value of the origin point
        yL          Local plane y-value of the destination point
        zL_0        Local plane z-value of the origin point
        zL          Local plane z-value of the destination point        
        -----------------------------------------------------------------------

        Calculated Values
        -----------------------------------------------------------------------
        x_0         Cartesian x-value of the origin point (if not in
                    parameters)
        x           Cartesian x-value of the destination point (if not in
                    parameters)
        y_0         Cartesian y-value of the origin point (if not in
                    parameters)
        y           Cartesian y-value of the destination point (if not in
                    parameters)
        z_0         Cartesian z-value of the origin point (if not in
                    parameters)
        z           Cartesian z-value of the destination point (if not in
                    parameters)
        xL_0        Local plane x-value of the origin point (always 0.0)
        xL          Local plane x-value of the destination point (if not in
                    parameters)
        yL_0        Local plane y-value of the origin point (always 0.0)
        yL          Local plane y-value of the destination point (if not in
                    parameters)
        zL_0        Local plane z-value of the origin point (always 0.0)
        zL          Local plane z-value of the destination point (if not in
                    parameters)
        -----------------------------------------------------------------------

        Usage
        -----------------------------------------------------------------------
        Instantiation of a cartesian_local object will perform calculations
        based on provided parameters.

        Use of the cartesian_local.data() method provides parameters and
        calculated values.
        
        Manual calculation is also available through object attribute
        assignment and use of the cartesian_local.build() method.
        -----------------------------------------------------------------------
        """

        self.params = kwargs
        self.__build()


    def __build(self):

        self.params.update(**self.__cartesian_to_local(**self.params))
        self.params.update(**self.__local_to_cartesian(**self.params))
        
        return


    def __cartesian_to_local(self, x_0=None, y_0=None, z_0=None, x=None,
                             y=None, z=None, xL_0=None, yL_0=None, zL_0=None,
                             xL=None, yL=None, zL=None, phi_0=None,
                             lambda__0=None, **kwargs):

        # Performs Cartesian-local transform for destination coordinate point

        if missing.all(xL_0, yL_0, zL_0):
            xL_0, yL_0, zL_0 = 0.0, 0.0, 0.0
            kwargs['xL_0'], kwargs['yL_0'], kwargs['zL_0'] = xL_0, yL_0, zL_0
        if missing.all(xL, yL, zL) and not missing.any(
            x_0, y_0, z_0, x, y, z, xL_0, yL_0, zL_0, phi_0, lambda__0):
            xL, yL, zL = self.__local_destination(
                x_0, y_0, z_0, x, y, z, xL_0, yL_0, zL_0, phi_0, lambda__0)
            kwargs['xL'], kwargs['yL'], kwargs['zL'] = xL, yL, zL

        return kwargs


    def __local_to_cartesian(self, x_0=None, y_0=None, z_0=None, x=None,
                             y=None, z=None, xL_0=None, yL_0=None, zL_0=None,
                             xL=None, yL=None, zL=None, phi_0=None,
                             lambda__0=None, **kwargs):
        
        # Performs Cartesian-local transform for destination coordinate point

        if missing.all(xL_0, yL_0, zL_0):
            xL_0, yL_0, zL_0 = 0.0, 0.0, 0.0
            kwargs['xL_0'], kwargs['yL_0'], kwargs['zL_0'] = xL_0, yL_0, zL_0
        if missing.all(x, y, z) and not missing.any(
            xL_0, yL_0, zL_0, xL, yL, zL, x_0, y_0, z_0, phi_0, lambda__0):
            x, y, z = self.__cartesian_destination(
                xL_0, yL_0, zL_0, xL, yL, zL, x_0, y_0, z_0, phi_0, lambda__0)
            kwargs['x'], kwargs['y'], kwargs['z'] = x, y, z

        return kwargs


    def __cartesian_destination(self, xL_0, yL_0, zL_0, xL, yL, zL, x_0, y_0,
                                z_0, phi_0, lambda__0):
        
        # Transforms a tangent-/secant-plane destination into a Cartesian
        # point based on geometric distance difference and geodetic origin
        # between tangent-/secant-plane origin and destination coordinates

        delta_xL, delta_yL, delta_zL = xL - xL_0, yL - yL_0, zL - zL_0
        sin_phi_0, cos_phi_0 = math.sin(phi_0), math.cos(phi_0)
        sin_lambda__0, cos_lambda__0 = math.sin(lambda__0), math.cos(lambda__0)
        x = (delta_xL * -sin_lambda__0 + 
             delta_yL * -cos_lambda__0 * sin_phi_0 +
             delta_zL * cos_phi_0 * cos_lambda__0 + 
             x_0)
        y = (delta_xL * cos_lambda__0 +
             delta_yL * -sin_phi_0 * sin_lambda__0 +
             delta_zL * cos_phi_0 * sin_lambda__0 +
             y_0)
        z = (delta_yL * cos_phi_0 +
             delta_zL * sin_phi_0 +
             z_0)

        return x, y, z


    def __local_destination(self, x_0, y_0, z_0, x, y, z, xL_0, yL_0, zL_0,
                            phi_0, lambda__0):
        
        # Transforms a Cartesian destination into a tangent-/secant-plane
        # point based on geometric distance difference and geodetic origin
        # between Cartesian origin and destination coordinates

        delta_x, delta_y, delta_z = x - x_0, y - y_0, z - z_0
        sin_phi_0, cos_phi_0 = math.sin(phi_0), math.cos(phi_0)
        sin_lambda__0, cos_lambda__0 = math.sin(lambda__0), math.cos(lambda__0)
        xL = (delta_x * -sin_lambda__0 +
              delta_y * cos_lambda__0 +
              xL_0)
        yL = (delta_x * -cos_lambda__0 * sin_phi_0 +
              delta_y * -sin_phi_0 * sin_lambda__0 +
              delta_z * cos_phi_0 +
              yL_0)
        zL = (delta_x * cos_phi_0 * cos_lambda__0 +
              delta_y * cos_phi_0 * sin_lambda__0 +
              delta_z * sin_phi_0 +
              zL_0)

        return xL, yL, zL


    def data(self):

        return self.params


    def build(self):

        self.__build()
        
        return
