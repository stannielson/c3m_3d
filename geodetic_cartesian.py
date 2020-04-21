"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          geodetic_cartesian
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Performs geodetic-cartesian and cartesian-geodetic
                transformation of origin and destination points.
-------------------------------------------------------------------------------
"""


import math
from c3m_3d.validate import missing


class geodetic_cartesian():


    def __init__(self, **kwargs):

        """
        Parameters
        -----------------------------------------------------------------------
        phi_0       Latitude of the origin point (in radians)
        phi         Latitude of the destination point (in radians)
        lambda__0   Longitude of the origin point (in radians)
        lambda_     Longitude of the origin point (in radians)
        h_0         Height above the ellipsoid of the origin point
        h           Height above the ellispoid of the destination point
        x_0         Cartesian x-value of the origin point
        x           Cartesian x-value of the destination point
        y_0         Cartesian y-value of the origin point
        y           Cartesian y-value of the destination point
        z_0         Cartesian z-value of the origin point
        z           Cartesian z-value of the destination point
        N_0         Radius of curvature for the origin point
        N           Radius of curvature for the destination point
        a           Semi-major axis of the ellipsoid
        b           Semi-minor axis of the ellipsoid
        e           First eccenstricity of the ellipsoid
        -----------------------------------------------------------------------

        Calculated Values
        -----------------------------------------------------------------------
        phi_0       Latitude of the origin point (in radians; if not in
                    parameters)
        phi         Latitude of the destination point (in radians; if not in
                    parameters)
        lambda__0   Longitude of the origin point (in radians; if not in
                    parameters)
        lambda_     Longitude of the origin point (in radians; if not in
                    parameters)
        h_0         Height above the ellipsoid of the origin point (if not in
                    parameters)
        h           Height above the ellispoid of the destination point (if not
                    in parameters)
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
        -----------------------------------------------------------------------

        Usage
        -----------------------------------------------------------------------
        Instantiation of a geodetic_cartesian object will perform calculations
        based on provided parameters.

        Use of the geodetic_cartesian.data() method provides parameters and
        calculated values.
        
        Manual calculation is also available through object attribute
        assignment and use of the geodetic_cartesian.build() method.
        -----------------------------------------------------------------------
        """

        self.params = kwargs
        self.__build()


    def __build(self):
        
        self.params.update(**self.__geodetic_to_cartesian(**self.params))
        self.params.update(**self.__cartesian_to_geodetic(**self.params))
        
        return


    def __geodetic_to_cartesian(self, phi_0=None, lambda__0=None, h_0=None,
                                N_0=None, phi=None, lambda_=None, h=None,
                                N=None, a=None, b=None, x_0=None, y_0=None,
                                z_0=None, x=None, y=None, z=None, **kwargs):
        
        # Performs geodetic-cartesian transforms for both origin and
        # destination coordinate points
        
        if missing.all(x_0, y_0, z_0) and not missing.any(
            phi_0, lambda__0, h_0, N_0, a, b):
            kwargs['x_0'], kwargs['y_0'], kwargs['z_0'] = (
                self.__cartesian_transform(phi_0, lambda__0, h_0, N_0, a, b))
        if missing.all(x, y, z) and not missing.any(
            phi, lambda_, h, N, a, b):
            kwargs['x'], kwargs['y'], kwargs['z'] = (
                self.__cartesian_transform(phi, lambda_, h, N, a, b))

        return kwargs


    def __cartesian_to_geodetic(self, x_0=None, y_0=None, z_0=None, x=None,
                                y=None, z=None,  N=None, a=None, b=None,
                                e=None, phi_0=None, lambda__0=None, h_0=None,
                                N_0=None, phi=None, lambda_=None, h=None,
                                **kwargs):
        
        # Performs cartesian-geodetic transforms for both origin and
        # destination coordinate points

        if missing.all(phi_0, lambda__0, h_0, N_0) and not missing.any(
            x_0, y_0, z_0, a, b, e):
            (kwargs['phi_0'], kwargs['lambda__0'], kwargs['h_0'],
             kwargs['N_0']) = (
                 self.__geodetic_transform(x_0, y_0, z_0, a, b, e))
        if missing.all(phi, lambda_, h, N) and not missing.any(
            x, y, z, a, b, e):
            kwargs['phi'], kwargs['lambda_'], kwargs['h'], kwargs['N'] = (
                self.__geodetic_transform(x, y, z, a, b, e))

        return kwargs


    def __geodetic_transform(self, x, y, z, a, b, e):
        
        # Transforms Cartesian (Earth-centered, Earth-fixed [ECEF]) coordinates
        # into geodetic coordinates

        phi, h = self.__phi_h(x, y, z, a, b)
        lambda_ = self.__lambda(x, y)
        N = self.__N(a, e, phi)

        return phi, lambda_, h, N



    def __phi_h(self, x, y, z, a, b):
        
        # Calculates geodetic latitude and height based on a cartesian point,
        # semi-major and semi-minor axes. Calculation uses Newtonian iteration
        # to produce latitude and height using the method developed by Lin and
        # Wang (1995)

        x2, y2, z2, a2, b2, a4, b4 = x**2, y**2, z**2, a**2, b**2, a**4, b**4
        p2 = x2 + y2
        p = math.sqrt(p2)
        ab, a2b2 = a * b, a2 * b2
        a2z2, b2p2, a4z2, b4p2 = a2 * z2, b2 * p2, a4 * z2, b4 * p2
        a2z2_b2p2, a4z2_b4p2 = a2z2 + b2p2, a4z2 + b4p2
        t1, t2, t3 = ab * a2z2_b2p2**1.5, a2b2 * a2z2_b2p2, 2 * a4z2_b4p2
        m = (t1 - t2) / t3
        counter = 0
        while True:
            a_2ma, b_2mb = a + 2 * m / a, b + 2 * m / b
            a_2ma2, b_2mb2 = a_2ma**2, b_2mb**2
            a_a_2ma3, b_b_2mb3 = a * a_2ma**3, b * b_2mb**3
            f_m = p2 / a_2ma2 + z2 / b_2mb2 - 1
            f0_m = -4 * (p2 / a_a_2ma3 + z2 / b_b_2mb3)
            m_n = m - f_m / f0_m
            counter += 1
            if m == m_n or counter == 1000:
                m = m_n
                break
            m = m_n
        _1_2ma2, _1_2mb2 = 1 + 2 * m / a**2, 1 + 2 * m / b**2
        p_e, z_e = p / _1_2ma2, z / _1_2mb2
        delta_p2, delta_z2 = (p - p_e)**2, (z - z_e)**2
        h = math.sqrt(delta_p2 + delta_z2)
        if p_e + abs(z) < p_e + abs(z_e):
            h *= -1
        a2z_e, b2p_e = a2 * z_e, b2 * p_e

        return math.atan2(a2z_e, b2p_e), h


    def __lambda(self, x, y):
        
        # Calculates geodetic longitude based on cartesian (x,y) values

        return math.atan2(y, x)


    def __h(self, x, y, phi, N):

        # Calculates geodetic height based on cartesian (x,y) values and
        # geodetic latitude

        x2, y2 = x**2, y**2
        cos_phi = math.cos(phi)

        return (math.sqrt(x2 + y2) / cos_phi) - N


    def __N(self, a, e, phi):
        
        # Calculates the radius of curvature for a point based on geodetic
        # latitude, along with the semi-major axis and first eccentricity of the
        # ellipsoid

        e2, sin_phi = e**2, math.sin(phi)
        sin2_phi = sin_phi**2
        _1_e2_sin2_phi = 1 - e2 * sin2_phi

        return a / math.sqrt(_1_e2_sin2_phi)
        


    def __cartesian_transform(self, phi, lambda_, h, N, a, b):

        # Transforms geodetic coordinates into Cartesian (Earth-centered,
        # Earth-fixed [ECEF]) coordinates

        sin_phi, cos_phi = math.sin(phi), math.cos(phi)
        sin_lambda_, cos_lambda_ = math.sin(lambda_), math.cos(lambda_)
        Nh, a2, b2 = N + h, a**2, b**2
        b2_a2 = b2 / a2
        N_b2_a2 = N * b2_a2
        N_b2_a2_h = N_b2_a2 + h
        x = Nh * cos_phi * cos_lambda_
        y = Nh * cos_phi * sin_lambda_
        z = N_b2_a2_h * sin_phi

        return x, y, z


    def data(self):
        
        return self.params


    def build(self):
        
        self.__build()

        return
