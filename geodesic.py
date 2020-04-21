"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          geodesic
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Performs Vincenty direct and inverse transformation of origin
                and destination points based on geodesic and geodetic
                measurements.
-------------------------------------------------------------------------------
"""


import math
from c3m_3d.validate import missing


class geodesic(object):


    def __init__(self, **kwargs):

        """
        Parameters
        -----------------------------------------------------------------------
        phi_0       Latitude of the origin point (in radians)
        phi         Latitude of the destination point (in radians)
        lambda__0   Longitude of the origin point (in radians)
        lambda_     Longitude of the origin point (in radians)
        s           Geodesic on the ellipsoid from origin to destination points
        alpha_0     Geodetic forward azimuth from the origin point
        alpha       Geodetic back azimuth to the origin point
        alpha_g     Geodetic mean azimuth
        a           Semi-major axis of the ellipsoid
        b           Semi-minor axis of the ellipsoid
        f           Flattening of the ellipsoid
        
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
        s           Geodesic on the ellipsoid from origin to destination points
                    (if not in parameters)
        alpha_0     Geodetic forward azimuth from the origin point (in radians;
                    if not in parameters)
        alpha       Geodetic back azimuth to the origin point (in radians; if
                    not in parameters)
        alpha_g     Geodetic mean azimuth (in radians; if not in parameters)
        -----------------------------------------------------------------------

        Usage
        -----------------------------------------------------------------------
        Instantiation of a geodesic object will perform calculations based on
        provided parameters.

        Use of the geodesic.data() method provides parameters and calculated
        values.
        
        Manual calculation is also available through object attribute
        assignment and use of the geodesic.build() method.
        -----------------------------------------------------------------------
        """
        
        self.params = kwargs
        self.__build()


    def __build(self):

        self.params.update(self.__direct(**self.params))
        self.params.update(self.__inverse(**self.params))
        self.params.update(self.__geodetic_azimuth(**self.params))
        
        return

   
    def __direct(self, phi_0=None, lambda__0=None, alpha_0=None, s=None,
                 a=None, b=None, f=None, phi=None, lambda_=None, alpha=None,
                 **kwargs):

        # Performs direct Vincenty transformation based on provided origin
        # point, geodetic forward azimuth, geodesic, and ellipsoid parameters

        if (not missing.any(phi_0, lambda__0, alpha_0, s, a, b, f) and
            missing.any(phi, lambda_, alpha)):
            beta_0 = math.atan(math.tan(phi_0) * (1 - f))
            sin_beta_0, cos_beta_0 = math.sin(beta_0), math.cos(beta_0)
            tan_beta_0 = math.tan(beta_0)
            cos_alpha_0 = math.cos(alpha_0)
            tan_sigma_0 = tan_beta_0 / cos_alpha_0
            sin_alpha_0, cos_alpha_0 = math.sin(alpha_0), math.cos(alpha_0)
            sigma_0 = math.atan2(math.tan(beta_0), math.cos(alpha_0))
            sin_alpha_e = cos_beta_0 * sin_alpha_0
            u2 = (1 - sin_alpha_e**2) * ((a**2 - b**2) / b**2)
            k_0 = (math.sqrt(1 + u2) - 1) / (math.sqrt(1 + u2) + 1)
            A, B = (1 + 0.25 * k_0**2) / (1 - k_0), k_0 * (1 - 0.375 * k_0**2)
            counter = 0
            sigma = s / (b * A)
            while True:
                _2sigma_m = 2 * sigma_0 + sigma
                cos_2sigma_m = math.cos(_2sigma_m)
                sin_sigma = math.sin(sigma)
                sin2_sigma = sin_sigma**2
                cos_sigma = math.cos(sigma)
                delta_sigma = B * sin_sigma * (
                    cos_2sigma_m + B / 4 * (
                        cos_sigma * (-1 + 2 * cos_2sigma_m**2)) -
                    B / 6 * cos_2sigma_m * (-3 + 4 * sin2_sigma) *
                    (-3 + 4 * cos_2sigma_m**2))
                sigma_iter = delta_sigma + s / (b * A)
                if sigma == sigma_iter or counter == 1000:
                    sigma = sigma_iter
                    break
                sigma = sigma_iter
                counter += 1
            if not phi:
                phi = math.atan2(
                    (sin_beta_0 * cos_sigma + cos_beta_0 * sin_sigma *
                     cos_alpha_0),
                    (1 - f) * math.sqrt(
                        sin_alpha_e**2 + (sin_beta_0 * sin_sigma -
                                          cos_beta_0 * cos_sigma *
                                          cos_alpha_0)**2))
            lambda__diff_s = math.atan2(
                sin_sigma * sin_alpha_0,
                (cos_beta_0 * cos_sigma - sin_beta_0 *
                 sin_sigma * cos_alpha_0))
            C = ((f / 16) * cos_alpha_0**2 *
                 (4 + f * (4 - 3 * (1 - sin_alpha_e**2))))
            lambda__diff = (lambda__diff_s - (1 - C) * f * sin_alpha_e *
                            (sigma + C * sin_sigma *
                             (cos_2sigma_m + C * cos_sigma *
                              (-1 + 2 * cos_2sigma_m**2))))
            if not lambda_:
                lambda_ = lambda__diff + lambda__0
            if not alpha:
                alpha = math.atan2(sin_alpha_e,
                                   (-sin_beta_0 * sin_sigma + cos_beta_0 *
                                    cos_sigma * cos_alpha_0)) + math.pi
        kwargs['phi'], kwargs['lambda_'], kwargs['alpha'] = phi, lambda_, alpha
        
        return kwargs
                 

    def __inverse(self, phi=None, lambda_=None, phi_0=None, lambda__0=None,
                  a=None, b=None, f=None, alpha_0=None, alpha=None, s=None,
                  **kwargs):

        # Performs inverse Vincenty transformation based on provided origin and
        # destination points, along with ellipsoid parameters
        
        if (not missing.any(phi, lambda_, phi_0, lambda__0, a, b, f) and
            missing.any(alpha_0, alpha, s)):        
            beta_0 = math.atan(math.tan(phi_0) * (1 - f))
            beta = math.atan(math.tan(phi) * (1 - f))
            lambda__diff = lambda_ - lambda__0
            counter = 0
            lambda__diff_s = lambda__diff
            while True:
                sin_beta_0, cos_beta_0 = math.sin(beta_0), math.cos(beta_0)
                tan_beta_0 = math.tan(beta_0)
                sin_beta, cos_beta = math.sin(beta), math.cos(beta)
                tan_beta = math.tan(beta)
                sin2_sigma = ((cos_beta * lambda__diff_s)**2 +
                              (cos_beta_0 * sin_beta - sin_beta_0 *
                               cos_beta * math.cos(lambda__diff_s))**2)
                sin_sigma = math.sqrt(sin2_sigma)
                cos_sigma = (math.sin(beta_0) * math.sin(beta) +
                             math.cos(beta_0) * math.cos(beta) *
                             math.cos(lambda__diff_s))
                sigma = math.atan2(sin_sigma, cos_sigma)
                sin_alpha_e = ((cos_beta_0 * cos_beta *
                                math.sin(lambda__diff_s)) / sin_sigma)
                cos2_alpha_e = 1 - sin_alpha_e**2
                cos_2sigma_m = (cos_sigma - 2 * math.sin(beta_0) *
                                math.sin(beta) / (1 - sin_alpha_e**2))
                C = ((f / 16) * cos2_alpha_e *
                     (4 + f * (4 - 3 * cos2_alpha_e)))
                lambda__diff_s_iter = (
                    lambda__diff + (1 - C) * f * sin_alpha_e *
                    (sigma + C * sin_sigma *
                     (cos_2sigma_m + C * math.cos(sigma) *
                      (-1 + 2 * cos_2sigma_m**2))))
                if lambda__diff_s == lambda__diff_s_iter or counter == 1000:
                    lambda__diff_s = lambda__diff_s_iter
                    break
                lambda__diff_s = lambda__diff_s_iter
                counter += 1
            u2 = (1 - sin_alpha_e**2) * ((a**2 - b**2) / b**2)
            k_0 = (math.sqrt(1 + u2) - 1) / (math.sqrt(1 + u2) + 1)
            A = (1 + 0.25 * k_0**2) / (1 - k_0)
            B = k_0 * (1 - 0.375 * k_0**2)
            delta_sigma = B * sin_sigma * (
                cos_2sigma_m + B / 4 * (
                    cos_sigma * (-1 + 2 * cos_2sigma_m**2)) -
                B / 6 * cos_2sigma_m * (-3 + 4 * sin2_sigma) *
                (-3 + 4 * cos_2sigma_m**2))
            if not s:
                s = b * A * (sigma - delta_sigma)
            tan_alpha_0 = (cos_beta * math.sin(lambda__diff_s) /
                           (cos_beta_0 * sin_beta - sin_beta_0 * cos_beta *
                            math.cos(lambda__diff_s)))
            tan_alpha = (cos_beta_0 * math.sin(lambda__diff_s) /
                         (-sin_beta_0 * cos_beta + cos_beta_0 * sin_beta *
                          math.cos(lambda__diff_s)))
            delta_phi_f, delta_lambda__f = phi - phi_0, lambda_ - lambda__0
            delta_phi_r, delta_lambda__r = phi_0 - phi, lambda__0 - lambda_
            if not alpha_0:
                alpha_0 = math.atan(tan_alpha_0)
                if delta_lambda__f < 0:
                    alpha_0 += 2 * math.pi
                else:
                    alpha_0 += math.pi
            if not alpha:
                alpha = math.atan(tan_alpha)
                if delta_lambda__r > 0:
                    alpha += math.pi
                if alpha < 0:
                    alpha +=  2 * math.pi
        kwargs['s'], kwargs['alpha_0'], kwargs['alpha'] = s, alpha_0, alpha
        
        return kwargs


    def __geodetic_azimuth(self, alpha_0=None, alpha=None, alpha_g=None,
                           **kwargs):
        
        # Produces geodetic mean azimuth based on forward and back geodetic
        # azimuths
        
        if not missing.any(alpha_0, alpha) and not alpha_g:
            alpha_g = (alpha_0 + alpha - math.pi) / 2
        kwargs['alpha_g'] = alpha_g

        return kwargs


    def data(self):

        return self.params


    def build(self):

        self.__build()
        
        return
