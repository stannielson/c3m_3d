"""
Cadastral Measurement Management and Maintenance in Three Dimensions
----------------------------------------------------------------------------
Title:          c3m_3d
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    This module performs geodetic calculations for the
                management and maintenacence of cadastral measurements,
                specific to cadastral benchmarks. The methods rely on
                established algorithms to produce angular and metric
                measurements to quantify cadastral data for surveying
                applications.

                The module functions based on provided parameters and use of
                objects that fill in additional parameters when their build
                methods are called.  With these calls, calculated results
                allow for use as parameters in other calls.  The primary
                objective of the module is to provide automated capabilities
                for cadastral calculations for use in operational cadastral
                surveying applications.

                WARNING: while the design of this module has limitations in
                computation due to hardware limits (including rounding
                errors at high precision), the precision of input parameters
                will affect the precision and accuracy of calculated values.
                The recommended use of this module is with high-precision
                parameters to mitigate inaccuracy.
----------------------------------------------------------------------------
"""
