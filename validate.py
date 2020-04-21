"""
Cadastral Measurement Management and Maintenance in Three Dimensions
-------------------------------------------------------------------------------
Title:          validate
Author:         Stanton K. Nielson, GIS Specialist
                BLM Wyoming High Desert District/Elmhurst College
Date:           April 21, 2020
Version:        1.0
Description:    Provides validation capabilities to test lists, tuples, sets,
                and dictionaries for missing values and simplify syntax for
                pervasive if testing
-------------------------------------------------------------------------------
"""


class missing(object):
    

    def any(*args):

        # Returns True if any value in a sequence is None
        
        return None in [i for i in args]
    

    def all(*args):

        # Returns True if all values in a sequence are None
        
        return [i for i in args] == [None for i in args]


    def any2(**kwargs):

        # Returns True if any value in a dictionary is None
        
        return None in [i for i in kwargs.values()]


    def all2(**kwargs):

        # Returns True if all values in a dictionary are None
        
        return [i for i in kwargs.values()] == [None for i in kwargs.values()]




    


    

    
