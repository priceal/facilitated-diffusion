"""
faciliated diffusion.
only include siding at this time.

"""

import math
import random
import FDFunctions as FDF

def slide():
    """
    Calculates a sliding event. Unit of time is the inverse of the off rate.
    Unit of distance is defined by the RMSD in one time unit.

    Returns
    -------
    dx : float
        The displacement during the sliding event.
    dt : float
        The residence time.

    """
    dt = -math.log(random.random())
    
    return math.sqrt(dt)*random.gauss(0,1), dt

def hop(ro,Rmax):
    """
    Parameters
    ----------
    ro : float
        Initial distance from central axis (capture surface at r=1).
    Rmax : float
        Radius of escape surface in units of capture radius.

    Returns
    -------
    x, y, z : tuple of floats
        Final coordinates of the hop (just inside/outside capture/escape surface)
    t : float
        Final time at which hop ends.
        
    """
    x, y, t = FDF.woc(ro, Rmax)
    z = math.sqrt(t)*FDF.boxmuller()[0]
    
    return x, y, z, t

