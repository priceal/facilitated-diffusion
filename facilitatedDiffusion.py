"""
faciliated diffusion.
only include siding at this time.

"""

import math
import random

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
 