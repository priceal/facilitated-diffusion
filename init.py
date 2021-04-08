# initialization lines for facilitated diffusion environment
# loads all modules needed in scripts
# 

randomSeed = 48572389572

try:
    math
except NameError:  
    print(' importing math ...')
    import math
    
try:
    plt
except NameError:  
    print(' importing pylab as plt ...')
    import pylab as plt
    
try:
    np
except NameError:
    print(' importing numpy as np ...')
    import numpy as np
    
try:
    fd
except NameError:
    print(' importing facilitated-diffusion as fd ...')
    import facilitatedDiffusion as fd
    
try:
    FDF
except NameError:
    print(' importing FDFunctions as FDF ...')
    import FDFunctions as FDF
    print(' seeding FDFunctions random number generator ...')
    FDF.randseed(randomSeed)
    FDF.viseed()
    
