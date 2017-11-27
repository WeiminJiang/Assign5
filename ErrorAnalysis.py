
import numpy as np
def L2ErrorNorm(phi, phiExact):
    '''
    Calculate the L2 error norm (RMS error)
    25836326
    '''

    # calculate the error and the error norms
    phiError = phi - phiExact
    L2 = np.sqrt(sum(phiError ** 2) / sum(phiExact ** 2))

    return L2, phiError

def total_variation(phi):
    '''
    Calculate the total variation
    25836326
    '''
    tv = 0
    
    for ix in range(len(phi)-1):
        tv = tv + np.abs(phi[ix+1] - phi[ix])
        
    return tv

def boundedness(phi):
    '''
    Test the boundeness
    25836326
    '''
    nx = len(phi)
    # if test = True: bounded, return 0
    # if test = False: unbouned, return 1
    test = True
    for ix in range(nx):
        if phi[ix] - 0. < 0. or phi[ix] - 1. > 0.: 
            test = False
            break

    if test == True:
        return 0
    else:
        return 1

def boundedness_step2(bb):
    
    nt = len(bb)
    test = True
    for it in range(nt):
        if bb[it] > 0.: 
            test = False
            break

    if test == True:
        return 'bounded'
    else:
        return 'unbounded'
    
    
    

