import numpy as np

def L2ErrorNorm(phi, phiExact):
    "Calculates the L2 error norm (RMS error)"
    '''
    25836326
    '''

    # calculate the error and the error norms
    phiError = phi - phiExact
    L2 = np.sqrt(sum(phiError ** 2) / sum(phiExact ** 2))

    return L2, phiError
    
def total_variation(phi):
    '''
    Return the Total Variation in phi
    25819903
    '''
    counter = 0
    
    for i in xrange(len(phi)-1):
        counter += np.abs(phi[i+1] - phi[i])
        
    return counter
    