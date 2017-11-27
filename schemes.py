# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 14:15:53 2017
"""

import numpy as np
from ErrorAnalysis import total_variation
from ErrorAnalysis import boundedness


def FTBS(phiOld, c, nt):
    '''
    25836326
    '''
    nx = len(phiOld)
    TV = np.zeros(nt + 1)
    TV[0] = total_variation(phiOld)
    bb = np.zeros(nt)

    phi = np.zeros(len(phiOld), dtype='float')

    for it in range(int(nt)):
        for j in range(0, nx):
            phi[j] = phiOld[j] - c * (phiOld[(j) % nx] - phiOld[(j - 1) % nx])
        phiOld = phi.copy()
        TV[it + 1] = total_variation(phiOld)
        bb[it] = boundedness(phiOld)
    return phiOld, TV, bb


def FTCS(phiOld, c, nt):
    '''
    23009181
    '''
    nx = len(phiOld)
    bb = np.zeros(nt)

    # New time-step array for phi#
    phi = np.zeros(len(phiOld), dtype='float')

    # FTCS for all time steps#
    for it in range(nt):
        for j in range(0, nx):
            phi[j] = phiOld[j] - 0.5 * c * (phiOld[(j + 1) % nx] - phiOld[(j - 1) % nx])
        phiOld = phi.copy()
        bb[it] = boundedness(phiOld)

    return phiOld, bb


def CTCS(phi, c, nt):
    '''
    Advection of profile in phi using CTCS using non-dimensional courant number
    coefficient, c
    
    23012207
    '''
    nx = len(phi)
    TV = np.zeros(nt + 1)
    bb = np.zeros(nt)
    TV[0] = total_variation(phi)

    # error handling
    if (nx <= 0):
        raise ValueError("The initial array must have a positive length!")

    try:
        value = float(c)
    except ValueError:
        print("c must be a valid number!")

    if nt != int(nt):
        raise ValueError("nt must be an integer!")

    # new time-step array for phiNew
    phiNew = np.zeros_like(phi)

    # FTCS for first time-step
    for x in range(nx):
        phiNew[x] = phi[x] - 0.5*c*(phi[(x+1)%nx] - phi[(x-1)%nx])

    phiOld = phi.copy()
    phi = phiNew.copy()
    
    # CTCS for all other time-steps
    for it in range(1, int(nt)):
        for x in range(nx):
            phiNew[x] = phiOld[x] - c*(phi[(x+1)%nx] - phi[(x-1)%nx])
    
        phiOld = phi.copy()
        phi = phiNew.copy()
        TV[it + 1] = total_variation(phi)
        bb[it] = boundedness(phi)

    return phi, TV, bb


def CTBS(phi, c, nt):
    '''
    CTBS Advection Scheme with intial conditions phi, Courant number c
    and number of time steps nt
    25819903
    '''

    # Exception Handling

    # Store length of phi to aviod recalculating each time
    nx = len(phi)
    bb = np.zeros(nt)

    if (nx == 0):
        raise Exception('Phi has length of 0. It needs to have a positive length.')

    if nt < 0:
        raise Exception('The number of timesteps is negative. It should be positive.')
    if nt != int(nt):
        raise Exception('The number of timesteps is not an integer.')

    flag = False
    for i in range(nx):
        if phi[i] != 0:
            flag = True
            break
    if flag == False:
        raise Exception('Phi is all zeros. This shows that it is either \
        uninitialized or the array is not the one intended to be used')

    # Define arrays to store calculated values
    phi_new = np.zeros_like(phi)
    phi_n = np.zeros_like(phi)
    phi_n_minus_1 = phi.copy()

    # FTCS for first time step
    for j in range(nx):
        phi_n[j] = phi_n_minus_1[j % nx] - 0.5 * c * (phi_n_minus_1[(j + 1) % nx]\
        - phi_n_minus_1[(j - 1) % nx])

    # CTBS for the rest of the time steos
    for n in range(nt):
        for j in range(nx):
            phi_new[j] = phi_n_minus_1[j % nx] - 2 * c * (phi_n[j % nx] - \
            phi_n[(j - 1) % nx])
        phi_n_minus_1 = phi_n.copy()
        phi_n = phi_new.copy()
        bb[n] = boundedness(phi_new)
    return phi_new, bb


def CTCS_art_diff(phi, c, d, nt):
    '''
    Advection of profile in phi using CTCS with artificial diffusion using
    non-dimensional courant number coefficient, c and diffusion number, d
    
    23012207
    '''

    nx = len(phi)
    bb = np.zeros(nt)

    # error handling
    if (nx <= 0):
        raise ValueError("The initial array must have a positive length!")

    try:
        value = float(c)
    except ValueError:
        print("c must be a valid number!")

    try:
        value = float(d)
    except ValueError:
        print("d must be a valid number!")

    if nt != int(nt):
        raise ValueError("nt must be an integer!")

    # new time-step array for phi
    phiNew = np.zeros_like(phi)

    # FTCS for advection and diffusion for first time-step
    for x in range(nx):
        phiNew[x] = phi[x] - 0.5*c*(phi[(x+1)%nx] - phi[(x-1)%nx]) \
        + 2*d*(phi[(x+1)%nx] - 2*phi[x%nx] + phi[(x-1)%nx])

    phiOld = phi.copy()
    phi = phiNew.copy()

    # CTCS for advection, FTCS for diffusion for all other time-steps
    for it in range(1, int(nt)):
        for x in range(nx):
            phiNew[x] = phiOld[x] - c*(phi[(x+1)%nx] - phi[(x-1)%nx]) \
            + 2*d*(phiOld[(x+1)%nx] - 2*phiOld[x%nx] + phiOld[(x-1)%nx])

        # set arrays ready for next time-step
        phiOld = phi.copy()
        phi = phiNew.copy()
        bb[it] = boundedness(phi)

    return phi, bb


def TVD(phi, c, nt):
    '''
    Total Variation Dimimnishing Advection Scheme with initial conditions phi,
    Courant number c and number of time steps nt
    25819903
    '''

    nx = len(phi)
    TV = np.zeros(nt + 1)
    bb = np.zeros(nt)

    # Error Handling
    if (nx == 0):
        raise ValueError('Phi has length of 0. It needs to have a positive\
        length.')

    if nt < 0:
        raise Exception('The number of timesteps is negative. It should be\
        positive.')
    if nt != int(nt):
        raise TypeError('The number of timesteps is not an integer.')

    flag = False
    for i in range(nx):
        if phi[i] != 0:
            flag = True
            break
    if flag == False:
        raise Exception('Phi is all zeros. This shows that it is either\
        uninitialized or the array is not the one intended to be used')

    phi_next = np.zeros_like(phi)

    TV[0] = total_variation(phi)
    for n in xrange(nt):
        for j in xrange(nx):
            lin_func_plus = limiter_function(phi, j)
            phi_h_plus = 0.5 * (1 + c) * phi[j] + 0.5 * (1 - c) * phi[(j + 1) % nx]
            phi_l_plus = 0
            if c >= 0:
                phi_l_plus = phi[j]
            else:
                phi_l_plus = phi[(j + 1) % nx]

            lin_func_minus = limiter_function(phi, (j - 1) % nx)
            phi_h_minus = 0.5 * (1 + c) * phi[(j - 1) % nx] + 0.5 * (1 - c) * phi[j]
            phi_l_minus = 0
            if c >= 0:
                phi_l_minus = phi[(j - 1) % nx]
            else:
                phi_l_minus = phi[(j)]

            phi_next[j] = phi[j] - c * \
                                   (
                                       lin_func_plus * phi_h_plus + \
                                       (1 - lin_func_plus) * phi_l_plus
                                       - (lin_func_minus * phi_h_minus + \
                                          (1 - lin_func_minus) * phi_l_minus)
                                   )

        phi = phi_next.copy()
        TV[n + 1] = total_variation(phi)
        bb[n] = boundedness(phi)

    return phi_next, TV, bb


def limiter_function(phi, j):
    '''
    Van Leer Limiter function for position j in phi
    25819903
    '''

    nx = len(phi)

    num = phi[j] - phi[(j - 1) % nx]
    denom = phi[(j + 1) % nx] - phi[j]

    r = 0
    if np.abs(denom) < 1e-6:
        r = 0
    else:
        r = num / denom

    return ((r + np.abs(r)) / (1 + np.abs(r)))


def WB(phiOld, c, nt):
    '''
    23009181
    '''
    nx = len(phiOld)
    bb = np.zeros(nt)

    phi = np.zeros(len(phiOld), dtype='float')

    # Warming and beam for every time step
    for it in range(nt):
        for j in range(0, nx):
            phi[j] = phiOld[j] - 0.5 * c * ((3 - c) * phiOld[j] -
                                            (1 - c) * phiOld[(j - 1) % nx] - (3 - c) * phiOld[(j - 1) % nx] +
                                            (1 - c) * phiOld[(j - 2) % nx])

        phiOld = phi.copy()
        bb[it] = boundedness(phiOld)
    return phiOld, bb


def Lax_Wendroff(phiOld, c, nt):
    '''
    25836326
    '''
    nx = len(phiOld)
    phi = np.zeros(len(phiOld), dtype='float')
    TV = np.zeros(nt + 1)
    bb = np.zeros(nt)
    TV[0] = total_variation(phiOld)

    # Testing Code
    if nx <= 0:
        raise ValueError('The length of Argument phi should be > 0.')

    test = False
    for ix in range(nx):
        if phiOld[ix] != 0:
            test = True
            break
    if test == False:
        raise Exception('Argument phi should not always be zero at all points.')

#    if not (isinstance(nt, int)):
#        raise TypeError('Argument nt should be an integer.')

#    if nt <= 0:
#        raise ValueError('Argument nt should be > 0.')

    if not (isinstance(c, float)):
        raise TypeError('Argument c should be a float.')

    # All time-steps
    for it in range(int(nt)):
        for j in range(nx):
            phi[j] = phiOld[j] - 0.5 * c * ((1.0 + c)*phiOld[(j) % nx] + (1.0 - c)*phiOld[(j + 1) % nx]\
                                            - (1.0 + c)*phiOld[(j - 1) % nx] - (1.0 - c)*phiOld[(j) % nx])

        phiOld = phi.copy()
        TV[it + 1] = total_variation(phi)
        bb[it] = boundedness(phiOld)
    return phiOld, TV, bb
