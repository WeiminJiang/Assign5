# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 13:55:14 2017
"""
from __future__ import division   

#remove previous workspace
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
        
 
from schemes import *
from initial_conditions import *
from ErrorAnalysis import *
import matplotlib.pyplot as plt   
        
def main():

    #constants
    x_min = 0
    x_max = 1
    nx = 40
    nt = 200*32
    dx = (x_max - x_min)/nx
    K = 1e-3
    u = 1
    T = 0.625*32
    
    #derived constants
    dt = T / nt
    c = dt * u / dx
    #print('dt=', dt, 'c=', c, 'nt =', nt)
    d = K*dt/dx**2 #Non-dimensional diffusion coefficent

    #initial conditions
    x = create_x(x_min, dx, nx)
    phi1 = initial_conditions_1(x)   
    phi2 = initial_conditions_2(x)

    #analytic solution
    phiAnalytic1 = initial_conditions_1((x-u*T)%(x_max-x_min))
    phiAnalytic2 = initial_conditions_2((x-u*T)%(x_max-x_min))
    
    #FTCS
    phiFTCS1 = FTCS(phi1, c, nt)
    phiFTCS2 = FTCS(phi2, c, nt)

    #CTCS    
    phiCTCS1 = CTCS(phi1, c, nt)
    phiCTCS2 = CTCS(phi2, c, nt)

    #FTBS
    phiFTBS1 = FTBS(phi1, c, nt)
    phiFTBS2 = FTBS(phi2, c, nt)

    #CTBS
    phiCTBS1 = CTBS(phi1, c, nt)
    phiCTBS2 = CTBS(phi2, c, nt)

    #CTCS with artificial diffusion
    phiCTCS_art_diff1 = CTCS_art_diff(phi1, c, d, nt)
    phiCTCS_art_diff2 = CTCS_art_diff(phi2, c, d, nt)

    #Warming and Beam
    phiWB1 = WB(phi1, c, nt)
    phiWB2 = WB(phi2, c, nt)
    
    #TVD
    phiTVD1 = TVD(phi1, c, nt)[0]
    phiTVD2 = TVD(phi2, c, nt)[0]

    #Lax_Wendroff
    phiLax_Wendroff1 = Lax_Wendroff(phi1, c, nt)
    phiLax_Wendroff2 = Lax_Wendroff(phi2, c, nt)
    
    

#    #plot all schemes with the first set of initial conditions
#    font = {'size': 10}
#    plt.rc('font', **font)
#    plt.figure(1)
#    plt.clf()
#    plt.ion()
#    plt.plot(x, phi1, label = 'Initial Conditions', color='b')
#    plt.plot(x, phiFTCS1, label='FTCS', color='r')
#    plt.plot(x, phiCTCS1, label='CTCS', color='g')
#    plt.plot(x, phiFTBS1, label='FTBS', color='c')
#    plt.plot(x, phiCTBS1, label='CTBS', color='m')
#    plt.plot(x, phiAnalytic1, label='Analytic Solution', color='k')
#    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    plt.axhline(0, linestyle=':', color='black')
#    plt.xlabel('$x$')
#    plt.ylabel(r'$\phi$')
#    plt.tight_layout()
#    #plt.savefig("C:\Users\Joshua\Desktop\initial_conditions_wave.pdf", \
#    #            bbox_extra_artists=(lgd,), bbox_inches='tight')
#    plt.show()
#    
#    #plot all schemes with the second set of initial conditions
#    plt.figure(2)
#    plt.clf()
#    plt.ion()
#    plt.plot(x, phi2, label = 'Initial Conditions', color='b')
#    plt.plot(x, phiFTCS2, label='FTCS', color='r')
#    plt.plot(x, phiCTCS2, label='CTCS', color='g')
#    plt.plot(x, phiFTBS2, label='FTBS', color='c')
#    #plt.plot(x, phiCTBS2, label='CTBS', color='m')
#    plt.plot(x, phiAnalytic2, label='Analytic Solution', color='k')
#    plt.xlim([0.0, 1.0])
#    #plt.ylim([-0.4 ,1.2])
#    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    plt.axhline(0, linestyle=':', color='black')
#    plt.xlabel('$x$')
#    plt.ylabel(r'$\phi$')
#    plt.tight_layout()
#    #plt.savefig("C:\Users\Joshua\Desktop\initial_conditions_square.pdf", \
#    #            bbox_extra_artists=(lgd,), bbox_inches='tight')
#    plt.show() 
#   
#    #plot all schemes from chap11 with first set of initial conditions
#    plt.figure(3)
#    plt.clf()
#    plt.plot(x, phiAnalytic1, label='Analytic Solution', color='k')
#    plt.plot(x, phiWB1, label = 'WB', color = 'r')
#    plt.plot(x, phiTVD1, label = 'TVD', color = 'b')
#    plt.plot(x, phiCTCS_art_diff1, label='CTCSAD', color='y')
#    plt.plot(x, phiLax_Wendroff1, label='LW', color='g')
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    plt.xlabel('$x$')
#    plt.ylabel(r'$\phi$')
#    #plt.savefig("C:\Users\Joshua\Desktop\initial_conditions_wave_part2.pdf", \
#    #            bbox_extra_artists=(lgd,), bbox_inches='tight')
#    plt.show()
#    
#   
#    #plot all schemes from chap11 with second set of initial conditions
#    font = {'size': 10}
#    plt.rc('font', **font)
#    plt.figure(1)
#    plt.clf()
#    plt.plot(x, phiAnalytic2, label='Analytic Solution', color='k')
#    #plt.plot(x, phiWB2, label = 'WB', color = 'r')
#    #plt.plot(x, phiTVD2, label = 'TVD', color = 'b')
#   # plt.plot(x, phiCTCS_art_diff2, label='CTCSAD', color='y')
#    plt.plot(x, phiLax_Wendroff2, label='LW', color='g')
#    plt.plot(x, phiCTCS2, label='CTCS', color='r')
#    plt.plot(x, phiFTBS2, label='FTBS', color='c')
#    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#    plt.xlabel('$x$')
#    plt.ylabel(r'$\phi$')
##   plt.savefig("C:\Users\Joshua\Desktop\initial_conditions_square_part2.pdf", \
##                bbox_extra_artists=(lgd,), bbox_inches='tight')
#    plt.show()
#    
#    
#    #plot for part 2 question 3 by 25819903 
#    plt.figure(2)
#    plt.clf()
#    plt.plot(TVD(phi1, c, nt)[1])
#    plt.xlabel('time step')
#    plt.ylabel('TV')
#    plt.savefig('3.pdf')
#    plt.show()
    
  
    
    
    # Experiment for part2 question 3 by 258362
    print('T = ', T)
    print ('Cumulative difference between scheme and analytic solution for the first set of initial conditions:')
    print ('LW = '+ str(L2ErrorNorm(phiLax_Wendroff1, phiAnalytic1)[0]))
    print ('CTCS = '+ str(L2ErrorNorm(phiCTCS1, phiAnalytic1)[0]))
    print ('FTBS = '+ str(L2ErrorNorm(phiFTBS1, phiAnalytic1)[0]))
    print ("")
    print ('Cumulative difference between scheme and analytic solution for the second set of initial conditions:')
    print ('LW = '+ str(L2ErrorNorm(phiLax_Wendroff2, phiAnalytic2)[0]))
    print ('CTCS = '+ str(L2ErrorNorm(phiCTCS2, phiAnalytic2)[0]))
    print ('FTBS = '+ str(L2ErrorNorm(phiFTBS2, phiAnalytic2)[0]))
    
  
    
    
    
    
    
'''
   
    #values for part 2 question 3 by 25819903
    print ('Cumulative difference between scheme and analytic solution for the first set of initial conditions:')
    print ('TVD = '+ str(L2ErrorNorm(phiTVD1, phiAnalytic1)[0]))
    print ('WB = ' + str(L2ErrorNorm(phiWB1, phiAnalytic1)[0]))
    print ('CTCS = '+ str(L2ErrorNorm(phiCTCS1, phiAnalytic1)[0]))
    print ('LW = '+ str(L2ErrorNorm(phiLax_Wendroff1, phiAnalytic1)[0]))
    print ("")
    print ('Cumulative difference between scheme and analytic solution for the second set of initial conditions:')
    print ('TVD = '+ str(L2ErrorNorm(phiTVD2, phiAnalytic2)[0]))
    print ('WB = ' + str(L2ErrorNorm(phiWB2, phiAnalytic2)[0]))
    print ('CTCS = '+ str(L2ErrorNorm(phiCTCS2, phiAnalytic2)[0]))
    print ('LW = '+ str(L2ErrorNorm(phiLax_Wendroff2, phiAnalytic2)[0]))
    print ("")
    
    #values for part 2 question 4 by 25819903
    print ('Integration of solution for the first set of initial conditions')
    print ('Initial Conditions = ' + str(np.trapz(phi1, x, dx)))
    print ('Analytic = ' + str(np.trapz(phiAnalytic1, x, dx)))
    print ('FCTS = ' + str(np.trapz(phiFTCS1, x, dx)))
    print ('CTCS = ' + str(np.trapz(phiCTCS1, x, dx)))  
    print ('FTBS = ' + str(np.trapz(phiFTBS1, x, dx)))
    print ('CTBS = ' + str(np.trapz(phiCTBS1, x, dx)))
    print ('WB = ' + str(np.trapz(phiWB1, x, dx)))
    print ('TVD = ' + str(np.trapz(phiTVD1, x, dx)))
    print ('CTCSAD = ' + str(np.trapz(phiCTCS_art_diff1, x, dx)))
    print ('LW = ' + str(np.trapz(phiLax_Wendroff1, x, dx)))
    print ("")
    print ('Integration of solution for the second set of initial conditions')
    print ('Initial Conditions = ' + str(np.trapz(phi2, x, dx)))
    print ('Analytic = ' + str(np.trapz(phiAnalytic2, x, dx)))
    print ('FCTS = ' + str(np.trapz(phiFTCS2, x, dx)))
    print ('CTCS = ' + str(np.trapz(phiCTCS2, x, dx)))  
    print ('FTBS = ' + str(np.trapz(phiFTBS2, x, dx)))
    print ('CTBS = ' + str(np.trapz(phiCTBS2, x, dx)))
    print ('WB = ' + str(np.trapz(phiWB2, x, dx)))
    print ('TVD = ' + str(np.trapz(phiTVD2, x, dx)))
    print ('CTCSAD = ' + str(np.trapz(phiCTCS_art_diff2, x, dx)))
    print ('LW = ' + str(np.trapz(phiLax_Wendroff2, x, dx)))
'''
main()