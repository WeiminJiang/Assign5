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
    
    x_min = 0
    x_max = 1
    nx = 40
    nt = 50
    dx = (x_max - x_min)/nx
    K = 1e-3
    u = 1
    T = 0.125
    
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
    phiFTCS1 = FTCS(phi1, c, nt)[0]
    phiFTCS2 = FTCS(phi2, c, nt)[0]

    #CTCS    
    phiCTCS1 = CTCS(phi1, c, nt)[0]
    phiCTCS2 = CTCS(phi2, c, nt)[0]

    #FTBS
    phiFTBS1 = FTBS(phi1, c, nt)[0]
    phiFTBS2 = FTBS(phi2, c, nt)[0]

    #CTBS
    phiCTBS1 = CTBS(phi1, c, nt)[0]
    phiCTBS2 = CTBS(phi2, c, nt)[0]

    #CTCS with artificial diffusion
    phiCTCS_art_diff1 = CTCS_art_diff(phi1, c, d, nt)[0]
    phiCTCS_art_diff2 = CTCS_art_diff(phi2, c, d, nt)[0]

    #Warming and Beam
    phiWB1 = WB(phi1, c, nt)[0]
    phiWB2 = WB(phi2, c, nt)[0]
    
    #TVD
    phiTVD1 = TVD(phi1, c, nt)[0]
    phiTVD2 = TVD(phi2, c, nt)[0]

    #Lax_Wendroff
    phiLax_Wendroff1 = Lax_Wendroff(phi1, c, nt)[0]
    phiLax_Wendroff2 = Lax_Wendroff(phi2, c, nt)[0]
     
    #plot all schemes with the first set of initial conditions
    font = {'size': 15}
    plt.rc('font', **font)
    plt.figure(figsize=(10,8))
    plt.clf()
    plt.ion()
    plt.plot(x, phi1, label = 'Initial Conditions', color='b')
    plt.plot(x, phiFTCS1, label='FTCS', color='r')
    plt.plot(x, phiCTCS1, label='CTCS', color='g')
    plt.plot(x, phiFTBS1, label='FTBS', color='c')
    plt.plot(x, phiCTBS1, label='CTBS', color='m')
    plt.plot(x, phiAnalytic1, label='Analytic Solution', color='k')
    lgd = plt.legend(loc='center left', bbox_to_anchor=(0.65, 0.8))
    plt.axhline(0, linestyle=':', color='black')
    plt.xlabel('$x$')
    plt.ylabel(r'$\phi$')
    plt.title('Initial Condition 1')
    plt.savefig('1a.jpg')
    plt.tight_layout()
    plt.show()
    
    #plot all schemes with the second set of initial conditions
    plt.figure(figsize=(10,8))
    plt.clf()
    plt.ion()
    plt.plot(x, phi2, label = 'Initial Conditions', color='b')
    plt.plot(x, phiFTCS2, label='FTCS', color='r')
    plt.plot(x, phiCTCS2, label='CTCS', color='g')
    plt.plot(x, phiFTBS2, label='FTBS', color='c')
    #plt.plot(x, phiCTBS2, label='CTBS', color='m')
    plt.plot(x, phiAnalytic2, label='Analytic Solution', color='k')
    plt.xlim([0.0, 1.0])
    #plt.ylim([-0.4 ,1.2])
    lgd = plt.legend(loc='center left', bbox_to_anchor=(0.65, 0.8))
    plt.axhline(0, linestyle=':', color='black')
    plt.xlabel('$x$')
    plt.ylabel(r'$\phi$')
    plt.title('Initial Condition 2')
    plt.tight_layout()
    plt.savefig('1b.jpg')
    plt.show() 
   
    #plot all schemes from chap11 with first set of initial conditions
    plt.figure(figsize=(10,8))
    plt.clf()
    plt.plot(x, phiAnalytic1, label='Analytic Solution', color='k')
    plt.plot(x, phiWB1, label = 'WB', color = 'r')
    plt.plot(x, phiTVD1, label = 'TVD', color = 'b')
    plt.plot(x, phiCTCS_art_diff1, label='CTCSAD', color='y')
    plt.plot(x, phiLax_Wendroff1, label='LW', color='g')
    plt.plot(x, phiCTCS1, label='CTCS', color='m')
    plt.plot(x, phiFTBS1, label='FTBS', color='c')
    plt.legend(loc='center left', bbox_to_anchor=(0.65, 0.8))
    plt.xlabel('$x$')
    plt.ylabel(r'$\phi$')
    plt.savefig('2a.jpg')
    plt.show()
    
   
    #plot all schemes from chap11 with second set of initial conditions
    font = {'size': 15}
    plt.rc('font', **font)
    plt.figure(figsize=(10,8))
    plt.clf()
    plt.plot(x, phiAnalytic2, label='Analytic Solution', color='k')
    plt.plot(x, phiWB2, label = 'WB', color = 'r')
    plt.plot(x, phiTVD2, label = 'TVD', color = 'b')
    plt.plot(x, phiCTCS_art_diff2, label='CTCSAD', color='y')
    plt.plot(x, phiLax_Wendroff2, label='LW', color='g')
    plt.plot(x, phiCTCS2, label='CTCS', color='m')
    plt.plot(x, phiFTBS2, label='FTBS', color='c')
    plt.legend(loc='center left', bbox_to_anchor=(0.63, 0.8))
    plt.xlabel('$x$')
    plt.ylabel(r'$\phi$')
    plt.savefig('2b.jpg')
    plt.show()
    
    
    # Experiment for part 2 question 3 by 25836326 
    font = {'size': 15}
    plt.rc('font', **font)
    plt.figure(figsize=(10,8))
    plt.clf()
    plt.plot(Lax_Wendroff(phi2, c, nt)[1], label='LW', color='g')
    plt.plot(CTCS(phi2, c, nt)[1], label='CTCS', color='y')
    plt.plot(FTBS(phi2, c, nt)[1], label='FTBS', color='r')
    plt.axhline(2, linestyle=':', color='black')
    plt.xlabel('time_steps')
    plt.ylabel('Total_variation')
    plt.legend()
    plt.savefig('3.jpg')
    plt.show()    

    print('Total variation of CTCS at last tiem step', total_variation(phiCTCS2))
    print('Total variation of LW at last tiem step', total_variation(phiLax_Wendroff2))

    
    # Experiment for part 2 question 3 by 25836326
    print ("")
    print ('Cumulative difference between scheme and analytic solution for initial condition 1:')
    print ('LW = '+ str(L2ErrorNorm(phiLax_Wendroff1, phiAnalytic1)[0]))
    print ('CTCS = '+ str(L2ErrorNorm(phiCTCS1, phiAnalytic1)[0]))
    print ('FTBS = '+ str(L2ErrorNorm(phiFTBS1, phiAnalytic1)[0]))
    print ("")
    print ('Cumulative difference between scheme and analytic solution for initial condition 2:')
    print ('LW = '+ str(L2ErrorNorm(phiLax_Wendroff2, phiAnalytic2)[0]))
    print ('CTCS = '+ str(L2ErrorNorm(phiCTCS2, phiAnalytic2)[0]))
    print ('FTBS = '+ str(L2ErrorNorm(phiFTBS2, phiAnalytic2)[0]))
    
    # Experiment for part2 question 4 by 25836326
    print ("")
    print('The boudeness of all schemes in inital condition 1')
    print ('FTCS is '+ str(boundedness_step2(FTCS(phi1, c, nt)[1])))
    print ('CTCS is '+ str(boundedness_step2(CTCS(phi1, c, nt)[2])))
    print ('FTBS is '+ str(boundedness_step2(FTBS(phi1, c, nt)[2])))
    print ('CTBS is '+ str(boundedness_step2(CTBS(phi1, c, nt)[1])))
    print ('WB is '+ str(boundedness_step2(WB(phi1, c, nt)[1])))
    print ('TVD is '+ str(boundedness_step2(TVD(phi1, c, nt)[2])))
    print ('CTCSAD is '+ str(boundedness_step2(CTCS_art_diff(phi1, c, d, nt)[1])))
    print ('LW is '+ str(boundedness_step2(Lax_Wendroff(phi1, c, nt)[2])))
  
    print ("")
    print('The boudeness of all schemes in inital condition 2')
    print ('FTCS is '+ str(boundedness_step2(FTCS(phi2, c, nt)[1])))
    print ('CTCS is '+ str(boundedness_step2(CTCS(phi2, c, nt)[2])))
    print ('FTBS is '+ str(boundedness_step2(FTBS(phi2, c, nt)[2])))
    print ('CTBS is '+ str(boundedness_step2(CTBS(phi2, c, nt)[1])))
    print ('WB is '+ str(boundedness_step2(WB(phi2, c, nt)[1])))
    print ('TVD is '+ str(boundedness_step2(TVD(phi2, c, nt)[2])))
    print ('CTCSAD is '+ str(boundedness_step2(CTCS_art_diff(phi2, c, d, nt)[1])))
    print ('LW is '+ str(boundedness_step2(Lax_Wendroff(phi2, c, nt)[2])))
    
    

    
    



main()