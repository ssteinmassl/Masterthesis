#!/usr/bin/env python

import os
import sys
from math import *
import string
import numpy as np
import thread
import time
import Tkinter
from scipy import constants as const
sys.path.append('/home/simonste/research/xrb/scripts/modules')
from my_functions import *

def usage():
    print __doc__


def FluxToMag(flux,filtcorr,lambdaeff):
    mag =  -2.5*log10(float(flux)*filtcorr)-2.402-5.0*log10(lambdaeff)
    return mag

def ABtoFnu(AB):
    fnu = 10**(-(AB+48.6)/2.5)
    return fnu


def ParseDefault(default):
    f = open(default, 'r')
    lines = f.readlines()
    f.close()

    parameters = {}
    parlist = []
    vallist = []
    for line in lines:
        try:
            par = line.split('=')[0]
            val = line.split('=')[1]
            if par =='COMMENT':
                continue
            if par =='BANDPASS':
                continue
            parameters[par]=val
        except:
            pass
    return parameters



def WriteParfile(dict, file):
    f = open(file, 'w')
    for par in dict:
        line = '%s= %s\n'%(par, dict[par])
        f.write(line)

    f.write('BANDPASS=       FILTER  g\n')
    f.write('BANDPASS=       FILTER  r\n')
    f.write('BANDPASS=       FILTER  i\n')
    f.write('BANDPASS=       FILTER  z\n')
    f.write('BANDPASS=       FILTER  J\n')
    f.write('BANDPASS=       FILTER  H\n')
    f.write('BANDPASS=       FILTER  K\n')

    f.close()

def XRBtoGROND(file, i, T_Disc, T_hs, Angle_hs,M_p):

    infile = open(file,'r').readlines()
    outfile = open(file+"mag",'w')

    outfile.write('Masse: %s    Neiungswinkel: %s    DiscTemperatur: %s    HSWinkel: %s    HSTemperatur: %s  \n'
                  %(M_p, i, T_Disc, Angle_hs, T_hs))

    outfile.write("# Phase \t g \t r \t i \t z \t J \t H \t K  \n")
    gAB = np.array([])
    rAB = np.array([])
    iAB = np.array([])
    zAB = np.array([])
    JAB = np.array([])
    HAB = np.array([])
    KAB = np.array([])
    phaseS = np.array([])
    
    LambdaMean = np.array([4586,6220,7641,8989,12398,16468,21706])
    Filtarea = np.ones(7)

    #print("Processing "+ file)
    for line in infile:
        try:
            K = FluxToMag( float( line.split()[26]),Filtarea[6],LambdaMean[6])
            Kratio = float( line.split()[27])/float( line.split()[28])
            H = FluxToMag( float( line.split()[22]),Filtarea[5],LambdaMean[5])
            Hratio = float( line.split()[23])/float( line.split()[24])
            J = FluxToMag( float( line.split()[18]),Filtarea[4],LambdaMean[4]) 
            Jratio = float( line.split()[19])/float( line.split()[20])
            z = FluxToMag( float( line.split()[14]),Filtarea[3],LambdaMean[3])
            zratio = float( line.split()[15])/float( line.split()[16])
            i = FluxToMag( float( line.split()[10]),Filtarea[2],LambdaMean[2]) 
            iratio = float( line.split()[11])/float( line.split()[12])
            r = FluxToMag( float( line.split()[6]),Filtarea[1],LambdaMean[1]) 
            rratio = float( line.split()[7])/float( line.split()[8])
            g = FluxToMag( float( line.split()[2]),Filtarea[0],LambdaMean[0]) 
            gratio = float( line.split()[3])/float( line.split()[4])
            phase = float( line.split()[0])
            
            phaseS = np.append(phaseS,phase)
            gAB = np.append(gAB,g)
            rAB = np.append(rAB,r)
            iAB = np.append(iAB,i)
            zAB = np.append(zAB,z)
            JAB = np.append(JAB,J)
            HAB = np.append(HAB,H)
            KAB = np.append(KAB,K)
            outfile.write(" %f  %f  %f  %f  %f  %f %f %f \n"% (phase, g, r, i, z, J, H, K))
            #return 
        except IndexError as err:
            pass
    outfile.close()
    Fmodnall = [phaseS]
    for filt in [gAB,rAB,iAB,zAB,JAB,HAB,KAB]:
        if np.size(filt)>1:
            #Fmodn = Normalize(filt)
            Fmodn = filt
            Fmodnall.append(Fmodn)
    return Fmodnall




    
    
def ExecXRB(dict, output, i, T_Disc, T_hs, Angle_hs,Dscale):
    parfile = '%s_tmp'%output
    LC_xrb_file  = '%s_tmp.LC'%output
    #LC_SysParFile ='%s_tmp.SysPa'%(output,)
    LC_GROND_file= '%s_tmp.LCmag'%output

    actual_dir = os.getcwd()
    workdir = '/home/simonste/sw/XRbinaryParalell/XRbpar1/'
    WriteParfile(dict, workdir+parfile)
    os.chdir(workdir)
    print '############### XRB output begin ############'
    os.system('%sa.out %s '%(workdir, parfile))
    print '############## XRB output end ##############'
    XRBtoGROND(LC_xrb_file, i, T_Disc, T_hs, Angle_hs,Dscale)

    os.chdir(actual_dir)
    os.system('cp %s ./%s'%(workdir+LC_GROND_file, output))
    os.system('cp %s_tmp.Sys* ./%s'%(workdir+output, output+'.SysPars'))
    return 


#### Physical stuff ##########################################

Ms = 1.9891e30 # solar mass in Kg
Rs = 9.955e8   # solar radius in m
Ts = 5778.0    # sun surface Temperature in Kelvin

G = const.G # gravitational constant in (m^3)/(Kg s^2)


### Secondary Mass
def M2(T):
    M2 = (T/Ts)**(4/2.5)
    #M2 = 0.26
    return M2
#%%
    
    
###Primary Mass 
### Radius Accretion disc Shakura & Sunyaev 1973
def rd(M_star, M_bh, P):
    q = M_star/M_bh
    return (0.6* a(M_star, M_bh, P))/(1+q)
#circulation radius king p.60 minimum outer radius
def R_circ(M_star, M_bh):
    q = M_star/M_bh
    return((1+q)*(0.5-0.227*log10(q))**4)
    
def R_L1(M_star, M_bh):
    q = M_star/M_bh
    
### Distance of objects
def a(M_star, M_bh, P):
    return (  ( (G*(M_star+M_bh)*Ms*(P*24*60*60)**2) /(4*pi**2) )**(1.0/3.0)  ) /Rs

### tidal disruption radius King
def rd_t(M_star, M_bh, P):
    q = M_star/M_bh
    return (0.9*0.49*q**(-2/3)*a(M_star, M_bh, P)/(0.6*q**(-2/3)+log(1+q**(-1/3))))
#%%
def R_in(M_bh):
    return(3*2*M_bh*Ms*G/const.c**2)

def ProduceLC(period, T_s,i, DiskOnOff,DiskrimOnOff, T_Disc, T_hs, Angle_hs,phaseoffset,Dscale, DiskSpotOnOff, HSwidth, Rinner,file ):
    workdir = '/home/simonste/sw/XRbinaryParalell/XRbpar1/'
    parameters =   ParseDefault(workdir+'default')

    T_pow = -0.75

    M_s = M2(T_s)
    #Mass Inklination relation
    Mlink = Massfunc()
    Mlink.f = fM(550)
    Mlink.f.fix = True
    Mlink.Mass_2 = M_s
    Mlink.Mass_2.fix = True
    M_p = Mlink(i)
    
    
    #consistency calculations
    
    q = M_s/M_p        

    #a(M_s, M_p, period) #distance objects [solar radii]
    #rd(M_s, M_p, period) #disc radius [solar radii]
    
    DSC = Diskminfunc()
    DSC.M_star = M_s
    DSC.M_star.fix = True
    DSC.P = period
    DSC.P.fix =True
    D_min = DSC(M_p)
    if Dscale < D_min:
        Dscale = D_min
    else:
        Dscale = Dscale
        
    
    A_Disc =   rd_t(M_s, M_p, period)/a(M_s, M_p, period)*Dscale    
    
    
    A_minDisc = Rinner*R_in(M_p)/3/(a(M_s,M_p,period)*Rs)
    Hedge = 0.1*(A_Disc-A_minDisc) #maximum
    
    #T_DiscEdge = L_Disc**T_pow
    T_DiscEdge = T_Disc
    
    
    
    
    #DISKRIM PArameters
    #arimwidth = 0.09
    #torusazero = A_Disc - 0.05
    #diskrim1 = diskrimstart -1
    #diskrim2 = diskrimstart + diskrimwidth
    
    
    #diskspot
    #Angle_hs2 = Angle_hs + HSwidth
    #Tratio = T_hs/ T_Disc
    #HSamin = A_Disc - HSawidth*(A_Disc-A_minDisc)
    
    
    
    
    #create parameter file
    parameters['INNERDISK']= 'OFF'
    parameters['DISK']= DiskOnOff
    parameters['DISKSPOTS']= DiskSpotOnOff
    parameters['DISKTORUS']= DiskrimOnOff
    #parameters['DISKTORUSAZERO']= '%s'%torusazero
    #parameters['DISKTORUSAWIDTH']= '%s'%arimwidth
    #parameters['DISKTORUSPARS']= 'POINT %s 0 %s \n DISKTORUSPARS= POINT %s %s %s \n DISKTORUSPARS= POINT %s 0 %s'%(diskrim1,T_DiscEdge,diskrimstart, diskrimheight,T_Rim,diskrim2,T_DiscEdge)
    #parameters['DISKSPOT']= '%s  %s %s %s %s' %(Angle_hs, Angle_hs2, HSamin,A_Disc, Tratio )
    parameters['DISKEDGET']= '%s  %s %s  %s'%(T_DiscEdge, T_hs, Angle_hs, HSwidth)
    parameters['MAINDISKT']= '%s  %s'%(T_Disc,T_pow)
    parameters['MAINDISKA']= '%s %s'%(A_minDisc,A_Disc)

    parameters['INCLINATION']= '%s'%i
    parameters['M1'] = '%s'%M_p
    parameters['STAR2TEMP']= '%s'%T_s
    parameters['IRRADIATION']= 'ON'
    parameters['MASSRATIO']= '%s'%q
    
    parameters['MAINDISKH']= '%s 1.125'%(Hedge)  #r**9/8 from king
    parameters['PERIOD'] = '%s'%period
    parameters['PHASEOFFSET']= '%s'%phaseoffset
    

    ExecXRB(parameters, file, i, T_Disc, T_hs, Angle_hs,Dscale)                    
########################################## Main Programm #################################################


## Specify
period = 2*0.1960797
phaseoffset = 0.5
DiskrimOnOff = 'OFF'
DiskOnOff = 'ON'
HsOnOff   = 'ON'
DiskSpotOnOff = 'OFF'  


#diskrimwidth_range = [40,80,120]
#diskrimstart_range = [45,225]
#diskrimheight_range = [0.05,0.1]
#T_Rim_range = [3000,5000,7000]
#HSawidthrange = [0.001,0.005,0.009]#,0.007,0.01,0.02,0.03,0.05,0.07,0.1]  #RATIO OF DISK
HSwidthrange = range(1,42,10)
NWrange=[60]#range(75,81,5)
Trange =range(1000,3200,300)
anglerange= [90]
THSrange=[3100]
Dscalerange= [0.725,0.825,0.925,0.999]
T_srange =[3800,4000,4200,4400,4600]
Rinner_range = [1000,3000,5000,7000,9000]

numofsim = len(NWrange)*len(Trange)*len(anglerange)*len(THSrange)*len(Dscalerange)*len(HSwidthrange)*len(T_srange)*len(Rinner_range)
count = 0

if DiskOnOff == 'ON':
    for i in NWrange:
        for T_Disc in Trange:
            for Dscale in Dscalerange:                
                for Angle_hs in anglerange:
                    if Angle_hs < 10:

                        AH = str(00)+str(Angle_hs)
                    elif Angle_hs < 100:
                        AH = str(0)+str(Angle_hs)
                    else:
                        AH = Angle_hs
                    for T_hs in THSrange:
                        TH = T_hs/10
                            
                        for HSwidth in HSwidthrange:
                            if HSwidth <10:
                                HSW = str(0)+str(HSwidth)
                            else: 
                                HSW = HSwidth                                         
                            for T_s in T_srange:
                                for R_inner in Rinner_range:
                                    file = '%s%s%s%s%s%s%s%s'%( i,T_Disc, AH, TH, int(Dscale*1000),HSW,T_s,R_inner)
                                    time1 = time.time()

                                    ProduceLC(period, T_s,i, DiskOnOff,DiskrimOnOff, T_Disc, T_hs, Angle_hs, phaseoffset, Dscale,DiskSpotOnOff, HSwidth,R_inner,file)

                                    os.system('mv %s ~/research/xrb/SimLC/HSRinner'%(file))
                                    os.system('mv %s ~/research/xrb/SimLC/HSRinner/SysParFiles'%(file+'.SysPars'))
                                    count += 1
                                    time2 = time.time()
                                    timediff = time2 -time1
                                    timeleft = timediff*(numofsim-count)/60
                                    print 'simulation ' + str(count) + ' of ' + str(numofsim) + ' complete!'                                
                                    print 'approx ' +str(timeleft) + ' minutes left'

                                                                                          


                                        

os.system('rm /home/simonste/sw/XRbinaryParalell/XRbpar1/*_tmp*')
#%%

