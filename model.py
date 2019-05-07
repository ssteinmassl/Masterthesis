def model(phase,Inkl,Dtemp,HSTemp,HSWidth,HSAwidth,T_sec,distfac,column):
    #Dsize,HSAngle,penalty
    path = '/home/simonste/research/xrb/SimLC/HSOnDiskSpotOn/'
    
    
    Inklcode = int(round(Inkl))
    if Dtemp < 1000:
        Dtempcode = '0' + str(int(round(Dtemp)))
    else:
        Dtempcode = int(round(Dtemp)) 
    #if HSAngle <100:
    #    AHScode = '0' + str(int(round(HSAngle)))
    #else:
    #   AHScode = int(round(HSAngle))
    AHScode = '090'
        
    THScode = int(round(HSTemp / 10.))
    
    #Dsizecode = int(round(Dsize*1000))
    Dsizecode ='800'
    
    if HSWidth <10:
        WHScode = '0' + str(int(round(HSWidth)))
    else:
        WHScode = int(round(HSWidth))
        
    AWHScode = int(round(HSAwidth*1000))
    T_seccode = int(round(T_sec))
    
    fil = '%s%s%s%s%s%s%s%s' %( Inklcode, Dtempcode,AHScode,THScode, Dsizecode,WHScode,AWHScode,T_seccode)
    ph, mag = np.loadtxt(path + fil, skiprows = 2, usecols = (0,column+1), unpack = True)
    #distfac is additive normalization parameter
    magd = mag+distfac
    ###Penalty for Chceck if Dsiye physical later in linking function
    #mag = penalty*mag
    #return mag    
    return magd[np.where(ph == np.round(phase,4))]    