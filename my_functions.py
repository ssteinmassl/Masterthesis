from threeML import *
from math import *
from math import *
from scipy import constants as const
from sympy.solvers import solve
from sympy import Symbol, pprint, S

Rs = 9.955e8 #[m]  
Ts = 5778.  
G = const.G
Ms = 1.9891e30

def fM(K):
    Ms = 1.9891e30
    P= 2 * 0.1960797
    f = K**3 * 10**9 * P *24*3600/ (2*const.pi*const.G)/Ms
    return f

def M2(T):
    M2 = (T/5778.)**(4/2.5)
    #M2 = 0.26
    return M2


def M1(ink):
    x = Symbol('x')
    f = fM(550)
    
    sol = solve(f*(1+M2(4200.)/x)**2 / sin(ink/180. *const.pi)**3 - x, x)[0]
    return np.array(sol)

def R_circ(M_star, M_bh):
    q = M_star/M_bh
    return((1+q)*(0.5-0.227*log10(q))**4)

### Distance of objects [solar radii]
def adistance(M_star, M_bh, P):
    adist = np.array([(  ( (G*(M_star+M_bh)*Ms*(P*24*60*60)**2) /(4*const.pi**2) )**(1.0/3.0)  ) /Rs])
    return adist

### tidal disruption radius King [solar radii]
def rd_t(M_star, M_bh, P):
    q = M_star/M_bh
    tidalradius = 0.9*0.49*q**(-2/3)*adistance(M_star, M_bh, P)/(0.6*q**(-2/3)+log(1+q**(-1/3)))
    return tidalradius
    #return 
    
def Radiusrange(M_bh):
    M_star = M2(4200.)
    P = 2 * 0.1960797
    Rratio = R_circ(M_star, M_bh)/(rd_t(M_star, M_bh, P)/adistance(M_star, M_bh, P))
    return Rratio

class Massfunc(Function1D):
        r"""
        description :
        
        parameters :
            f: 
                desc : func value
                initial value : 1
                
            Mass_2: 
                desc : sec Mass
                initial value : 1
        """

        __metaclass__ = FunctionMeta
           
        def _set_units(self,x_unit, y_unit):
            self.f.unit = astropy_units.dimensionless_unscaled
            self.Mass_2.unit = astropy_units.dimensionless_unscaled
        def evaluate(self,x,f,Mass_2):
            m = Symbol('m')
            xrad = x/180. *np.pi
            xx = Sin(f= 1.0 / (2 * np.pi), K = 1., phi =  0)
            sol = solve(f*(1+Mass_2/m)**2/xx(xrad)**3 -m,m)[0]
            #sol = solve(f*(1+Mass_2/m)**2 / sin(x/180. *const.pi)**3 - m, m)[0]
            return sol

    



class Diskminfunc(Function1D):
        r"""
        description : Function for Diskmin as ratio of r_circ to r_tidal
        
        parameters :
            M_star:
                desc :
                initial value : 0.6
                
            P:
                desc:
                initial value: 0.19
                
        """

        __metaclass__ = FunctionMeta
           
        def _set_units(self,x_unit, y_unit):
            self.M_star.unit = astropy_units.dimensionless_unscaled
            self.P.unit = astropy_units.dimensionless_unscaled
        def evaluate(self,x,M_star,P):
            Rratio = R_circ(M_star, x)/(rd_t(M_star, x, P)/adistance(M_star, x, P))
            return Rratio
        
        
def model(phase,Inkl,Dtemp,Dsize,HSTemp,HSWidth,T_sec,distfac,column):
    #Dsize,HSAngle,penalty,HSAwidth,
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
    
    Dsizecode = int(round(Dsize*1000))
    #Dsizecode ='800'
    
    if HSWidth <10:
        WHScode = '0' + str(int(round(HSWidth)))
    else:
        WHScode = int(round(HSWidth))
        
    #AWHScode = int(round(HSAwidth*1000))
    AWHScode = '5'
    T_seccode = int(round(T_sec))
    
    fil = '%s%s%s%s%s%s%s%s' %( Inklcode, Dtempcode,AHScode,THScode, Dsizecode,WHScode,AWHScode,T_seccode)
    ph, mag = np.loadtxt(path + fil, skiprows = 2, usecols = (0,column+1), unpack = True)
    #distfac is additive normalization parameter
    magd = mag+distfac
    ###Penalty for Chceck if Dsiye physical later in linking function
    #mag = penalty*mag
    #return mag    
    return magd[np.where(ph == np.round(phase,4))]   


def modelRim(phase,Inkl,Dtemp,HSTemp,HSWidth,HSAwidth,T_sec,distfac,Drim_width,Drim_start,Drim_height,Drim_Temp,column):
    #Dsize,HSAngle,penalty
    path = '/home/simonste/research/xrb/SimLC/HSRIM/'
    
    
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
    Drim_widthcode= int(round(Drim_width))
    Drim_startcode= int(round(Drim_start))
    Drim_heightcode= int(round(Drim_height*100))
    Drim_Tempcode= int(round(Drim_Temp))
    
    fil = '%s%s%s%s%s%s%s%s%s%s%s%s' %( Inklcode, Dtempcode,AHScode,THScode, Dsizecode,WHScode,AWHScode,T_seccode,Drim_widthcode,Drim_startcode,Drim_heightcode,Drim_Tempcode)
    ph, mag = np.loadtxt(path + fil, skiprows = 2, usecols = (0,column+1), unpack = True)
    #distfac is additive normalization parameter
    magd = mag+distfac
    ###Penalty for Chceck if Dsiye physical later in linking function
    #mag = penalty*mag
    #return mag    
    return magd[np.where(ph == np.round(phase,4))]

def modelRinner(phase,Inkl,Dtemp,HSTemp,HSWidth,T_sec,distfac,Dsize,R_in,column):
    
    path = '/home/simonste/research/xrb/SimLC/HSRinner/'
    
    
    Inklcode = int(round(Inkl))
    if Dtemp < 1000:
        Dtempcode = '0' + str(int(round(Dtemp)))
    else:
        Dtempcode = int(round(Dtemp)) 
    
    AHScode = '090'
        
    THScode = int(round(HSTemp / 10.))
    
    Dsizecode = int(round(Dsize*1000))
    
    if HSWidth <10:
        WHScode = '0' + str(int(round(HSWidth)))
    else:
        WHScode = int(round(HSWidth))
        
    
    T_seccode = int(round(T_sec))
    R_in_code = int(round(R_in))
    
    fil = '%s%s%s%s%s%s%s%s' %( Inklcode, Dtempcode,AHScode,THScode, Dsizecode,WHScode,T_seccode,R_in_code)
    ph, mag = np.loadtxt(path + fil, skiprows = 2, usecols = (0,column+1), unpack = True)
    #distfac is additive normalization parameter
    magd = mag+distfac
    
    return magd[np.where(ph == np.round(phase,4))]
