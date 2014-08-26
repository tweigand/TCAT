from proteus import *
from proteus.default_p import *
from Classic_Transport_parameters import *
from Classic_Transport import *

name="Classic_Transport"

L=(88.9,1.0,1.0)

nc = 2
nd = 1
T = 2.0e4

mass_frac_in = 0.01254

coefficients = Classic_Transport()

### Initial conditions ###
class Mass_Fraction_IC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        if (x[0] <= 0.0):
          return mass_frac_in
        else:
          return 0.0

class Pressure_IC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        if (x[0] == 0.0):
          return density(mass_frac_in)*grav*(x[0]-L[0])
        else:
          return density(0.0)*grav*(x[0]-L[0])


initialConditions = {0:Pressure_IC(mass_frac_in),1:Mass_Fraction_IC(mass_frac_in)}

def getDBC_mass_frac(x,flag):
    if (x[0]==0.0):
    	return lambda x,t: mass_frac_in
    elif (x[0]==L[0]):
        return lambda x,t: 0.0

def getDBC_pressure(x,flag):
    if (x[0]==L[0]):
    	return lambda x,t: 0.0

def mass_fraction_noflux(x,flag):
    if (x[0]==L[0]):
         return lambda x,t: 0.0

def pressure_flux(x,flag):
    if (x[0]==0):
         return lambda x,t: -Q/A*density(mass_frac_in)

def no_flux(x,flag):
    return lambda x,t: 0.0 


dirichletConditions = {0:getDBC_pressure,1:getDBC_mass_frac}
fluxBoundaryConditions = {0:'setFlow'}
advectiveFluxBoundaryConditions =  {0:no_flux,1:mass_fraction_noflux}
diffusiveFluxBoundaryConditions = {0:{0:pressure_flux},1:{1:mass_fraction_noflux}}