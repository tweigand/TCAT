from proteus import *
from proteus.default_p import *
from Dilute_parameters import *
from Dilute_Transport import *

name="Dilute_Transport"

L=(88.9,1.0,1.0)

nc = 2
nd = 1
T = 2.0e2

mass_frac_in = 1.e-6

coefficients = Dilute_Transport()

### Initial conditions ###
class Mass_Fraction_IC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        if (x[0] <= 0.0):
          return mass_frac_in
        else:
          return 0.0

class Activity_IC:
    def __init__(self,val=0.0):
        self.val = val
    def uOfXT(self,x,t):
        if (x[0] <= 0.0):
          return Activity(mass_frac_in)
        else:
          return Activity(0.0)

initialConditions = {0:Mass_Fraction_IC(mass_frac_in),1:Activity_IC(mass_frac_in)}

def getDBC_mass_frac(x,flag):
    if (x[0]==0.0):
    	return lambda x,t: mass_frac_in

def getDBC_activity(x,flag):
    if (x[0]==0.0):
    	return lambda x,t: Activity(mass_frac_in)


def mass_fraction_noflux(x,flag):
    if (x[0]==L[0]):
         return lambda x,t: 0.0

def activity_noflux(x,flag):
    if (x[0]==L[0]):
         return lambda x,t: 0.0


dirichletConditions = {0:getDBC_mass_frac,1:getDBC_activity}
#fluxBoundaryConditions = {0:'setFlow',1:'setFlow'}
advectiveFluxBoundaryConditions =  {0:mass_fraction_noflux,1:activity_noflux}
diffusiveFluxBoundaryConditions = {0:{0:mass_fraction_noflux},1:{1:activity_noflux}}