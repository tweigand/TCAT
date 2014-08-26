from math import *
import numpy
import proteus
from Dilute_parameters import *
from proteus.TransportCoefficients import *

class Dilute_Transport(TC_base):
    """
    (Mu)_t + \deld (Bu - A \grad phi) + C u = 0 
    """
    
    def __init__(self,nc=2,nd=1):
        mass={}
        advection={}
        diffusion={}
        potential={}
        reaction={}
        hamiltonian={}
        self.nd=nd
        variableNames=['mass_frac','activity']
        mass = {0:{0:'nonlinear'},
                1:{1:'constant'}}
        advection = {0:{0:'nonlinear'},
                     1:{1:'constant'}}
        diffusion = {0:{0:{0:'nonlinear'},
                        1:{0:'nonlinear',
                           1:'constant'}},
                     1:{0:{0:'constant'},
                        1:{1:'constant'}}}
        potential = {0:{0:'u'},
                     1:{1:'u'}}
        reaction = {0:{0:'constant'},
                    1:{0:'nonlinear',
                       1:'linear'}}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian)

    
    def evaluate(self,t,c):

        space_dim = c[('f',0)].shape[-1]
        grad_mass_frac_all = numpy.reshape(c['grad(u)',0],(len(c[('u',0)].flat),space_dim))

        for k in range(len(c[('u',0)].flat)):
            
            mass_frac = c[('u',0)].flat[k]
            dummy_act = c[('u',1)].flat[k]
            grad_mass_frac = grad_mass_frac_all.flat[k]

            rho = density(mass_frac)
            d_rho = d_density(mass_frac)

            mole_frac = mol_frac(mass_frac)
            d_mole_frac = d_mol_frac(mass_frac)
            
            MW_w = MW_fluid(mass_frac)
            d_MW_w = d_MW_fluid(mass_frac)

            D_abwe = D_eff(mass_frac,grad_mass_frac)
            d_D_abwe = d_D_eff(mass_frac,grad_mass_frac)

            act = Activity(mass_frac)
            d_act = d_Activity(mass_frac)

            vol = Vol(mass_frac)
            d_vol = d_Vol(mass_frac)

## Transport ##

            denom = visc(mass_frac)-c1*perm*grad_density(mass_frac,grad_mass_frac)
            v_advec = perm/denom*grav*rho*rho
            d_v_advec = perm*(d_visc(mass_frac)-c1*perm*d_grad_density(mass_frac,grad_mass_frac))/(denom*denom)*rho*rho*grav + 2.0*perm/denom*rho*d_rho*grav       

            print perm/denom*grav*rho

            c[('m',0)].flat[k] = vol_frac*density(mass_frac)*mass_frac
            c[('dm',0,0)].flat[k] = vol_frac*(d_density(mass_frac)*mass_frac + density(mass_frac))

            
            c[('f',0)].flat[k] = v_advec*mass_frac
            c[('df',0,0)].flat[k] = d_v_advec*mass_frac + v_advec          

            for j in range(space_dim):
              
                c[('a',0,0)].flat[k*space_dim**2+j*space_dim+j] = vol_frac*density(mass_frac)*D_e
                c[('da',0,0,0)].flat[k*space_dim**2+j*space_dim+j] = vol_frac*D*d_density(mass_frac)

                term1 = (1.0-mole_frac)/(1.0-mass_frac)
                d_term1 = ( (mass_frac-1.0)*d_mole_frac-mole_frac+1.0 ) / ( (mass_frac-1.0)*(mass_frac-1.0) )
                
                term2 = MW_w*MW_w/(MW_a*MW_a*MW_b)
                d_term2 = 2.0*MW_w*d_MW_w/(MW_a*MW_a*MW_b)

                c[('a',0,0)].flat[k*space_dim**2+j*space_dim+j] = 0.0#vol_frac*R*temp*term1*term2*rho*D_abwe
                c[('da',0,0,0)].flat[k*space_dim**2+j*space_dim+j] = 0.0#vol_frac*R*temp*(d_term1*term2*rho*D_abwe + term1*d_term2*rho*D_abwe + term1*term2*d_rho*D_abwe + term1*d_term2*rho*d_D_abwe)

                term2 = MW_w/(MW_a*MW_a)
                d_term2 = d_MW_w/(MW_a*MW_a)

                term3 = 1.0/act
                d_term3 = -d_act/(act*act)

                c[('a',0,1)].flat[k*space_dim**2+j*space_dim+j] = 0.0#vol_frac*R*temp*term1*term2*term3*rho*D_abwe*mass_frac
                c[('da',0,1,0)].flat[k*space_dim**2+j*space_dim+j] = 0.0#vol_frac*R*temp*(d_term1*term2*term3*rho*D_abwe*mass_frac + term1*d_term2*term3*rho*D_abwe*mass_frac + term1*term2*d_term3*rho*D_abwe*mass_frac + term1*term2*term3*d_rho*D_abwe*mass_frac + term1*term2*term3*rho*d_D_abwe*mass_frac + term1*term2*term3*rho*D_abwe)


### Dummy Equation for Activity ##

            c[('r',1)].flat[k] = dummy_act - act
            c[('dr',1,0)].flat[k] = -d_act
            c[('dr',1,1)].flat[k] = 1.0       

