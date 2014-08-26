from math import *
import numpy
import proteus
from Classic_Transport_parameters import *
from proteus.TransportCoefficients import *

class Classic_Transport(TC_base):
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
        mass = {0:{0:'nonlinear',
                   1:'constant'},
                1:{1:'nonlinear'}}
        advection = {0:{0:'constant',
                        1:'nonlinear'},
                     1:{1:'nonlinear'}}
        diffusion = {0:{0:{0:'constant',
                           1:'nonlinear'},
                        1:{0:'constant',
                           1:'constant'}},
                     1:{0:{0:'constant',
                           1:'nonlinear'},
                        1:{0:'constant',
                           1:'nonlinear'}}}
        potential = {0:{0:'u'},
                     1:{1:'u'}}
        reaction = {0:{0:'constant'},
                    1:{1:'constant'}}
        TC_base.__init__(self,
                         nc,
                         mass,
                         advection,
                         diffusion,
                         potential,
                         reaction,
                         hamiltonian)
        self.variableNames=['pressure','mass_frac']
    
    def evaluate(self,t,c):

        space_dim = c[('f',0)].shape[-1]
        grad_press_all = numpy.reshape(c['grad(u)',0],(len(c[('u',0)].flat),space_dim))
        for k in range(len(c[('u',0)].flat)):

            grad_press = grad_press_all.flat[k]            

            mass_frac = c[('u',1)].flat[k]

            rho = density(mass_frac)
            d_rho = d_density(mass_frac)
  
            mu = visc(mass_frac)
            d_mu = d_visc(mass_frac)

## Flow ##

            c[('m',0)].flat[k] = vol_frac*rho
            c[('dm',0,1)].flat[k] = vol_frac*d_rho

            c[('f',0)].flat[k] = perm/mu*grav*rho*rho
            c[('df',0,1)].flat[k] = perm*grav*( -d_mu/(mu*mu)*rho*rho + 2.0*mu*rho*d_rho)

            for j in range(space_dim):
              
                c[('a',0,0)].flat[k*space_dim**2+j*space_dim+j] = perm/mu*rho
                c[('da',0,0,1)].flat[k*space_dim**2+j*space_dim+j] = perm*( -d_mu/(mu*mu)*rho + d_rho/mu)

## Transport ##

            c[('m',1)].flat[k] = vol_frac*rho*mass_frac
            c[('dm',1,1)].flat[k] = vol_frac*(d_rho*mass_frac + rho)

            c[('f',1)].flat[k] = perm/mu*grav*rho*rho*mass_frac
            c[('df',1,1)].flat[k] = perm*grav*( -d_mu/(mu*mu)*rho*rho*mass_frac + 2.0/mu*rho*d_rho*mass_frac + rho*rho/mu )

            for j in range(space_dim):
              
                c[('a',1,0)].flat[k*space_dim**2+j*space_dim+j] = perm/mu*rho*mass_frac
                c[('da',1,0,1)].flat[k*space_dim**2+j*space_dim+j] = perm*( -d_mu/(mu*mu)*rho*mass_frac + d_rho/mu*mass_frac + rho/mu)

                c[('a',1,1)].flat[k*space_dim**2+j*space_dim+j] = vol_frac*rho*D_e
                c[('da',1,1,1)].flat[k*space_dim**2+j*space_dim+j] = vol_frac*d_rho*D_e