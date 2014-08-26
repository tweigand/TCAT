import numpy


# Units are cm, grams, sec


vol_frac = 0.3341
perm = 5.04e-6   #cm/sec

grav = -980.

Q = 9.469e-3 # Pump Rate into Column
A = 4.9087   # Cross-Sectional Area of Column

R = 8.3145e7
temp = 298.
D_e = 8.36e-4

def density(m_f):
  a0 = 0.9982
  a1 = 0.8411
  a2 = 0.5038
  a3 = 0.8250
  return a0 + a1*m_f + a2*m_f*m_f + a3*m_f*m_f*m_f

def d_density(m_f):
  if (m_f < 0.0):
    return 0.0
  else:
    a0 = 0.9982
    a1 = 0.8411
    a2 = 0.5038
    a3 = 0.8250
    return a1 + 2.0*a2*m_f + 3.0*a3*m_f*m_f

def grad_density(m_f,grad_m_f):
  a0 = 0.9982
  a1 = 0.8411
  a2 = 0.5038
  a3 = 0.8250
  return a1*grad_m_f + 2.0*a2*m_f*grad_m_f + 3.0*a3*m_f*m_f*grad_m_f

def d_grad_density(m_f,grad_m_f):
  a0 = 0.9982
  a1 = 0.8411
  a2 = 0.5038
  a3 = 0.8250
  return 2.0*a2*grad_m_f + 6.0*a3*m_f*grad_m_f

def visc(m_f):

  a0 = -0.03221 
  a1 = 1.861
  a2 = -6.262
  a3 = 17.18
  conv = 0.01
  return conv * numpy.exp( a0 + a1*m_f + a2*m_f*m_f + a3*m_f*m_f*m_f )

def d_visc(m_f):
  if (m_f < 0.0):
    return 0.0
  else:
    a0 = -0.03221 
    a1 = 1.861
    a2 = -6.262
    a3 = 17.18
    conv = 0.01
    return conv * (a1 + 2.0*a2*m_f + 3.0*a3*m_f*m_f) * numpy.exp( a0 + a1*m_f + a2*m_f*m_f + a3*m_f*m_f*m_f )
