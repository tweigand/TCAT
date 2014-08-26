import numpy


# Units are cm, grams, sec


vol_frac = 0.3341
perm = 5.04e-6   #cm/sec
c1 = 5.2689e-2
c2 = 5.3906e11
grav = 980.

R = 8.3145e7
temp = 298.
D_e = 8.36e-4
D = 7.48e-11

MW_a = 199.8 # g/mol
MW_b = 18.015 # g/mol

def density(m_f):
  a0 = 0.9982
  a1 = 0.8411
  a2 = 0.5038
  a3 = 0.8250
  return a0 + a1*m_f + a2*m_f*m_f + a3*m_f*m_f*m_f

def d_density(m_f):
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
  a0 = -0.03221 
  a1 = 1.861
  a2 = -6.262
  a3 = 17.18
  conv = 0.01
  return conv * (a1 + 2.0*a2*m_f + 3.0*a3*m_f*m_f) * numpy.exp( a0 + a1*m_f + a2*m_f*m_f + a3*m_f*m_f*m_f )

def mol_frac(m_f):
  return MW_a*m_f/(MW_a*m_f+MW_b*(1.0-m_f))

def d_mol_frac(m_f):
  return (MW_a*MW_b)/pow(MW_b*(m_f-1.0)-MW_a*m_f,2.0)

def MW_fluid(m_f):
  x = mol_frac(m_f)
  return MW_a*x + MW_b*(1.0-x)

def d_MW_fluid(m_f):
  dx = d_mol_frac(m_f)
  return MW_a*dx - MW_b*dx

def Vol(m_f):
  a0 = 0.9982
  a1 = 0.8411
  a2 = 0.5038
  a3 = 0.8250
  den = density(m_f)
  d_p_w =  -(a1 + 2.0*a2*m_f + 3.0*a3*m_f*m_f)/(den*den)
  return 1.0/den + (1.0 - m_f)*d_p_w

def d_Vol(m_f):
  a0 = 0.9982
  a1 = 0.8411
  a2 = 0.5038
  a3 = 0.8250
  den = density(m_f)
  term1 = d_density(m_f)
  d_p_w =  -term1/(den*den)
  d_p_w_2 = ( -(2.0*a2 + 6.0*a3*m_f)*den + 2.0*term1*term1 ) / (den*den*den)
  return d_p_w + d_p_w_2 - m_f*d_p_w_2 - d_p_w

def Molality(m_f):
  return 1000.*m_f/( (1-m_f)*MW_a)

def d_Molality(m_f):
  return 1000./( MW_a*(m_f - 1.0)*(m_f - 1.0) )

def Ionic(m_f):
  return 3.0*Molality(m_f)

def d_Ionic(m_f):
  dm = d_Molality(m_f)
  return 3.0*dm

def Activity(m_f):
  b0 = 2.3525
  b1 = 1.793589
  b2 = 3.244255e-1
  b3 = 2.086318e-1
  b4 = -5.6619089e-2
  b5 = 1.21564998e-2
  b6 = -1.2930306e-3
  b7 = 4.84942655e-5 
  m = Molality(m_f)
  I = Ionic(m_f)

  if (m_f < 1.e-12):
    return 1.0
  else:
    return numpy.exp( -b0*pow(I,0.5)/(1.+b1*pow(I,0.5)) + b2*m + b3*m*m + b4*m*m*m + b5*m*m*m*m + b6*m*m*m*m*m + b7*m*m*m*m*m*m )

def d_Activity(m_f):
  b0 = 2.3525
  b1 = 1.793589
  b2 = 3.244255e-1
  b3 = 2.086318e-1
  b4 = -5.6619089e-2
  b5 = 1.21564998e-2
  b6 = -1.2930306e-3
  b7 = 4.84942655e-5 
  
  m = Molality(m_f)
  dm = d_Molality(m_f)
  
  I = Ionic(m_f)
  dI = d_Ionic(m_f)
  if (m_f < 1.e-12):
    return 0.0
  else:
    arg = numpy.exp( -b0*pow(I,0.5)/(1.+b1*pow(I,0.5)) + b2*m + b3*m*m + b4*m*m*m + b5*m*m*m*m + b6*m*m*m*m*m + b7*m*m*m*m*m*m )
    deriv = -b0*dI/(2.0*pow(I,0.5)*(b1*pow(I,0.5) + 1.0)*(b1*pow(I,0.5) + 1.0) ) + b2*dm + 2.0*b3*m*dm + 3.0*b4*m*m*dm + 4.0*b5*m*m*m*dm + 5.0*b6*m*m*m*m*dm + 6.0*b7*m*m*m*m*m*dm 
    return arg*deriv

def Diff(m_f):
  MW_w = MW_fluid(m_f) 
  return D*MW_a*MW_a/(R*temp*MW_w)

def d_Diff(m_f):
  MW_w = MW_fluid(m_f)
  d_MW_w = d_MW_fluid(m_f) 
  return -D*MW_a*MW_a/(R*temp)*d_MW_w/(MW_w*MW_w)


def D_eff(m_f,grad_m_f):
  D_w = D
  den = density(m_f)
  x = Molality(m_f)
  MW_w = MW_fluid(m_f) 
  return D_w/(1.0-(vol_frac*c2*grad_m_f/den)*(1.0-x)/(1.0-m_f)*(MW_w/MW_a)*D_w)

def d_D_eff(m_f,grad_m_f):
  D_w = D
  d_D_w = 0.0
  d_grad_m_f = d_grad_density(m_f,grad_m_f)
  den = density(m_f)
  d_den = d_density(m_f)
  x = Molality(m_f)
  d_x = d_Molality(m_f)
  MW_w = MW_fluid(m_f) 
  d_MW_w = d_MW_fluid(m_f)
  
  h = 1.0-(vol_frac*c2*grad_m_f/den)*(1.0-x)/(1.0-m_f)*(MW_w/MW_a)*D_w  
  
  term1 = vol_frac*c2*grad_m_f/den
  d_term1 = vol_frac*c2*(den*d_grad_m_f - d_den*grad_m_f)/(den*den)
  
  term2 = (1.0-x)/(1.0-m_f)
  d_term2 = ((m_f-1)*d_x-x+1)/( (m_f-1.0)*(m_f-1.0) )

  term3 = MW_w/MW_a
  d_term3 = d_MW_w/MW_a

  return (d_D_w*h - D_w*(d_term1*term2*term3*D_w + term1*d_term2*term3*D_w + term1*term2*d_term3*D_w + term1*term2*term3*d_D_w) )/(h*h)


