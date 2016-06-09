#This is used to estimate what the timestep should be for HERMES given the initial conditions.

#Global params
M = 1
HSR = 1
HSR_frac = 0.1
safe_factor = 4

#planet 1
m1 = 6e-6
a1 = 5.0
e1 = 0
r1_hill = a1*(m1/3.)**(1./3.)

#planet 2
m2 = 3e-8
a2 = a1 + r1_hill
e2 = 0.2

#calc timestep
v_rel = ((M/a2)*(1+e2)/(1-e2))**0.5 - ((M/a1)*(1+e1)/(1-e1))**0.5
dr = HSR_frac*HSR*r1_hill
dt = dr/v_rel
print dt/safe_factor
