#This is used to estimate what the timestep should be for HYBARID given the initial conditions.

#Global params
M = 1
HSR = 1
HSR_frac = 0.1

#planet 1
m1 = 1e-3
a1 = 2.515
e1 = 0
r1_hill = a1*(m1/3.)**(1./3.)

#planet 2
m2 = 1e-3
a2 = a1 + r1_hill
e2 = 0.5

#calc timestep
v_rel = ((M/a2)*(1+e2)/(1-e2))**0.5 - ((M/a1)*(1+e1)/(1-e1))**0.5
dr = HSR_frac*HSR*r1_hill
dt = dr/v_rel
print dt/4.