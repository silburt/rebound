import rebound
import numpy as np
import matplotlib.pyplot as plt
import sys

setup_checkpoint = 1
do_reverse = 0

if setup_checkpoint == 1:
    Np = 3
    tmax = 50
    times = np.linspace(tmax/Np,tmax,Np)
    print times
    prad = 2e-4                                 #planet radius
    mp = 1e-10
    
    sim = rebound.Simulation()
    sim.add(m=1.,r=0.0046)
    sim.add(m=1e-3,a=1.,r=0.5*prad)

    ev = 2*np.sqrt(2*sim.particles[1].m/prad)  #escape velocity
    
    for i in xrange(0,Np):
        phi = np.random.random()*np.pi*2
        x = 1.*sim.particles[1].x + prad*np.cos(phi)
        y = 1.*sim.particles[1].y + prad*np.sin(phi)
        vx = ev*np.cos(phi)
        vy = ev*np.sin(phi)
        sim.add(m=mp, x=x, y=y, vx=vx, vy=vy,r=1e-5)
        tsteps = np.linspace(sim.t,times[i],10)
        for t in tsteps:
            sim.integrate(t)
            plt.scatter(sim.particles[0].x,sim.particles[0].y,color='yellow',s=20)
            plt.scatter(sim.particles[1].x,sim.particles[1].y,color='red',s=8)
            for p in sim.particles[2:]:
                plt.scatter(p.x,p.y,color='blue',s=5)
        sys.stdout.write("\r" + str(sim.t))
        sys.stdout.flush()
    sim.save("ReverseCollisions.bin")
    del sim
    plt.show()
    plt.clf()

if do_reverse == 1:
    sim = rebound.Simulation.from_file("ReverseCollisions.bin")
    print(sim.particles[1].r)
    sim.integrator='hybarid'
    sim.collisions_track_dE = 1
    sim.collision = 'direct'
    sim.collision_resolve = 'merge'
    sim.testparticle_type = 1
    sim.dt = 0.0001
    sim.ri_hybarid.switch_radius = 6
    sim.ri_hybarid.CE_radius = 20
    sim.N_active = 2
    E0 = sim.calculate_energy()
    dE = np.zeros(0)
    times = np.logspace(np.log(sim.t-0.1),np.log(0.1),1000,base=np.e)
    sim.status()
    for t in times:
        sim.integrate(t)
        Ei = sim.calculate_energy() + sim.collisions_dE
        dE = np.append(dE,(Ei-E0)/E0)
        sys.stdout.write("\r" + str(sim.t))
        sys.stdout.flush()
    sim.status()

    plt.plot(times[dE>0],abs(dE[dE>0]),'o',ms=3, markeredgecolor='none',color='blue',label='positive values')
    plt.plot(times[dE<0],abs(dE[dE<0]),'o',ms=3, markeredgecolor='none',color='green',label='negative values')
    plt.xlim([max(times),0.1])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('time')
    plt.ylabel('signed dE/E')
    plt.legend(loc='lower right',prop={'size':10})
    plt.show()



