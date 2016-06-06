//Example problem for ERESS 2016
//Agenda: Add additional force (rebx and custom, r->usleep=0), add second planet, add 500 planetesimals, push to github.

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);
double E0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    //Parameter List
    int add_gr = 0;
    int N_planets = 2;
    int N_planetesimals = 500;
    
	//Setup - Simulation
	r->integrator	= REB_INTEGRATOR_HYBARID;
    r->heartbeat	= heartbeat;
    r->ri_hybarid.switch_ratio = 3;
    r->ri_hybarid.CE_radius = 20.;
    r->testparticle_type = 1;
    r->dt = 0.01;
    
    //Setup - Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->collisions_track_dE = 1;
    
    //Setup - Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 6;
    reb_configure_box(r,boxsize,2,2,1);
    
    srand(14);

    //Add Particles
    //Star
    struct reb_particle star = {0};
    star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
    reb_add(r, star);
    
    double m=5e-4, e=0.5, inc=reb_random_normal(0.00001), rp=0.00042;
    {//Planet 1
        double a=1;
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p.r = rp;       //radius of planet (AU)
        p.id = r->N;
        reb_add(r, p);
    }
    double a_outer = 3;
    if(N_planets>1){//Planet 2
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a_outer, 0, inc, 0, 0, 0);
        p.r = rp;       //radius of planet (AU)
        p.id = r->N;
        reb_add(r, p);
    }
    r->N_active = r->N;
    
    //Add Planetesimal Disk
    double amin = a_outer, amax = a_outer+2;
    {
        double planetesimal_mass = 1e-9;
        double powerlaw = 0;
        while(r->N<N_planetesimals + r->N_active){
            r->usleep = 1000;
            struct reb_particle pt = {0};
            double a    = reb_random_powerlaw(amin,amax,powerlaw);
            double e    = reb_random_rayleigh(5e-3);   //rayleigh dist
            double inc  = reb_random_rayleigh(5e-3);
            double Omega = reb_random_uniform(0,2.*M_PI);
            double apsis = reb_random_uniform(0,2.*M_PI);
            double phi 	= reb_random_uniform(0,2.*M_PI);
            pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
            pt.r 		= 0.00000934532;
            pt.id = r->N;
            reb_add(r, pt);
        }
    }

    //Final setup parameters
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    system("rm -f energy.txt");
    
    //Set GR Potential
    if(add_gr){
        struct rebx_extras* rebx = rebx_init(r);
        double c = 173.26203208556151;
        int source_index = 0;
        rebx_add_gr_potential(rebx, source_index, c);
    }
    
    //Integrate!
    reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    
    //Example force
    //r->particles[1].vx += 1e-4;
    //r->usleep = 1000;
    
    //output energy
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0+r->collisions_dE)/E0);
        reb_output_timing(r, 0);
        printf("dE = %e",relE);
    }
}