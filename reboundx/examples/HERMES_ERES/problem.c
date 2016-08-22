//Example problem for ERESS 2016
//Agenda: Add additional force (rebx and custom, r->usleep=0), add second planet, add 500 planetesimals, compare to WHFAST.

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);
double E0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    //Parameter List
    int N_planets = 1;
    int N_planetesimals = 0;
    
    //Set GR Potential
//    struct rebx_extras* rebx = rebx_init(r);
//    rebx_add_gr_potential(rebx, 0, 173.26203208556151);
    
	//Setup - Simulation
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;
    r->ri_hermes.solar_switch_factor = 20.;
    r->testparticle_type = 1;
    r->dt = 0.01;
    
    //Setup - Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->collision_resolve_keep_sorted = 1;
    r->track_energy_offset = 1;
    
    srand(14);

    //Star
    struct reb_particle star = {0};
    star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
    reb_add(r, star);
    
    //Planet 1
    struct reb_particle p = {0};
    p = reb_tools_orbit_to_particle(r->G, star, 5e-4, 1, 0.5, 0, 0, 0, 0);
    p.r = 0.00042;       //radius of planet (AU)
    reb_add(r, p);

    if(N_planets > 1){
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, 5e-4, 3, 0, reb_random_normal(0.00001), 0, 0, 0);
        p2.r = 0.00042;       //radius of planet (AU)
        reb_add(r, p2);
    }
    r->N_active = r->N;
    
    //Planetesimal Disk
    while(r->N<N_planetesimals + r->N_active){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(3,5,0);
        double e    = reb_random_rayleigh(5e-3);   //rayleigh dist
        double inc  = reb_random_rayleigh(5e-3);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, 1e-9, a, e, inc, Omega, apsis, phi);
        pt.r 		= 0.00000934532;
        reb_add(r, pt);
        r->usleep = 1000;
    }

    //Final setup parameters
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
    //Integrate!
    reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    //output energy
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("dE = %e",relE);
    }
}
