/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
double E0;

double tout = 1;
void heartbeat(struct reb_simulation* r){
    if(r->t > tout){
        tout *= 1.01;
        //reb_integrator_synchronize(r);
        double dE = (reb_tools_energy(r) - E0)/E0;
        reb_output_timing(r, 0);
        printf("%e",dE);
        
        FILE* f = fopen("test.txt", "a");
        fprintf(f,"%f,%e\n",r->t,dE);
        fclose(f);
    }
}

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_create_simulation();
    r->integrator = REB_INTEGRATOR_WHFAST;
    //r->usleep = 1000;
	r->dt = 0.1;
	r->heartbeat = heartbeat;
	r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    r->ri_whfasthelio.safe_mode = 1; // Turn off safe mode. Need to call reb_integrator_synchronize() before outputs.

    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;     //switch to track the energy from collisions/ejections
    r->collision_resolve_keep_sorted = 1;
    
	struct reb_particle star = {0};
	star.m = 1.;
	reb_add(r, star);
	
    {
        double m=1e-4, a=1, e=0;
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, reb_random_uniform(0,2.*M_PI));
        p.r = 0.000467;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    {
        double m=1e-5, a=2, e=0.1;
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, reb_random_uniform(0,2.*M_PI));
        p.r = 0.000467;
        p.hash = r->N;
        reb_add(r, p);
    }

    r->N_active = r->N;
    
    int N_planetesimals = 0;
    double planetesimal_mass = 1e-8;
    //double amin = 0.4, amax = 0.6;        //for planetesimal disk
    double amin = 0.95, amax = 1.05;
    double powerlaw = 0;
    while(r->N<N_planetesimals + r->N_active){
        struct reb_particle pt = {0};
        double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.00001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
        pt.r 		= 4e-5;
        pt.hash = r->N;
        reb_add(r, pt);
    }
    
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
    system("rm -v test.txt");
	reb_integrate(r,1e6);
}

