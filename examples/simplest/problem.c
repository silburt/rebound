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

void heartbeat(struct reb_simulation* r){
    double dE = (reb_tools_energy(r) - E0)/E0;
	printf("%f, %e\n",r->t,dE);
}

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_create_simulation();
    r->integrator = REB_INTEGRATOR_WHFASTHELIO;
    r->usleep = 1000;
	r->dt = 0.1;
	r->heartbeat = heartbeat;
	r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.

	struct reb_particle p1 = {0};
	p1.m = 1.;
	reb_add(r, p1);
	
	struct reb_particle p2 = {0};
    p2.m = 1e-3;
	p2.x = 1;
	p2.vy = 1;
	reb_add(r, p2);

    r->N_active = r->N;
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
	reb_integrate(r,1e6);
}

