//An example for massive body collisions. 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double E0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_IAS15;
    r->ri_hybarid.switch_ratio = 5;  //units of Hill radii
    r->ri_hybarid.CE_radius = 15.;   //X*radius
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    r->usleep = 50000;
    
    r->dt = 0.001;
    double tmax = INFINITY;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->collisions_track_dE = 1;
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    //planet 1
    {
        double e,a,m;
        a=0.5,m=5e-5,e=0;
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, 0);
        p.r = 1.6e-4;              //radius of particle is in AU!
        p.id = r->N;
        reb_add(r, p);
    }
    
    //massive body collision - com issues during collision I think, and energy scaling issues?
    {
        double e,a,m,f;
        a=0.55,m=1e-5;e=0.4,f=-1.252;    //for collisions with planetesimals
        //a=0.55,m=5e-5,e=0.4,f=-1.254;
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, f);
        p.r = 1.6e-4;              //radius of particle is in AU!
        p.id = r->N;
        reb_add(r, p);
    }
    
    {
        double e,a,m,f;
        a=0.7,m=1e-5;e=0,f=0;    //for collisions with planetesimals
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, f);
        p.r = 1.6e-4;              //radius of particle is in AU!
        p.id = r->N;
        reb_add(r, p);
    }
    
    r->N_active = r->N;

    {//planetesimal
        double e,a,m,f;
        a=0.71,m=1e-9;e=0,f=0.1;    //for collisions with planetesimals
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, 0, 0, 0, f);
        p.r = 1.6e-4;              //radius of particle is in AU!
        p.id = r->N;
        reb_add(r, p);
    }
    
    reb_move_to_com(r);
    
    system("rm -f energy.txt");
    
    //calculate initial energy
    E0 = reb_tools_energy(r);
    
    //Integrate!
    reb_integrate(r, tmax);
    
}

double tout = 0.1;
void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        tout *=1.01;
        double E = reb_tools_energy(r) + r->collisions_dE;
        double dE = fabs((E-E0)/E0);
        FILE* f = fopen("energy.txt","a+");
        int N_mini = 0;
        if (r->ri_hybarid.mini_active){
            N_mini = r->ri_hybarid.mini->N;
        }
        fprintf(f,"%e,%e,%d,%e,%d,%d,0\n",r->t,dE,N_mini,r->collisions_dE,r->N,N_mini);
        fclose(f);
    }
    
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r) + r->collisions_dE;
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        int N_mini = 0;
        if (r->ri_hybarid.mini_active){
            N_mini = r->ri_hybarid.mini->N;
        }
        printf("dE=%e,N_mini=%d",dE,N_mini);
    }
}
