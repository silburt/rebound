//This is essentially the Kirsch example, using it as a vanilla example, need to vary various parameters and make sure the scaling relations of Armitage work.  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r, int index);
double draw_ainv_powerlaw(double min, double max);

double E0;
int N_prev;
char output_name[100] = {0};
char removed[200] = {0};
char* argv4;
double log_constant, tlog_output, lin_constant, tlin_output;
clock_t start_t;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    int seed = 10;
    srand(seed);
    strcat(output_name,"output/test");
    argv4 = "output/test";
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.hill_switch_factor = 6;         //Hill radii
    r->ri_hermes.adaptive_hill_switch_factor = 1;
    r->ri_hermes.solar_switch_factor = 20.;     //X*radius
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    double tmax = 1e5 * 6.283;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
//    r->boundary	= REB_BOUNDARY_OPEN;
//    const double boxsize = 25;
//    reb_configure_box(r,boxsize,2,2,1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
    
    //does order of adding bodies matter? Yes, apparently it does says Hanno. Jacobi coords.
    
    //planets
    double m_neptune = 5e-5, r_neptune = 1.6e-4;
    double m_earth = 3e-6, r_earth = 0.00004258689;
    {
        double a=1, m=m_neptune, e=0, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = r_neptune;
        p2.hash = r->N;
        reb_add(r, p2);
    }
    
    //planet
    {
        double a=2, m=m_neptune, e=0, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = r_neptune;
        p2.hash = r->N;
        reb_add(r, p2);
    }

    r->N_active = r->N;
    r->dt = pow(1,1.5)/50;
    
    //planetesimals-what's a reasonable ini? Perfectly cold disk seems unlikely...
    int N_planetesimals = 10;
    double planetesimal_mass = 1e-8;
    double amin = 0.85, amax=1.15;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,0);
        double e = reb_random_rayleigh(0.2);   //rayleigh dist
        double inc = reb_random_rayleigh(0.005);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e, inc, Omega, apsis, phi);
		pt.r 		= 0.00000934532;
        pt.hash = r->N;
		reb_add(r, pt);
    }
    
    int n_output = 5000;
    log_constant = pow(tmax + 1, 1./(n_output - 1));
    tlog_output = r->dt;

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    N_prev = r->N;
    
    //naming stuff
    char timeout[200] = {0};
    strcat(timeout,output_name);
    strcat(removed,output_name); strcat(removed,"_removed.txt");
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    system(syss);
    strcat(output_name,".txt");
    
    //Integrate!
    start_t = clock();
    reb_integrate(r, tmax);
    
    double time = (double)(clock() - start_t)/ CLOCKS_PER_SEC;
    strcat(timeout,"_elapsedtime.txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    //Center if ejected particle
    if(r->N < N_prev){
        N_prev = r->N;
        double E = reb_tools_energy(r);
        reb_move_to_com(r);
        r->energy_offset += E - reb_tools_energy(r);
    }
    
    if(r->t > tlog_output){//log output or linear output!!
        tlog_output = r->t*log_constant;
        
        double E = reb_tools_energy(r);
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        int N_mini = 0;
        if(r->integrator == REB_INTEGRATOR_HERMES) N_mini = r->ri_hermes.mini->N - r->ri_hermes.mini->N_active;
        
        double total_t = (double)(clock() - start_t) / CLOCKS_PER_SEC;
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%e,%e,%f,%f,%d,%d,%llu\n",r->t,dE,total_t,r->ri_hermes.hill_switch_factor,r->N,N_mini,r->ri_hermes.steps_miniactive);
        fclose(append);
        
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
    }
}

double calc_a(struct reb_simulation* r, int index){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[index]; //output planet only.
    const double mu = r->G*(com.m + p.m);
    const double dvx = p.vx-com.vx;
    const double dvy = p.vy-com.vy;
    const double dvz = p.vz-com.vz;
    const double dx = p.x-com.x;
    const double dy = p.y-com.y;
    const double dz = p.z-com.z;
    
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt( dx*dx + dy*dy + dz*dz );    //distance
    const double dinv = 1./d;
    const double a = -mu/( v2 - 2.*mu*dinv );
    
    return a;
}

//returns value randomly drawn from P(x) = 1/x distribution
double draw_ainv_powerlaw(double min, double max){
    double y = reb_random_uniform(0., 1.);
    return exp(y*log(max/min) + log(min));
}


