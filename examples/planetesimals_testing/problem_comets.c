/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double E0, LA0, LL0, t_output, t_log_output, xyz_t = 0;
int xyz_counter = 0, numdt = 20;
char* mercury_dir; char* swifter_dir;

double HSRbase;

//temp
int output_xyz = 0; //switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};
int *in_mini;
int N_CE = 0;
int L_CE = 0;
struct reb_particle comr0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    r->usleep = 10000;
    
    double tmax = INFINITY;
    int N_planetesimals = 100;
    int seed = 2;
    strcat(output_name,argv[1]); strcat(output_name,".txt"); argv4=argv[1];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->ri_hermes.solar_switch_factor = 50.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;        //Hill radii
    r->dt = 5;
    
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;     //switch to track the energy from collisions/ejections
    r->collision_resolve_keep_sorted = 1;
    
    //For many ejected planetesimals this leads to error jumps.
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 300;
    reb_configure_box(r,boxsize,1,1,1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(seed);
    t_log_output = 1.00048;
    t_output = r->dt;
    
    //planet 1
    {
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, 5e-5, 30, 0.01, 0, 0, 0, 0);
        p.r = 1.6e-4;              //radius of particle is in AU!
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    
    //planetesimals
    double planetesimal_mass = 1e-14;
    double amin = 30, amax = 40;        //for planetesimal disk
    double powerlaw = 0.5;
    double e = 0.9;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(1);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, e, inc, Omega, apsis,phi);
		pt.r 		= 1e-7;
        
        reb_add(r, pt);
    }
    
    //com
    reb_move_to_com(r);
    comr0 = reb_get_com(r);
    
    //energy & mom
    E0 = reb_tools_energy(r);
    system("rm -v output/energy.txt");
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    
    //in_mini
    in_mini = calloc(sizeof(int),r->N);
    
    //Integrate!
    reb_integrate(r, tmax);
    
}

void heartbeat(struct reb_simulation* r){
    //change hill switch radius
    
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        double E = reb_tools_energy(r);// + r->ri_hermes.com_dE;
        //double dE = fabs((E-E0)/E0);
        double dE = (E-E0)/E0;
        reb_output_timing(r, 0);
        printf("    dE=%e,N_mini=%d",dE,r->ri_hermes.mini->N);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;

        FILE *append;
        append = fopen("output/energy.txt", "a");
        fprintf(append, "%.16f,%.16f,%d,%d\n",r->t,dE,r->N,r->ri_hermes.mini->N);
        fclose(append);
        
    }
}
